import click
import spglib
import phonopy
from phonopy.interface.vasp import (read_vasp_from_strings,
                                    get_vasp_structure_lines)
from phonopy.structure.atoms import PhonopyAtoms
from aiida.manage.configuration import load_profile
from aiida.orm import QueryBuilder, Int, Float, Bool, Str, load_node
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import run, submit
from aiida_phonopy.common.utils import phonopy_atoms_to_structure
load_profile()


def launch_aiida():

    Dict = DataFactory('dict')
    unitcell_str = """ Sr Ti O
   1.0
     3.9050000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    3.9050000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000    3.9050000000000000
 Sr Ti O
   1   1   3
Direct
   0.0000000000000000  0.0000000000000000  0.0000000000000000
   0.5000000000000000  0.5000000000000000  0.5000000000000000
   0.5000000000000000  0.0000000000000000  0.5000000000000000
   0.5000000000000000  0.5000000000000000  0.0000000000000000
   0.0000000000000000  0.5000000000000000  0.5000000000000000"""

    lat, pos, num = spglib.refine_cell(
        read_vasp_from_strings(unitcell_str).to_tuple())
    cell = read_vasp_from_strings(unitcell_str)
    cell = PhonopyAtoms(cell=lat, scaled_positions=pos, numbers=num)
    cell = read_vasp_from_strings(
        '\n'.join(get_vasp_structure_lines(cell)))
    structure = phonopy_atoms_to_structure(cell)

    base_incar_dict = {
        'PREC': 'Accurate',
        'IBRION': -1,
        'EDIFF': 1e-8,
        'NELMIN': 5,
        'NELM': 100,
        'ENCUT': 420,
        'IALGO': 38,
        'ISMEAR': 0,
        'SIGMA': 0.01,
        'GGA': 'PS',
        'LREAL': False,
        'lcharg': False,
        'lwave': False,
    }

    base_config = {'code_string': 'vasp544mpi@stern',
                   'potential_family': 'PBE.54',
                   'potential_mapping': {'O': 'O',
                                         'Ti': 'Ti_pv',
                                         'Sr': 'Sr_sv'},
                   'options': {'resources': {'parallel_env': 'mpi*',
                                             'tot_num_mpiprocs': 16},
                               'max_wallclock_seconds': 3600 * 10}}
    base_parser_settings = {'add_energies': True,
                            'add_forces': True,
                            'add_stress': True}
    forces_config = base_config.copy()
    forces_config.update({'kpoints_mesh': [4, 4, 4],  # k-point density,
                          'parser_settings': base_parser_settings,
                          'parameters': base_incar_dict})
    nac_config = base_config.copy()
    nac_parser_settings = {'add_born_charges': True,
                           'add_dielectrics': True}
    nac_parser_settings.update(base_parser_settings)
    nac_incar_dict = {'lepsilon': True}
    nac_incar_dict.update(base_incar_dict)
    nac_config.update({'kpoints_mesh': [8, 8, 8],  # k-point density,
                       'parser_settings': nac_parser_settings,
                       'parameters': nac_incar_dict})

    PhonopyIterHA = WorkflowFactory('phonopy.iter_ha')
    builder = PhonopyIterHA.get_builder()
    builder.structure = structure
    builder.calculator_settings = Dict(dict={'forces': forces_config,
                                             'nac': nac_config})
    builder.run_phonopy = Bool(False)
    builder.remote_phonopy = Bool(True)
    builder.code_string = Str('phonopy@stern')
    builder.phonon_settings = Dict(
        dict={'mesh': 50.0,
              'supercell_matrix': [2, 2, 2],
              'distance': 0.01,
              'is_nac': True,
              'fc_calculator': 'alm'})
    builder.symmetry_tolerance = Float(1e-5)
    builder.options = Dict(dict=base_config['options'])
    builder.metadata.label = "SrTiO3 iterative phonon 2x2x2"
    builder.metadata.description = "SrTiO3 iterative phonon 2x2x2"
    builder.max_iteration = Int(20)
    builder.number_of_snapshots = Int(50)
    builder.temperature = Float(300.0)

    # Chose how to run the calculation
    run_by_deamon = True
    if not run_by_deamon:
        result = run(builder)
        print(result)
    else:
        future = submit(builder)
        print(future)
        print('Running workchain with pk={}'.format(future.pk))


def dump_phonopy(pk, pk_nac):
    ph = get_phonon(pk, pk_nac)
    settings = {'force_sets': True,
                'displacements': True,
                'force_constants': False,
                'born_effective_charge': True,
                'dielectric_constant': True}
    filename = "phonopy_params-%d.yaml" % pk
    ph.save(filename=filename,
            settings=settings)
    print("%s was made for PK=%d." % (filename, pk))


def get_phonon(pk, pk_nac):
    n = load_node(pk)
    unitcell = n.inputs.structure.get_ase()
    smat = n.outputs.phonon_setting_info['supercell_matrix']
    ph = phonopy.load(unitcell=unitcell,
                      supercell_matrix=smat,
                      primitive_matrix='auto')
    force_sets = n.outputs.force_sets.get_array('force_sets')
    dataset = n.outputs.phonon_setting_info['displacement_dataset']
    ph.dataset = dataset
    ph.forces = force_sets

    n_nac = load_node(pk_nac)
    if 'nac_params' in n_nac.outputs:
        borns = n_nac.outputs.nac_params.get_array('born_charges')
        epsilon = n_nac.outputs.nac_params.get_array('epsilon')
        nac_params = {'born': borns,
                      'factor': 14.399652,
                      'dielectric': epsilon}
        ph.nac_params = nac_params

    return ph


def search_pk():
    qb = QueryBuilder()
    IterHarmonicApprox = WorkflowFactory('phonopy.iter_ha')
    qb.append(IterHarmonicApprox)
    qb.order_by({IterHarmonicApprox: {'ctime': 'desc'}})
    qb.first()
    uuid = qb.first()[0].uuid

    qb = QueryBuilder()
    qb.append(
        IterHarmonicApprox,
        tag='iter_ha',
        filters={'uuid': {'==': uuid}})

    PhonopyWorkChain = WorkflowFactory('phonopy.phonopy')
    qb.append(PhonopyWorkChain, with_incoming='iter_ha')
    qb.order_by({PhonopyWorkChain: {'ctime': 'asc'}})

    pks = [n[0].pk for n in qb.all() if n[0].is_finished_ok]

    return pks


def band_phonopy(pk, pk_nac):
    ph = get_phonon(pk, pk_nac)
    ph.produce_force_constants(fc_calculator='alm')
    filename = "band-%d.yaml" % pk
    ph.auto_band_structure(write_yaml=True, filename=filename)
    print("%s was made for PK=%d." % (filename, pk))


@click.command()
@click.argument('pks', nargs=-1)
@click.option('--dump', default=False, is_flag=True)
@click.option('--band', default=False, is_flag=True)
@click.option('--list-pks', 'list_pks', default=False, is_flag=True)
@click.option('--pk-nac', 'pk_nac', default=18077, type=int)
def main(pks, dump, band, list_pks, pk_nac):
    """

    pk_nac is the node that contains NAC params.

    Usage
    -----
    If no option is suppied, IterHarmonicApprox calculation starts.

    The command options are used to analyize the (on-going) results.

    To list PKs of PhonopyWorkChain calculations
    % python launch_phonon_SrTiO3.py --show-pks

    To write band-28158.yaml
    % python launch_phonon_SrTiO3.py --band 28158 ...

    To write all band-*.yaml
    % python launch_phonon_SrTiO3.py --band

    band-*.yaml can be visualized by phonopy-band command,
    % phonopy-bandplot --legacy --legend band-*.yaml
    Legacy mode (--legacy) is necessary to plot multiple band-*.yaml
    simultaneously.

    To write phonopy_params-28158.yaml ...
    % python launch_phonon_SrTiO3.py --dump 28158 ...

    To write all phonopy_params-*.yaml
    % python launch_phonon_SrTiO3.py --dump

    """

    if list_pks:  # List PKs of PhonopyWorkChain
        print(",".join(["%d" % pk for pk in search_pk()]))
    elif dump:  # write phonopy_params.yaml
        if pks:
            for pk in pks:
                dump_phonopy(int(pk), pk_nac)
        else:
            for pk in search_pk():
                dump_phonopy(int(pk), pk_nac)
    elif band:  # write band.yaml
        if pks:
            for pk in pks:
                band_phonopy(int(pk), pk_nac)
        else:
            for pk in search_pk():
                band_phonopy(int(pk), pk_nac)
    else:
        launch_aiida()


if __name__ == '__main__':
    main()

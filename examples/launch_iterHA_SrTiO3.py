import click
import spglib
import phonopy
import numpy as np
from phonopy.interface.vasp import (read_vasp_from_strings,
                                    get_vasp_structure_lines)
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.harmonic.displacement import get_displacements_and_forces
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
    builder.metadata.label = "SrTiO3 iterative phonon 2x2x2 1000K"
    builder.metadata.description = "SrTiO3 iterative phonon 2x2x2 1000K"
    builder.max_iteration = Int(40)
    builder.number_of_snapshots = Int(50)
    builder.temperature = Float(1000.0)
    builder.number_of_steps_for_fitting = Int(8)

    # Chose how to run the calculation
    run_by_deamon = True
    if not run_by_deamon:
        result = run(builder)
        print(result)
    else:
        future = submit(builder)
        print(future)
        print('Running workchain with pk={}'.format(future.pk))


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


def find_latest_uuid():
    IterHarmonicApprox = WorkflowFactory('phonopy.iter_ha')
    qb = QueryBuilder()
    qb.append(IterHarmonicApprox)
    qb.order_by({IterHarmonicApprox: {'ctime': 'desc'}})
    qb.first()
    return qb.first()[0].uuid


def search_pk(uuid):
    IterHarmonicApprox = WorkflowFactory('phonopy.iter_ha')

    if uuid is None:
        _uuid = find_latest_uuid()
    else:
        _uuid = uuid

    qb = QueryBuilder()
    qb.append(
        IterHarmonicApprox,
        tag='iter_ha',
        filters={'uuid': {'==': _uuid}})

    PhonopyWorkChain = WorkflowFactory('phonopy.phonopy')
    qb.append(PhonopyWorkChain, with_incoming='iter_ha')
    qb.order_by({PhonopyWorkChain: {'ctime': 'asc'}})

    pks = [n[0].pk for n in qb.all() if n[0].is_finished_ok]

    return pks


def get_num_prev(uuid):
    if uuid is None:
        _uuid = find_latest_uuid()
    else:
        _uuid = uuid
    n = load_node(_uuid)
    return n.inputs.number_of_steps_for_fitting.value


def bunch_phonons(pks, pk_nac):
    phonons = [get_phonon(pk, pk_nac) for pk in pks]
    all_forces = []
    all_disps = []
    for ph in phonons:
        disps, forces = get_displacements_and_forces(ph.dataset)
        all_disps.append(disps)
        all_forces.append(forces)
    d = np.concatenate(all_disps, axis=0)
    f = np.concatenate(all_forces, axis=0)
    ph = phonons[0]
    ph.dataset = {'displacements': d, 'forces': f}
    return ph


def bunch_band_phonopy(pks, pk_nac):
    ph = bunch_phonons(pks, pk_nac)
    ph.produce_force_constants(fc_calculator='alm')
    pks_str = get_pks_str(pks)
    filename = "bunch-band-%s.yaml" % pks_str
    ph.auto_band_structure(write_yaml=True, filename=filename)
    print("%s was made for PKs=%s." % (filename, pks_str))


def band_phonopy(pk, pk_nac):
    ph = get_phonon(pk, pk_nac)
    ph.produce_force_constants(fc_calculator='alm')
    filename = "band-%d.yaml" % pk
    ph.auto_band_structure(write_yaml=True, filename=filename)
    print("%s was made for PK=%d." % (filename, pk))


def get_bunch_tp_phonopy(pks, pk_nac, temperature):
    ph = bunch_phonons(pks, pk_nac)
    ph.produce_force_constants(fc_calculator='alm')
    ph.run_mesh(mesh=100.0, shift=[0.5, 0.5, 0.5])
    ph.run_thermal_properties(temperatures=[temperature, ])
    free_energy = ph.get_thermal_properties_dict()['free_energy'][0]
    return free_energy


def get_tp_phonopy(pk, pk_nac, temperature):
    ph = get_phonon(pk, pk_nac)
    ph.produce_force_constants(fc_calculator='alm')
    ph.run_mesh(mesh=100.0, shift=[0.5, 0.5, 0.5])
    ph.run_thermal_properties(temperatures=[temperature, ])
    free_energy = ph.get_thermal_properties_dict()['free_energy'][0]
    return free_energy


def bunch_dump_phonopy(pks, pk_nac):
    ph = bunch_phonons(pks, pk_nac)
    settings = {'force_sets': True,
                'displacements': True,
                'force_constants': False,
                'born_effective_charge': True,
                'dielectric_constant': True}
    pks_str = get_pks_str(pks)
    filename = "phonopy_params-%s.yaml" % pks_str
    ph.save(filename=filename, settings=settings)
    print("%s was made for PKs=%s." % (filename, pks_str))


def dump_phonopy(pk, pk_nac):
    ph = get_phonon(pk, pk_nac)
    settings = {'force_sets': True,
                'displacements': True,
                'force_constants': False,
                'born_effective_charge': True,
                'dielectric_constant': True}
    filename = "phonopy_params-%d.yaml" % pk
    ph.save(filename=filename, settings=settings)
    print("%s was made for PK=%d." % (filename, pk))


def get_pks_str(pks):
    if len(pks) > 4:
        pks_str = "%d-%dpks-%d" % (pks[0], len(pks) - 2, pks[-1])
    else:
        pks_str = "-".join(["%d" % pk for pk in pks])
    return pks_str


def get_pks_list(pks, num_prev):
    pks_list = []
    if len(pks) >= num_prev:
        for i in range(num_prev - 1):
            pks_list.append(pks[:(i + 1)])
        for i in range(len(pks) - num_prev + 1):
            pks_list.append(pks[i:(i + num_prev)])
    else:
        for i in range(len(pks)):
            pks_list.append(pks[:(i + 1)])
    return pks_list


def get_temperature(uuid):
    if uuid is None:
        _uuid = find_latest_uuid()
    else:
        _uuid = uuid
    n = load_node(_uuid)
    return n.inputs.temperature.value


@click.command()
@click.argument('pks', nargs=-1)
@click.option('--uuid', default=None)
@click.option('--dump', default=False, is_flag=True)
@click.option('--bunch-dump', 'bunch_dump', default=False, is_flag=True)
@click.option('--band', default=False, is_flag=True)
@click.option('--bunch-band', 'bunch_band', default=False, is_flag=True)
@click.option('--fe', 'free_energy', default=False, is_flag=True)
@click.option('--bunch-fe', 'bunch_free_energy', default=False, is_flag=True)
@click.option('--list-pks', 'list_pks', default=False, is_flag=True)
@click.option('--show-uuid', 'show_uuid', default=False, is_flag=True)
@click.option('--pk-nac', 'pk_nac', default=None, type=int)
@click.option('--num-prev', 'num_prev', default=None, type=int)
def main(uuid, pks, dump, bunch_dump, band, bunch_band,
         free_energy, bunch_free_energy,
         list_pks, show_uuid, pk_nac, num_prev):
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

    if pk_nac is None:
        _pk_nac = search_pk(uuid)[0]
    else:
        _pk_nac = pk_nac

    if num_prev is None:
        _num_prev = get_num_prev(uuid)
    else:
        _num_prev = num_prev

    if pks:
        _pks = [int(pk) for pk in pks]
    else:
        _pks = search_pk(uuid)

    if free_energy or bunch_free_energy:
        temperature = get_temperature(uuid)

    if show_uuid:  # Show latest IterHarmonicApprox node uuid
        print(find_latest_uuid())
    elif list_pks:  # List PKs of PhonopyWorkChain
        print(",".join(["%d" % pk for pk in search_pk(uuid)]))
        if len(_pks) > 0:
            _pks = _pks[1:]  # No.0 is omitted.
            print(get_pks_list(_pks, get_num_prev(uuid)))
    elif dump:  # write phonopy_params.yaml
        for pk in _pks:
            dump_phonopy(int(pk), _pk_nac)
    elif bunch_dump:
        if not pks:
            if len(_pks) > 0:
                _pks = _pks[1:]  # No.0 is omitted.
        if _pks:
            for pk_set in get_pks_list(_pks, _num_prev):
                bunch_dump_phonopy(pk_set, _pk_nac)
    elif band:  # write band.yaml
        for pk in _pks:
            band_phonopy(pk, _pk_nac)
    elif bunch_band:
        if not pks:
            if len(_pks) > 0:
                _pks = _pks[1:]  # No.0 is omitted.
        if _pks:
            for pk_set in get_pks_list(_pks, _num_prev):
                bunch_band_phonopy(pk_set, _pk_nac)
    elif free_energy:  # Show free energy
        for i, pk in enumerate(_pks):
            fe = get_tp_phonopy(pk, _pk_nac, temperature)
            print("%d %f" % (i, fe))

    elif bunch_free_energy:
        if not pks:
            if len(_pks) > 0:
                _pks = _pks[1:]  # No.0 is omitted.
        if _pks:
            for i, pk_set in enumerate(get_pks_list(_pks, _num_prev)):
                fe = get_bunch_tp_phonopy(pk_set, _pk_nac, temperature)
                print("%d %f" % (i + 1, fe))
    else:
        launch_aiida()


if __name__ == '__main__':
    main()

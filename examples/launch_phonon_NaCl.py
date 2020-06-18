from phonopy.interface.vasp import read_vasp_from_strings
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import submit
from aiida_phonopy.common.utils import phonopy_atoms_to_structure
from aiida.orm import Float, Bool, Str
from aiida import load_profile
load_profile()


def launch_aiida():

    Dict = DataFactory('dict')
    unitcell_str = """ Na Cl
   1.00000000000000
     5.6903014761756712    0.0000000000000000    0.0000000000000000
     0.0000000000000000    5.6903014761756712    0.0000000000000000
     0.0000000000000000    0.0000000000000000    5.6903014761756712
   4   4
Direct
  0.0000000000000000  0.0000000000000000  0.0000000000000000
  0.0000000000000000  0.5000000000000000  0.5000000000000000
  0.5000000000000000  0.0000000000000000  0.5000000000000000
  0.5000000000000000  0.5000000000000000  0.0000000000000000
  0.5000000000000000  0.5000000000000000  0.5000000000000000
  0.5000000000000000  0.0000000000000000  0.0000000000000000
  0.0000000000000000  0.5000000000000000  0.0000000000000000
  0.0000000000000000  0.0000000000000000  0.5000000000000000"""

    cell = read_vasp_from_strings(unitcell_str)
    structure = phonopy_atoms_to_structure(cell)

    base_incar_dict = {
        'PREC': 'Accurate',
        'IBRION': -1,
        'EDIFF': 1e-8,
        'NELMIN': 5,
        'NELM': 100,
        'ENCUT': 520,
        'IALGO': 38,
        'ISMEAR': 0,
        'SIGMA': 0.01,
        'LREAL': False,
        'lcharg': False,
        'lwave': False,
    }

    base_config = {'code_string': 'vasp544mpi@nancy',
                   'kpoints_density': 0.5,  # k-point density,
                   'potential_family': 'PBE.54',
                   'potential_mapping': {'Na': 'Na_pv', 'Cl': 'Cl'},
                   'options': {'resources': {'num_machines': 1,
                                             'parallel_env': 'mpi*',
                                             'tot_num_mpiprocs': 24},
                               'max_wallclock_seconds': 3600 * 10}}
    base_parser_settings = {'add_energies': True,
                            'add_forces': True,
                            'add_stress': True}
    forces_config = base_config.copy()
    forces_config.update({'parser_settings': base_parser_settings,
                          'parameters': base_incar_dict})
    nac_config = base_config.copy()
    nac_parser_settings = {'add_born_charges': True,
                           'add_dielectrics': True}
    nac_parser_settings.update(base_parser_settings)
    nac_incar_dict = {'lepsilon': True}
    nac_incar_dict.update(base_incar_dict)
    nac_config.update({'parser_settings': nac_parser_settings,
                       'parameters': nac_incar_dict})

    PhononPhonopy = WorkflowFactory('phonopy.phonopy')
    builder = PhononPhonopy.get_builder()
    builder.structure = structure
    builder.calculator_settings = Dict(dict={'forces': forces_config,
                                             'nac': nac_config})
    builder.run_phonopy = Bool(True)
    builder.remote_phonopy = Bool(True)
    builder.code_string = Str('phonopy@nancy')
    builder.phonon_settings = Dict(
        dict={'mesh': 50.0,
              'supercell_matrix': [2, 2, 2],
              'distance': 0.01,
              'is_nac': True})
    builder.symmetry_tolerance = Float(1e-5)
    builder.options = Dict(dict=base_config['options'])
    builder.metadata.label = "NaCl 2x2x2 phonon example"
    builder.metadata.description = "NaCl 2x2x2 phonon example"

    future = submit(builder)
    print(future)
    print('Running workchain with pk={}'.format(future.pk))


if __name__ == '__main__':
    launch_aiida()

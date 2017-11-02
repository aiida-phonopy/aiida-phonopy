from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.orm import CalculationFactory, DataFactory
from aiida.work.run import run, submit, async
from aiida.orm.data.structure import StructureData
from aiida.orm.data.base import Str, Float, Bool

VaspCalculation = CalculationFactory('vasp.vasp')
PwCalculation = CalculationFactory('quantumespresso.pw')
PhonopyCalculation = CalculationFactory('phonopy')

KpointsData = DataFactory("array.kpoints")
ParameterData = DataFactory('parameter')


# Define structure

import numpy as np

cell = [[ 3.1900000572, 0,           0],
        [-1.5950000286, 2.762621076, 0],
        [ 0.0,          0,           5.1890001297]]

structure = StructureData(cell=cell)

scaled_positions=[(0.6666669,  0.3333334,  0.0000000),
                  (0.3333331,  0.6666663,  0.5000000),
                  (0.6666669,  0.3333334,  0.3750000),
                  (0.3333331,  0.6666663,  0.8750000)]

symbols=['Ga', 'Ga', 'N', 'N']

positions = np.dot(scaled_positions, cell)

for i, scaled_position in enumerate(scaled_positions):
    structure.append_atom(position=np.dot(scaled_position, cell).tolist(),
                          symbols=symbols[i])


# PHONOPY settings
ph_settings = ParameterData(dict={'supercell': [[2, 0, 0],
                                                [0, 2, 0],
                                                [0, 0, 2]],
                                  'primitive': [[1.0, 0.0, 0.0],
                                                [0.0, 1.0, 0.0],
                                                [0.0, 0.0, 1.0]],
                                  'distance': 0.01,
                                  'mesh': [20, 20, 20],
                                  'symmetry_precision': 1e-5
                                  # 'code': 'phonopy@boston'  # comment to use local phonopy
                                  })

# VASP SPECIFIC
if True:   # Set TRUE to use VASP or FALSE to use Quantum Espresso
    incar_dict = {
        'NELMIN' : 5,
        'NELM'   : 100,
        'ENCUT'  : 400,
        'ALGO'   : 38,
        'ISMEAR' : 0,
        'SIGMA'  : 0.01,
        'GGA'    : 'PS'
    }

    es_settings = ParameterData(dict=incar_dict)

    from pymatgen.io import vasp as vaspio

    potcar = vaspio.Potcar(symbols=['Ga', 'N'],
                           functional='PBE')

    # custom k-points
    # supported_modes = Enum(("Gamma", "Monkhorst", "Automatic", "Line_mode", "Cartesian", "Reciprocal"))
    kpoints_dict = {'type': 'Monkhorst',
                    'points': [2, 2, 2],
                    'shift': [0.0, 0.0, 0.0]}

    settings_dict = {'code': {'optimize': 'vasp544mpi@boston',
                              'forces': 'vasp544mpi@boston',
                              'born_charges': 'vasp544mpi@boston'},
                     'parameters': incar_dict,
                     #'kpoints': kpoints_dict,
                     'kpoints_per_atom': 100,  # k-point density
                     'pseudos': potcar.as_dict()}

    # pseudos = ParameterData(dict=potcar.as_dict())
    es_settings = ParameterData(dict=settings_dict)


# QE SPECIFIC
if False:
    parameters_dict = {
        'CONTROL': {'calculation': 'scf',
                    'tstress': True,  #  Important that this stays to get stress
                    'tprnfor': True,},
        'SYSTEM': {'ecutwfc': 30.,
                   'ecutrho': 200.,},
        'ELECTRONS': {'conv_thr': 1.e-6,}
    }

    # Kpoints
    #kpoints_mesh = 2
    #kpoints = KpointsData()
    #kpoints.set_kpoints_mesh([kpoints_mesh, kpoints_mesh, kpoints_mesh])
    #code = Code.get_from_string('pw@stern_outside')

    pseudos = Str('pbe_ps')

    settings_dict = {'code': 'pw@stern_outside',
                     'parameters': parameters_dict,
                     'kpoints_per_atom': 1000,  # k-point density
                     'pseudos_family': 'pbe_ps'}

    es_settings = ParameterData(dict=settings_dict)


# LAMMPS SPECIFIC
if False:
    # GaN Tersoff potentials parameters (can be readed from a file if needed)
    tersoff_gan = {
        'Ga Ga Ga': '1.0 0.007874 1.846 1.918000 0.75000 -0.301300 1.0 1.0 1.44970 410.132 2.87 0.15 1.60916 535.199',
        'N  N  N': '1.0 0.766120 0.000 0.178493 0.20172 -0.045238 1.0 1.0 2.38426 423.769 2.20 0.20 3.55779 1044.77',
        'Ga Ga N': '1.0 0.001632 0.000 65.20700 2.82100 -0.518000 1.0 0.0 0.00000 0.00000 2.90 0.20 0.00000 0.00000',
        'Ga N  N': '1.0 0.001632 0.000 65.20700 2.82100 -0.518000 1.0 1.0 2.63906 3864.27 2.90 0.20 2.93516 6136.44',
        'N  Ga Ga': '1.0 0.001632 0.000 65.20700 2.82100 -0.518000 1.0 1.0 2.63906 3864.27 2.90 0.20 2.93516 6136.44',
        'N  Ga N ': '1.0 0.766120 0.000 0.178493 0.20172 -0.045238 1.0 0.0 0.00000 0.00000 2.20 0.20 0.00000 0.00000',
        'N  N  Ga': '1.0 0.001632 0.000 65.20700 2.82100 -0.518000 1.0 0.0 0.00000 0.00000 2.90 0.20 0.00000 0.00000',
        'Ga N  Ga': '1.0 0.007874 1.846 1.918000 0.75000 -0.301300 1.0 0.0 0.00000 0.00000 2.87 0.15 0.00000 0.00000'}

    potential = {'pair_style': 'tersoff',  # potential type
                 'data': tersoff_gan}  # potential data

    parameters = {'relaxation': 'tri',  # iso/aniso/tri
                  'pressure': 0.0,  # kbars
                  'vmax': 0.000001,  # Angstrom^3
                  'energy_tolerance': 1.0e-25,  # eV
                  'force_tolerance': 1.0e-25,  # eV angstrom
                  'max_evaluations': 1000000,
                  'max_iterations': 500000}

    settings_dict = {'code': {'optimize': 'lammps_optimize@boston',
                              'forces': 'lammps_force@boston'},
                     'parameters': parameters,
                     'potential': potential}

    es_settings = ParameterData(dict=settings_dict)


# CODE INDEPENDENT
machine_dict = {'resources': {'num_machines': 1,
                              'parallel_env': 'mpi*',
                              'tot_num_mpiprocs': 16},
                'max_wallclock_seconds': 30 * 60,
                }

machine = ParameterData(dict=machine_dict)

from aiida.workflows.wc_phonon import PhononPhonopy

if True:
    result = run(PhononPhonopy,
                 structure=structure,
                 machine=machine,
                 es_settings=es_settings,
                 ph_settings=ph_settings,
                 # Optional settings
                 pressure=Float(0.0),
                 optimize=Bool(True),
                 )

    print (result)
else:
    future = submit(PhononPhonopy,
                    structure=structure,
                    machine=machine,
                    es_settings=es_settings,
                    ph_settings=ph_settings,
                    # Optional settings
                    pressure=Float(0),
                    optimize=Bool(True),
                    )

    print future
    print('Running workchain with pk={}'.format(future.pid))
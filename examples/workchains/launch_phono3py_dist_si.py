# Works run by the daemon (using submit)

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()


from aiida.orm import load_node, DataFactory, WorkflowFactory
from aiida.work.run import run, submit, async

from aiida.orm.data.base import Str, Float, Bool, Int

import numpy as np

ForceConstantsData = DataFactory('phonopy.force_constants')
ForceSetsData = DataFactory('phonopy.force_sets')
PhononDosData = DataFactory('phonopy.phonon_dos')
NacData = DataFactory('phonopy.nac')

ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')


# Silicon structure
a = 5.404
cell = [[a, 0, 0],
        [0, a, 0],
        [0, 0, a]]

symbols=['Si'] * 8
scaled_positions = [(0.875,  0.875,  0.875),
                    (0.875,  0.375,  0.375),
                    (0.375,  0.875,  0.375),
                    (0.375,  0.375,  0.875),
                    (0.125,  0.125,  0.125),
                    (0.125,  0.625,  0.625),
                    (0.625,  0.125,  0.625),
                    (0.625,  0.625,  0.125)]

structure = StructureData(cell=cell)
for i, scaled_position in enumerate(scaled_positions):
    structure.append_atom(position=np.dot(scaled_position, cell).tolist(),
                          symbols=symbols[i])

# Machine
machine_dict = {'resources': {'num_machines': 1,
                              'parallel_env': 'mpi*',
                              'tot_num_mpiprocs': 16},
                'max_wallclock_seconds': 3600 * 10,
                }

machine = ParameterData(dict=machine_dict)

# PHONOPY settings
ph_settings = ParameterData(dict={'supercell': [[2, 0, 0],
                                                [0, 2, 0],
                                                [0, 0, 2]],
                                  'primitive': [[0.0, 0.5, 0.5],
                                                [0.5, 0.0, 0.5],
                                                [0.5, 0.5, 0.0]],
                                  'distance': 0.01,
                                  'mesh': [20, 20, 20],
                                  'symmetry_precision': 1e-5,
                                  'code': 'phono3py@stern_in',
                                  'machine': machine_dict})

Phono3pyDist = WorkflowFactory('phonopy.phono3py_dist')

# Chose how to run the calculation
run_by_deamon = False
if not run_by_deamon:

    result = run(Phono3pyDist,
                 structure=structure,
                 parameters=ph_settings,
                 data_sets=load_node(81481),  # load phonon3 WorkChain output data_set
                 gp_chunks=Int(8)
                 )

    print (result)
else:
    future = submit(Phono3pyDist,
                    structure=structure,
                    parameters=ph_settings,
                    data_sets=load_node(81481),  # load phonon3 WorkChain output data_set
                    gp_chunks=Int(8)
                    )
    print future
    print('Running workchain with pk={}'.format(future.pid))
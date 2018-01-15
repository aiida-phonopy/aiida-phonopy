# Calculate DOS Band structure and thermal properties
from aiida import load_dbenv
load_dbenv()

from aiida.orm import Code, DataFactory, load_node
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')

import numpy as np
import os

codename = 'phonopy_tot@boston_in'
code = Code.get_from_string(codename)


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
positions = np.dot(scaled_positions, cell)

for i, scaled_position in enumerate(scaled_positions):
    structure.append_atom(position=np.dot(scaled_position, cell).tolist(),
                          symbols=symbols[i])

# Machine
machine_dict = {'resources': {'num_machines': 1,
                              'parallel_env': 'mpi*',
                              'tot_num_mpiprocs': 16},
                'max_wallclock_seconds': 3600 * 10,
                }


# PHONOPY settings
parameters = ParameterData(dict={'supercell': [[2, 0, 0],
                                                [0, 2, 0],
                                                [0, 0, 2]],
                                  'primitive': [[0.0, 0.5, 0.5],
                                                [0.5, 0.0, 0.5],
                                                [0.5, 0.5, 0.0]],
                                  'distance': 0.01,
                                  'mesh': [20, 20, 20],
                                  'symmetry_precision': 1e-5,

                                  })

calc = code.new_calc(max_wallclock_seconds=3600,
                     resources={"num_machines": 1,
                                "parallel_env":"mpi*",
                                "tot_num_mpiprocs": 6})


calc.label = "test phonopy calculation"
calc.description = "A much longer description"

calc.use_structure(load_node(79024))
calc.use_code(code)
calc.use_parameters(parameters)

# Chose to use forces or force constants
if True:
    calc.use_data_sets(load_node(81106))  # This node should contain a ForceSetsData object
else:
    calc.use_force_constants(load_node(25806))  # This node should contain a ForceConstantsData object

# Set bands (optional)
from aiida_phonopy.workchains.phonon import get_primitive, get_path_using_seekpath
primitive = get_primitive(structure, parameters)['primitive_structure']
bands = get_path_using_seekpath(primitive, band_resolution=30)
calc.use_bands(bands)

if False:
    subfolder, script_filename = calc.submit_test()
    print "Test_submit for calculation (uuid='{}')".format(calc.uuid)
    print "Submit file in {}".format(os.path.join(
        os.path.relpath(subfolder.abspath),
        script_filename))
else:
    calc.store_all()
    print "created calculation with PK={}".format(calc.pk)
    calc.submit()
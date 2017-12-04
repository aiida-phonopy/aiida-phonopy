# Calculate DOS Band structure and thermal properties
# Additionally if DATA_SETS (forces) are used force constants are also calculated
from aiida import load_dbenv
load_dbenv()

from aiida.orm import Code, DataFactory, load_node
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
ForceSetsData = DataFactory('phonopy.force_sets')
ForceConstantsData = DataFactory('phonopy.force_constants')

import numpy as np
import os

codename = 'phonopy_tot@boston_in'
code = Code.get_from_string(codename)

cell = [[ 4.2119998932, 0.0,          0.0],
        [ 0.0,          4.2119998932, 0.0],
        [ 0.0,          0.0,          4.2119998932]]

structure = StructureData(cell=cell)

scaled_positions=[(0.0000000,  0.0000000,  0.0000000),
                  (0.0000000,  0.5000000,  0.5000000),
                  (0.5000000,  0.0000000,  0.5000000),
                  (0.5000000,  0.5000000,  0.0000000),
                  (0.5000000,  0.5000000,  0.5000000),
                  (0.5000000,  0.0000000,  0.0000000),
                  (0.0000000,  0.5000000,  0.0000000),
                  (0.0000000,  0.0000000,  0.5000000)]

symbols=['Mg', 'Mg', 'Mg', 'Mg', 'O', 'O', 'O', 'O']


for i, scaled_position in enumerate(scaled_positions):
    structure.append_atom(position=np.dot(scaled_position, cell).tolist(),
                  symbols=symbols[i])

parameters = ParameterData(dict={'supercell': [[2, 0, 0],
                                               [0, 2, 0],
                                               [0, 0, 2]],
                                 'primitive': [[1.0, 0.0, 0.0],
                                               [0.0, 1.0, 0.0],
                                               [0.0, 0.0, 1.0]],
                                 'distance': 0.01,
                                 'mesh': [40, 40, 40],
                                 'symmetry_precision': 1e-5}
                          )

calc = code.new_calc(max_wallclock_seconds=3600,
                     resources={"num_machines": 1,
                                "parallel_env":"mpi*",
                                "tot_num_mpiprocs": 6})


calc.label = "test phonopy calculation"
calc.description = "A much longer description"

calc.use_structure(structure)
calc.use_code(code)
calc.use_parameters(parameters)

# Chose to use forces or force constants
if True:
    # Use forces
    force_sets = ForceSetsData()
    force_sets.read_from_phonopy_file('FORCE_SETS')
    calc.use_data_sets(force_sets)
else:
    # Use force constants
    force_constants = ForceConstantsData()
    force_constants.read_from_phonopy_file('FORCE_CONSTANTS')
    calc.use_force_constants(force_constants)

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
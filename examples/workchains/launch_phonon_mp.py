# This script gets a crystal structure from Materials Project and performs a phonon calculation:
# To use this script you neef to define an environment variable with your Materials Project API key
# in your shell .profile/.bashrc as: export PMG_MAPI_KEY="Y0urK3yH3r3"

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.orm import CalculationFactory, DataFactory, WorkflowFactory
from aiida.work.run import run, submit
from aiida.orm.data.structure import StructureData
from aiida.orm.data.base import Str, Float, Bool, Int

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen, os

KpointsData = DataFactory("array.kpoints")
ParameterData = DataFactory('parameter')

import numpy as np

def get_supercell_size(structure, max_atoms=100, crystal_system=None):

    def axis_symmetry(axis, crystal_system):
        symmetry_dict = {'cubic':      [1, 1, 1],
                         'hexagonal':  [1, 1, 2],
                         'tetragonal': [1, 1, 2],
                         'monoclinic': [1, 1, 2],
                         'trigonal':   [1, 1, 2]}

        try:
            return np.where(np.array(symmetry_dict[crystal_system]) == symmetry_dict[crystal_system][axis])[0]
        except KeyError:
            # If symmetry not defined in symmetry_dict or is None
            return np.array([0, 1, 2])

    cell = np.array(structure.cell)
    num_atoms = len(structure.sites)

    supercell_size = [1, 1, 1]
    while True:

        test_cell = np.dot(cell.T, np.diag(supercell_size)).T
        norm = np.linalg.norm(test_cell, axis=1)
        index = np.argmin(norm)
        supercell_size_test = list(supercell_size)

        for i in axis_symmetry(index, crystal_system):
            supercell_size_test[i] += 1

        anum_atoms_supercell = num_atoms * np.prod(supercell_size_test)
        if anum_atoms_supercell > max_atoms:
            atoms_minus =  num_atoms * np.prod(supercell_size)
            atoms_plus = num_atoms * np.prod(supercell_size_test)

            if max_atoms - atoms_minus < atoms_plus - max_atoms:
                return supercell_size
            else:
                return supercell_size_test
        else:
            supercell_size = supercell_size_test

# Set Materials Project structure ID here:
structure_id = 'mp-1265'

rester = pymatgen.MPRester(os.environ['PMG_MAPI_KEY'])

pmg_structure = rester.get_structure_by_material_id(structure_id)
pmg_band = rester.get_bandstructure_by_material_id(structure_id)

material_name = pmg_structure.formula.replace('1','').replace(' ','')

spa = SpacegroupAnalyzer(pmg_structure)

conventional = spa.get_conventional_standard_structure()
primitive = spa.get_primitive_standard_structure()

print ('Conventional cell')
print conventional
primitive_matrix = np.dot(np.linalg.inv(conventional.lattice.matrix), primitive.lattice.matrix)
primitive_matrix = np.round(primitive_matrix, decimals=6).tolist()

structure = StructureData(pymatgen=conventional).store()

print ('Primitive matrix')
print primitive_matrix

crystal_system = spa.get_crystal_system()
print 'Crystal system: {}'.format(crystal_system)

supercell_size = get_supercell_size(structure, crystal_system=crystal_system)
supercell = np.diag(supercell_size).tolist()
print ('Supercell shape: {}'.format(supercell_size))


# Machine
machine_dict = {'resources': {'num_machines': 1,
                              'parallel_env': 'mpi*',
                              'tot_num_mpiprocs': 16},
                'max_wallclock_seconds': 3600 * 10,
                }

# PHONOPY settings
ph_settings = ParameterData(dict={'supercell': supercell,
                                  'primitive': primitive_matrix,
                                  'mesh': [80, 80, 80],
                                  'distance': 0.01,
                                  'symmetry_precision': 1e-5,
                                  # Uncomment the following lines to use phonopy remotely
                                  # 'code': 'phonopy_tot@boston_in', # this uses phonopy.force_constants plugin
                                  # 'machine': machine_dict
                                  })

# VASP SPECIFIC
incar_dict = {
    'NELMIN' : 5,
    'NELM'   : 100,
    'ENCUT'  : 400,
    'ALGO'   : 38,
    'ISMEAR' : 0,
    'SIGMA'  : 0.01,
    'GGA'    : 'PS'
}

settings_dict = {'code': 'vasp@stern_in',
                 'parameters': incar_dict,
                 'kpoints_density': 0.5,  # k-point density,
                 'pseudos_family': 'pbe_test_family',
                 'family_folder': '/Users/abel/VASP/test_paw/',
                 'machine': machine_dict
                 }

es_settings = ParameterData(dict=settings_dict)

PhononPhonopy = WorkflowFactory('phonopy.phonon')

# Chose how to run the calculation
run_by_deamon = False
if not run_by_deamon:
    result = run(PhononPhonopy,
                 structure=structure,
                 es_settings=es_settings,
                 ph_settings=ph_settings,
                 # Optional settings
                 # pressure=Float(0),
                 optimize=Bool(True),
                 use_nac=Bool(True)
                 )

    print (result)
else:
    future = submit(PhononPhonopy,
                    structure=structure,
                    es_settings=es_settings,
                    ph_settings=ph_settings,
                    # Optional settings
                    # pressure=Float(0),
                    # optimize=Bool(False),
                    use_nac=Bool(True)
                    )

    print future
    print('Running workchain with pk={}'.format(future.pid))
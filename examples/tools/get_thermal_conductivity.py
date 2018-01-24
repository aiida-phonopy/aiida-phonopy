########################################################
# This is an example of how to use the data from phonon3
# WorkChain to calculate the thermal conductivity locally
# with phono3py using the information obtained by phonon3
# Workchain. To run this script use a phono3 workchain
# pk number:
# $ python get_thermal_conductivity.py pknumber
#########################################################

from aiida import load_dbenv
load_dbenv()

from aiida.orm import load_node, load_workflow, DataFactory

from phonopy.harmonic.force_constants import show_drift_force_constants
from phono3py.phonon3.fc3 import show_drift_fc3
from phono3py.phonon3 import Phono3py

from aiida_phonopy.workchains.phonon import phonopy_bulk_from_structure

import sys


ParameterData = DataFactory('parameter')

if len(sys.argv) < 2:
    print ('use: python get_thermal_conductivity.py {pk_number}')
    exit()

# Set WorkChain phonon3 PK number
################################
wc = load_node(int(sys.argv[1]))
################################

# get phono3py data from database
force_sets = wc.out.force_sets
structure = wc.out.final_structure
ph_settings = wc.inp.ph_settings

# If modifications in phonon parameters are necessary
ph_settings_dict = ph_settings.get_dict()
ph_settings_dict['mesh'] = [40, 40, 40]
ph_settings = ParameterData(dict=ph_settings_dict)

# Create phono3py and calculate as usual
phono3py = Phono3py(phonopy_bulk_from_structure(structure),
                    supercell_matrix=ph_settings.dict.supercell,
                    primitive_matrix=ph_settings.dict.primitive,
                    symprec=ph_settings.dict.symmetry_precision,
                    mesh=ph_settings.dict.mesh,
                    log_level=1)

phono3py.produce_fc3(force_sets.get_forces3(),
                     displacement_dataset=force_sets.get_data_sets3(),
                     is_translational_symmetry=True,
                     is_permutation_symmetry=True,
                     is_permutation_symmetry_fc2=True)

fc3 = phono3py.get_fc3()
fc2 = phono3py.get_fc2()

if True:
    from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5
    print ('Writing FC2 and FC3 into HDF5 files')
    write_fc2_to_hdf5(fc3, filename='fc2.hdf5')
    write_fc3_to_hdf5(fc3, filename='fc3.hdf5')

show_drift_fc3(fc3)
show_drift_force_constants(fc2, name='fc2')

# Use NAC if available
use_nac = False
if 'nac_data' in wc.get_outputs():
    primitive = phono3py.get_phonon_primitive()
    nac_params = wc.out.nac_data.get_born_parameters_phonopy(primitive_cell=primitive.get_cell())
    phono3py.set_phph_interaction(nac_params=nac_params)


phono3py.run_thermal_conductivity(temperatures=range(0, 1001, 10),
                                  boundary_mfp=1e6,  # This is to avoid divergence of phonon life time.
                                  write_kappa=True)

# Conductivity_RTA object (https://git.io/vVRUW)
cond_rta = phono3py.get_thermal_conductivity()


########################################################
# This is an example of how to use the data from phonon3
# WorkChain to obtain the input files necessary for running
# a phono3py calculation "manually".
# To run this script use a phono3py workchain pk number:
# $ python get_phono3py_inputs.py pknumber
#########################################################
from aiida import load_dbenv
load_dbenv()

from aiida.plugins import load_node, DataFactory
from aiida_phonopy.common.utils import phonopy_atoms_from_structure
from aiida_phonopy.common.raw_parsers import get_poscar_txt, get_BORN_txt

from phono3py.file_IO import write_disp_fc3_yaml, write_FORCES_FC3
from phonopy.structure.cells import get_supercell

import sys


ParameterData = DataFactory('parameter')

if len(sys.argv) < 2:
    print ('use: python get_phonon3py_inputs.py {pk_number}')
    exit()

# Set WorkChain phonon3 PK number
################################
wc = load_node(int(sys.argv[1]))
################################

# get phono3py data from database
force_sets = wc.out.force_sets
structure = wc.out.final_structure
ph_settings = wc.inp.ph_settings

supercell = get_supercell(phonopy_atoms_from_structure(structure),
                          ph_settings.dict.supercell,
                          symprec=ph_settings.dict.symmetry_precision)

write_disp_fc3_yaml(force_sets.get_data_sets3(),
                    supercell,
                    filename='disp_fc3.yaml')


write_FORCES_FC3(force_sets.get_data_sets3(), force_sets.get_forces3(), filename='FORCES_FC3')

with open('POSCAR-unitcell', mode='w') as f:
    f.writelines(get_poscar_txt(structure))

if 'nac_data' in wc.get_outputs_dict():
    nac_data = wc.out.nac_data
    print ('Writing BORN file')
    with open('BORN', mode='w') as f:
        f.writelines(get_BORN_txt(nac_data,
                                  structure=structure,
                                  parameters=ph_settings))

if 'force_constants_2order' and 'force_constants_3order' in wc.get_outputs_dict():
    from phono3py.file_IO import write_fc2_to_hdf5, write_fc3_to_hdf5

    print ('Writing FC2 and FC3 into HDF5 files')
    fc2 = wc.out.force_constants_2order.get_data()
    write_fc2_to_hdf5(fc2, filename='fc2.hdf5')
    fc3 = wc.out.force_constants_3order.get_data()
    write_fc3_to_hdf5(fc3, filename='fc3.hdf5')

########################################################
# This is an example of how to use the data from phonon3
# WorkChain to obtain the input files necessary for running
# a phono3py calculation.
# To run this script use a phono3py workchain pk number:
# $ python get_phono3py_inputs.py pknumber
#########################################################
from aiida import load_dbenv
load_dbenv()

from aiida.orm import load_node, DataFactory

from phono3py.file_IO import write_disp_fc3_yaml, write_FORCES_FC3
from phonopy.structure.cells import get_supercell
from aiida_phonopy.workchains.phonon import phonopy_bulk_from_structure
from aiida_phonopy.common.raw_parsers import get_poscar_txt

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

supercell = get_supercell(phonopy_bulk_from_structure(structure),
                          ph_settings.dict.supercell,
                          symprec=ph_settings.dict.symmetry_precision)

write_disp_fc3_yaml(force_sets.get_data_sets3(),
                    supercell,
                    filename='disp_fc3.yaml')


write_FORCES_FC3(force_sets.get_data_sets3(), force_sets.get_forces3(), filename='FORCES_FC3')

with open('POSCAR-unitcell', mode='w') as f:
    f.writelines(get_poscar_txt(structure))

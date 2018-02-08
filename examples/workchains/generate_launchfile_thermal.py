# This file is used to generate input launch file for thermal conductivity calculations using VASP
# Y.-W. FANG, created Feb. 8, 2018
# in 3 functions: generate_vasp_params(), generate_qe_params() and generate_lammps_params().
# generate_inputs() function at the end of the file decides which function to use according to the plugin

#Import the modules needed to run this script
import os
from os.path import exists
import numpy as np

#define spaces: four spaces == indent_N (N=1,2,3..)
indent_1, indent_2, indent_3, indent_4 = (' '*4, ' '*8, ' '*12, ' '*16)

title = input("Enter the title for the output python script: ")
if title.endswith('.py'): #TRUE if input file has a file extension of .py
    pass
else:
    title = title + '.py'
#Remove spaces from the title.
title = title.replace(" ", "_")
print(title,'\n')

#Ask whether to overwirte it if file exists
if  exists(title):
    overwrite = input("The file has been already existed. Overwrite it? y=yes, n=no\n")
    if overwrite.lower() == 'n':
        print('Program exits because you donnot want to overwrite it')
        exit(1)
    else:
        pass

with open(title,'wt') as f:
    f.write('#############################################################\n')
    f.write('# This script is for generating the launch scirpt used for  #\n')
    f.write('# calculating the thermal conductivity using VASP.          #\n')
    f.write('#############################################################\n')
    f.write('\n'*2)
    f.write('from aiida import load_dbenv, is_dbenv_loaded\n')
    f.write('if not is_dbenv_loaded():\n')
    f.write(indent_1 + 'load_dbenv()\n')
    f.write('\n')
    f.write('from aiida.orm import CalculationFactory, DataFactory, WorkflowFactory, load_node\n')
    f.write('from aiida.work.run import run, submit\n')
    f.write('from aiida.orm.data.structure import StructureData\n')
    f.write('from aiida.orm.data.base import Str, Float, Bool, Int\n')
    f.write("KpointsData = DataFactory('array.kpoints')\n")
    f.write("ParameterData = DataFactory('parameter')\n")
    f.write("import numpy as np\n")
    f.write("# Define structure\n")

    with open('POSCAR', 'rt', encoding='utf-8') as f1:
        line_number=0
        for line in f1:
            line_number=line_number+1
            if line_number==2:
                a = float(line)
                print(a,type(a))
            if line_number==3:
                lattice_vector_1=line.split()
                lattice_vector_1=list(map(float, lattice_vector_1))
                lattice_vector_1=np.array(lattice_vector_1)
                print(lattice_vector_1[:],type(lattice_vector_1))


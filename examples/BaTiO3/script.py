# -*- coding: utf-8 -*-
"""Simple script for running a PhonopyCalculation."""
from aiida import load_profile, orm
from aiida.engine import submit
import phonopy

from aiida_phonopy.calculations.phonopy import PhonopyCalculation
from aiida_phonopy.data import PhonopyData, PreProcessData

load_profile()


def main():
    """Run a simple Phonopy post-processing using PhonopyCalculation."""
    ph = phonopy.load('./phonopy.yaml', force_sets_filename='./FORCE_SETS')

    preprocess_data = PreProcessData(phonopy_atoms=ph.unitcell, supercell_matrix=[3, 3, 3])
    preprocess_data.set_displacements_from_dataset(dataset=ph.dataset)

    phonopy_data = PhonopyData(preprocess_data=preprocess_data)
    phonopy_data.set_forces(sets_of_forces=ph.forces)

    code = orm.load_code('phonopy@localhost')
    parameters = {'IRREPS': [0, 0, 0], 'FC_SYMMETRY': True, 'WRITE_FORCE_CONSTANTS': True}

    inputs = {
        'code': code,
        'parameters': orm.Dict(parameters),
        'phonopy_data': phonopy_data,
        'metadata': {
            'options': {
                'resources': {
                    'num_machines': 1
                },
            }
        }
    }

    submit(PhonopyCalculation, **inputs)


if __name__ == '__main__':
    main()

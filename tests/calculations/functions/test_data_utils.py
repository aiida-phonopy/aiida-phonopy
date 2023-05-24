# -*- coding: utf-8 -*-
"""Tests for :mod:`calculations.functions.data_utils`."""
from aiida import orm


def test_generate_preprocess_data(generate_structure):
    """Test the `generate_preprocess_data` calcfunction."""
    from aiida_phonopy.calculations.functions.data_utils import generate_preprocess_data

    inputs = {
        'structure': generate_structure(),
        'displacement_generator': orm.Dict({'distance': 0.001}),
        'symprec': orm.Float(1e-4),
        'is_symmetry': orm.Bool(False),
        'distinguish_kinds': orm.Bool(True),
        'supercell_matrix': orm.List([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        'primitive_matrix': orm.List([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    }
    generate_preprocess_data(**inputs)

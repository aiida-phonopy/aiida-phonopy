# -*- coding: utf-8 -*-
"""Tests for :mod:`calculations.functions.data_utils`."""
from aiida import orm
import numpy as np
import pytest


@pytest.fixture
def get_structure():
    """Return a Si supercell structure."""

    def _get_structure():
        """Return a Si supercell structure."""
        from ase.build import bulk, make_supercell

        matrix = [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]
        atoms = bulk('Si')
        atoms = make_supercell(atoms, matrix)

        return orm.StructureData(ase=atoms)

    return _get_structure


def test_phonopy_primitive(generate_structure):
    """Test the primitive matrix behaviour of Phonopy."""
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    ase = generate_structure().get_ase()
    ucell = PhonopyAtoms(symbols=ase.get_chemical_symbols(), cell=ase.cell, scaled_positions=ase.get_scaled_positions())
    ph = Phonopy(ucell, primitive_matrix='auto')
    assert ph.primitive_matrix.tolist() == np.eye(3).tolist()


def test_generate_preprocess_data(get_structure):
    """Test the `generate_preprocess_data` calcfunction."""
    from aiida_phonopy.calculations.functions.data_utils import generate_preprocess_data

    inputs = {
        'structure': get_structure(),
        'displacement_generator': orm.Dict({'distance': 0.001}),
        'symprec': orm.Float(1e-4),
        'is_symmetry': orm.Bool(True),
        'distinguish_kinds': orm.Bool(True),
        'supercell_matrix': orm.List([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        'primitive_matrix': orm.Str('auto'),
    }
    preprocess = generate_preprocess_data(**inputs)

    assert preprocess.symprec == 1e-4
    assert preprocess.is_symmetry == True
    assert preprocess.distinguish_kinds == True
    assert preprocess.supercell_matrix == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert preprocess.primitive_matrix != [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


def test_generate_preprocess_data_with_no_sym(get_structure):
    """Test the `generate_preprocess_data` calcfunction `is_symmetry` input behaviour."""
    from aiida_phonopy.calculations.functions.data_utils import generate_preprocess_data

    inputs = {
        'structure': get_structure(),
        'is_symmetry': orm.Bool(False),
    }
    preprocess = generate_preprocess_data(**inputs)

    assert preprocess.symprec == 1e-5
    assert preprocess.is_symmetry == False
    assert preprocess.supercell_matrix == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert preprocess.primitive_matrix == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

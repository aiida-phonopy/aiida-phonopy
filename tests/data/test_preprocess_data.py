# -*- coding: utf-8 -*-
"""Tests for the :mod:`~aiida_phonopy.data` module."""
import pytest


@pytest.mark.usefixtures('aiida_profile')
def test_preprocess_attributes(generate_preprocess_data):
    """Test PreProcessData` attributes."""
    preprocess_data = generate_preprocess_data()
    assert preprocess_data.displacement_dataset == {
        'first_atoms': [{
            'displacement': [0.0, 0.007071067811865475, 0.007071067811865475],
            'number': 0
        }],
        'natom': 2
    }
    assert preprocess_data.displacements == [[0, 0.0, 0.007071067811865475, 0.007071067811865475]]


@pytest.mark.usefixtures('aiida_profile')
def test_preprocess_with_no_symmetries(generate_preprocess_data):
    """Test `PreProcessData` `is_symmetry` behaviour.

    Here we expect that `is_symmetry=False` gives `primitive_matrix = diag(1,1,1)`.
    """
    from numpy import eye

    raw_data = generate_preprocess_data(inputs={'is_symmetry': False})
    assert not raw_data.is_symmetry
    assert raw_data.primitive_matrix.tolist() == eye(3).tolist()


@pytest.mark.usefixtures('aiida_profile')
@pytest.mark.parametrize(
    'pbc', (
        (True, True, True),
        (True, True, False),
        (True, False, True),
        (False, True, True),
        (True, False, False),
        (False, True, False),
        (False, False, True),
    )
)
def test_pbc(generate_preprocess_data, pbc):
    """Test `PreProcessData` using different PBC."""
    preprocess_data = generate_preprocess_data(pbc=pbc)

    assert preprocess_data.pbc == pbc
    assert preprocess_data.get_unitcell().pbc == pbc
    assert preprocess_data.get_primitive_cell().pbc == pbc
    assert preprocess_data.get_supercell().pbc == pbc

    supercells = preprocess_data.get_supercells_with_displacements()

    for structure in supercells.values():
        assert structure.pbc == pbc


@pytest.mark.usefixtures('aiida_profile')
def test_generate_preprocess(generate_structure):
    """Test PreProcessData` methods."""
    from aiida_phonopy.data.preprocess import PreProcessData
    structure = generate_structure()
    PreProcessData.generate_preprocess_data(structure=structure)

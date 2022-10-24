# -*- coding: utf-8 -*-
"""Tests for the :mod:`data` module."""
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
def test_generate_preprocess(generate_structure):
    """Test PreProcessData` methods."""
    from aiida_phonopy.data.preprocess import PreProcessData
    structure = generate_structure()
    PreProcessData.generate_preprocess_data(structure=structure)
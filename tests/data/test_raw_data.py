# -*- coding: utf-8 -*-
"""Tests for the :mod:`data` module."""

import numpy as np
import pytest


@pytest.mark.usefixtures('aiida_profile')
def test_raw_attributes(generate_raw_data):
    """Test `RawData` attributes."""
    raw_data = generate_raw_data()
    param = 5.43
    cell = cell = [[0.0, param / 2.0, param / 2.0], [param / 2.0, 0, param / 2.0], [param / 2.0, param / 2.0, 0.0]]

    assert raw_data.numbers.tolist() == [14, 14]
    assert raw_data.masses.tolist() == [28.0855, 28.0855]
    assert raw_data.positions.tolist() == [[0., 0., 0.], [param / 4.0, param / 4.0, param / 4.0]]
    assert raw_data.cell.tolist() == cell
    assert raw_data.magnetic_moments == None
    assert raw_data.symbols == ['Si', 'Si']
    assert raw_data.names == ['Si', 'Si']
    assert raw_data.supercell_matrix is not None
    assert raw_data.primitive_matrix is not None
    assert raw_data.symprec == 1e-5
    assert raw_data.is_symmetry == True
    assert raw_data.kinds_map == None
    assert raw_data.distinguish_kinds == True
    assert raw_data.dielectric == None
    assert raw_data.born_charges == None
    assert raw_data.has_nac_parameters() == False


@pytest.mark.usefixtures('aiida_profile')
def test_raw_attributes_with_kinds(generate_raw_data):
    """Test `RawData` attributes."""
    raw_data = generate_raw_data(structure_id='silicon-with-names')
    param = 5.43
    cell = cell = [[0.0, param / 2.0, param / 2.0], [param / 2.0, 0, param / 2.0], [param / 2.0, param / 2.0, 0.0]]

    assert raw_data.numbers.tolist() == [1, 2]
    assert raw_data.masses.tolist() == [28.0855, 28.0855]
    assert raw_data.positions.tolist() == [[0., 0., 0.], [param / 4.0, param / 4.0, param / 4.0]]
    assert raw_data.cell.tolist() == cell
    assert raw_data.magnetic_moments == None
    assert raw_data.symbols == ['Si', 'Si']
    assert raw_data.names == ['A', 'B']
    assert raw_data.kinds_map is not None


@pytest.mark.usefixtures('aiida_profile')
def test_raw_methods(generate_raw_data):
    """Test `RawData` methods."""
    from phonopy import Phonopy
    raw_data = generate_raw_data(structure_id='silicon-with-names')
    phonopy_istance = raw_data.get_phonopy_instance()

    assert isinstance(phonopy_istance, Phonopy)


@pytest.mark.usefixtures('aiida_profile')
def test_invalid_nacs(generate_raw_data):
    """Test `RawData` invalid nac attributes."""
    raw_data = generate_raw_data()

    message = 'the array is not of the correct shape'

    with pytest.raises(ValueError) as exception:
        raw_data.set_dielectric([1, 1, 1])

    assert message in str(exception.value)

    with pytest.raises(ValueError) as exception:
        raw_data.set_born_charges(np.eye(3))

    assert message in str(exception.value)


@pytest.mark.usefixtures('aiida_profile')
def test_valid_nacs(generate_raw_data):
    """Test `RawData` valid nac attributes."""
    raw_data = generate_raw_data()
    raw_data.set_dielectric(np.eye(3))
    bcs = np.zeros((2, 3, 3))
    raw_data.set_born_charges(bcs)

    assert raw_data.has_nac_parameters() == True

    raw_data.get_phonopy_instance()
    raw_data.get_phonopy_instance(symmetrize_nac=False)

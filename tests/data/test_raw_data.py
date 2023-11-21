# -*- coding: utf-8 -*-
"""Tests for the :mod:`aiida_phonopy.data` module."""

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
    assert raw_data.magnetic_moments is None
    assert raw_data.symbols == ['Si', 'Si']
    assert raw_data.names == ['Si', 'Si']
    assert raw_data.pbc == (True, True, True)
    assert raw_data.supercell_matrix is not None
    assert raw_data.primitive_matrix is not None
    assert raw_data.symprec == 1e-5
    assert raw_data.is_symmetry
    assert raw_data.kinds_map is None
    assert raw_data.distinguish_kinds
    assert raw_data.dielectric is None
    assert raw_data.born_charges is None
    assert not raw_data.has_nac_parameters()


@pytest.mark.usefixtures('aiida_profile')
def test_raw_with_no_symmetries(generate_raw_data):
    """Test `RawData` `is_symmetry` behaviour.

    Here we expect that `is_symmetry=False` gives `primitive_matrix = diag(1,1,1)`.
    """
    from numpy import eye

    raw_data = generate_raw_data(inputs={'is_symmetry': False})
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
def test_pbc(generate_raw_data, pbc):
    """Test `RawData` with different PBC."""
    raw_data = generate_raw_data(pbc=pbc)

    assert raw_data.pbc == pbc
    assert raw_data.get_unitcell().pbc == pbc
    assert raw_data.get_primitive_cell().pbc == pbc
    assert raw_data.get_supercell().pbc == pbc


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

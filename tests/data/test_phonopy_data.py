# -*- coding: utf-8 -*-
"""Tests for the :mod:`data` module."""

import numpy as np
import pytest


@pytest.mark.usefixtures('aiida_profile')
def test_phonopy_attributes(generate_phonopy_data):
    """Test PhonopyData` attributes."""
    phonopy_data = generate_phonopy_data()
    assert phonopy_data.forces.tolist() == [[[1., 0., 0.], [-1., 0., 0.]]]


@pytest.mark.usefixtures('aiida_profile')
def test_phonopy_invalid_forces(generate_phonopy_data):
    """Test PreProcessData` invalid forces."""

    message = 'the array is not of the correct shape'

    with pytest.raises(ValueError) as exception:
        generate_phonopy_data(forces=[1, 0, 0])

    assert message in str(exception.value)

    with pytest.raises(ValueError) as exception:
        phonopy_data = generate_phonopy_data()
        phonopy_data.set_residual_forces(forces=[1, 0, 0])

    assert message in str(exception.value)


@pytest.mark.usefixtures('aiida_profile')
def test_phonopy_methods(generate_phonopy_data):
    """Test PreProcessData` methods."""
    phonopy_data = generate_phonopy_data()

    phonopy_data.get_phonopy_instance()

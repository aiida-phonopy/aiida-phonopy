# -*- coding: utf-8 -*-
"""Tests for the :mod:`aiida_phonopy.utils.validators` module."""
from aiida.orm import Int
import numpy as np
import pytest

from aiida_phonopy.utils.validators import validate_matrix, validate_nac, validate_positive_integer


@pytest.mark.parametrize(('value', 'message'), (
    ([2, 2, 2], None),
    ([[1, 0, 0], [0, 1, 0], [0, 0, 1]], None),
    (np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]), None),
    ([2, 2], 'need exactly 3 diagonal elements or 3x3 arrays.'),
    ([[1, 0, 0], [0, 0], [0, 0, 1]], 'matrix need to have 3x1 or 3x3 shape.'),
    ([[1, 0, 0], [0, 1, 0], ['a', 0, 1]
      ], "type `<class 'str'>` of a is not an accepted type in matrix; only `int` and `float` are valid."),
))
def test_validate_matrix(value, message):
    """Test the `validate_matrix` function."""
    assert validate_matrix(value, None) == message


def test_validate_positive_integer():
    """Test the `validate_positive_integer` function."""
    assert validate_positive_integer(Int(2), None) is None
    assert validate_positive_integer(Int(-2), None) is not None

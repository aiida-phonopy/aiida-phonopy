# -*- coding: utf-8 -*-
"""Validators for the `aiida-phonopy` plugin."""
from aiida.orm import List

__all__ = ['validate_matrix', 'validate_positive_integer', 'validate_nac']


def validate_matrix(value, _):
    """Validate the `supercell_matrix` and `primitive_matrix` inputs."""
    import numpy as np

    if not isinstance(value, (list, List, np.ndarray)):
        return 'value is not of the right type; only `list`, `aiida.orm.List` and `numpy.ndarray`'

    if isinstance(value, np.ndarray):
        value = value.tolist()

    if not len(value) == 3:
        return 'need exactly 3 diagonal elements or 3x3 arrays.'

    for row in value:
        if isinstance(row, list):
            if not len(row) in [0, 3]:
                return 'matrix need to have 3x1 or 3x3 shape.'
            for element in row:
                if not isinstance(element, (int, float)):
                    return (
                        f'type `{type(element)}` of {element} is not an accepted '
                        'type in matrix; only `int` and `float` are valid.'
                    )


def validate_positive_integer(value, _):
    """Validate that `value` is positive."""
    if not value.value > 0:
        return f'{value} is not positive.'


def validate_nac(value, _):
    """Validate that `value` is positive."""
    try:
        value.get_array('dielectric')
        value.get_array('born_charges')
    except KeyError:
        return 'data does not contain `dieletric` and/or `born_charges` arraynames.'

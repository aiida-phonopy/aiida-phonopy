# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for :class:`aiida_phonopy.calculations.PhonopyCalculation`."""
from aiida import orm
from aiida.common import datastructures
import pytest

from aiida_phonopy.utils.resources import get_default_options


@pytest.fixture
def generate_inputs(
    fixture_code,
    generate_phonopy_data,
):
    """Return a minimal input for the `PhonopyCalculation`."""

    def _generate_inputs(parameters=None, phonopy_options={}):
        """Return a minimal input for the `PhonopyCalculation`."""
        phonopy_data = generate_phonopy_data(**phonopy_options)

        # Generating parameters
        if parameters is None:
            parameters = orm.Dict({})

        # Building input dict to return
        ret_dic = {
            'code': fixture_code('phonopy'),
            'phonopy_data': phonopy_data,
            'parameters': parameters,
            'metadata': {
                'options': get_default_options()
            },
        }

        return ret_dic

    return _generate_inputs


def test_phonopy_default(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test a default `PhonopyCalculation`."""
    entry_point_name = 'phonopy.phonopy'
    inputs = generate_inputs()

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    raw_inputs = ['phonopy.yaml', 'aiida.in']
    retrieve_temporary_list = [
        'phonopy.yaml',
        'force_constants.hdf5',
        'band.hdf5',
        'qpoints.hdf5',
        'mesh.hdf5',
        'irreps.yaml',
        'thermal_properties.yaml',
        'thermal_displacements.yaml',
        'thermal_displacement_matrices.yaml',
        'modulation.yaml',
        'total_dos.dat',
        'projected_dos.dat',
    ]
    retrieve_list = ['aiida.out']

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert len(calc_info.local_copy_list) == 0

    assert len(calc_info.retrieve_list) == len(retrieve_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    assert len(calc_info.retrieve_temporary_list) == len(retrieve_temporary_list)
    assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(raw_inputs)
    assert len(calc_info.codes_info) == 1


def test_phonopy_with_fc(fixture_sandbox, generate_calc_job, generate_inputs, generate_force_constants):
    """Test a `PhonopyCalculation` with force constants."""
    entry_point_name = 'phonopy.phonopy'
    inputs = generate_inputs()
    inputs.pop('phonopy_data')
    inputs['force_constants'] = generate_force_constants()

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    raw_inputs = ['phonopy.yaml', 'aiida.in']
    retrieve_temporary_list = [
        'phonopy.yaml',
        'band.hdf5',
        'qpoints.hdf5',
        'mesh.hdf5',
        'irreps.yaml',
        'thermal_properties.yaml',
        'thermal_displacements.yaml',
        'thermal_displacement_matrices.yaml',
        'modulation.yaml',
        'total_dos.dat',
        'projected_dos.dat',
    ]
    retrieve_list = ['aiida.out']

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert len(calc_info.local_copy_list) == 0

    assert len(calc_info.retrieve_list) == len(retrieve_list)
    assert sorted(calc_info.retrieve_list) == sorted(retrieve_list)

    assert len(calc_info.retrieve_temporary_list) == len(retrieve_temporary_list)
    assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(raw_inputs + ['force_constants.hdf5'])
    assert len(calc_info.codes_info) == 1


def test_phonopy_cmdsline(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test a `PhonopyCalculation` with user-defined cmdline from `parameters`."""
    entry_point_name = 'phonopy.phonopy'
    parameters = orm.Dict({'band': 'auto', 'writedm': False})
    inputs = generate_inputs(parameters=parameters)

    cmd = ['aiida.in']

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    assert len(calc_info.codes_info) == 1
    assert sorted(calc_info.codes_info[0].cmdline_params) == sorted(cmd)

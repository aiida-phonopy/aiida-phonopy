# -*- coding: utf-8 -*-
"""Tests for the `PhonopyPpParser` class."""

from aiida import orm
import pytest

from aiida_phonopy.data import preprocess


@pytest.fixture
def generate_phonopy_calculation_inputs(fixture_code, generate_example_phonopy_data):
    """Return BTO phonopy calculation inputs."""

    def _generate_phonopy_calculation_inputs(parameters=None):
        """Return BTO phonopy calculation inputs."""
        entry_point_name = 'phonopy.phonopy'

        if parameters is None:
            parameters = {}

        inputs = {
            'code': fixture_code(entry_point_name),
            'parameters': orm.Dict(dict=parameters),
            'phonopy_data': generate_example_phonopy_data(),
            'metadata': {
                'options': {
                    'resources': {
                        'num_machines': 1
                    },
                }
            }
        }

        return inputs

    return _generate_phonopy_calculation_inputs


@pytest.fixture
def generate_alas_phonopy_data():
    """Return AlAs phonopy data."""

    def _generate_alas_phonopy_data():
        """Return AlAs phonopy data."""
        import numpy as np

        from aiida_phonopy.data import PhonopyData, PreProcessData

        param = 5.678702301378800
        cell = np.diag((param, param, param)).tolist()

        structure = orm.StructureData(cell=cell)
        structure.append_atom(position=(0., 0., 0.), symbols='Al')
        structure.append_atom(position=(0., 0.5, 0.5), symbols='Al')
        structure.append_atom(position=(0.5, 0., 0.5), symbols='Al')
        structure.append_atom(position=(0.5, 0.5, 0.), symbols='Al')
        structure.append_atom(position=(0.25, 0.25, 0.25), symbols='As')
        structure.append_atom(position=(0.25, 0.75, 0.75), symbols='As')
        structure.append_atom(position=(0.75, 0.75, 0.25), symbols='As')
        structure.append_atom(position=(0.75, 0.25, 0.75), symbols='As')

        preprocess_data = PreProcessData(structure=structure)
        preprocess_data.set_displacements()
        phonopy_data = PhonopyData(preprocess_data=preprocess_data)

        return phonopy_data

    return _generate_alas_phonopy_data


@pytest.fixture
def generate_minimal_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return {
        'parameters': orm.Dict(dict={}),
    }


@pytest.fixture
def generate_band_inputs(generate_alas_phonopy_data):
    """Return only those inputs that the parser will expect to be there."""
    return {
        'phonopy_data': generate_alas_phonopy_data(),
        'parameters': orm.Dict(dict={'BAND': 'auto'}),
    }


@pytest.fixture
def generate_dos_inputs(generate_alas_phonopy_data):
    """Return only those inputs that the parser will expect to be there."""
    return {
        'phonopy_data': generate_alas_phonopy_data(),
        'parameters': orm.Dict(dict={
            'DOS': True,
            'MESH': 50,
            'WRITE_MESH': False
        }),
    }


@pytest.fixture
def generate_pdos_inputs(generate_alas_phonopy_data):
    """Return only those inputs that the parser will expect to be there."""
    return {
        'phonopy_data': generate_alas_phonopy_data(),
        'parameters': orm.Dict(dict={
            'PDOS': 'Al',
            'MESH': 50,
            'WRITE_MESH': False
        }),
    }


def test_phonopy_minimal(generate_calc_job_node, generate_parser, generate_minimal_inputs, tmpdir):
    """Test a `phonopy` calculation where files will be parsed from temporary directory."""
    entry_point_calc_job = 'phonopy.phonopy'
    entry_point_parser = 'phonopy.phonopy'

    attributes = {'retrieve_temporary_list': ['phonopy.yaml']}

    node = generate_calc_job_node(
        entry_point_calc_job,
        test_name='default_minimal',
        inputs=generate_minimal_inputs,
        attributes=attributes,
        retrieve_temporary=(tmpdir, attributes['retrieve_temporary_list']),
    )

    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results


def test_phonopy_bands(generate_calc_job_node, generate_parser, generate_band_inputs, tmpdir):
    """Test a `phonopy` calculation where bands will be parsed from temporary directory."""
    entry_point_calc_job = 'phonopy.phonopy'
    entry_point_parser = 'phonopy.phonopy'

    attributes = {'retrieve_temporary_list': ['phonopy.yaml', 'band.hdf5']}

    node = generate_calc_job_node(
        entry_point_calc_job,
        test_name='default_bands',
        inputs=generate_band_inputs,
        attributes=attributes,
        retrieve_temporary=(tmpdir, attributes['retrieve_temporary_list']),
    )

    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'phonon_bands' in results


def test_phonopy_dos(generate_calc_job_node, generate_parser, generate_dos_inputs, tmpdir):
    """Test a `phonopy` calculation where dos will be parsed from temporary directory."""
    entry_point_calc_job = 'phonopy.phonopy'
    entry_point_parser = 'phonopy.phonopy'

    attributes = {'retrieve_temporary_list': ['phonopy.yaml', 'total_dos.dat']}

    node = generate_calc_job_node(
        entry_point_calc_job,
        test_name='default_dos',
        inputs=generate_dos_inputs,
        attributes=attributes,
        retrieve_temporary=(tmpdir, attributes['retrieve_temporary_list']),
    )

    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'total_phonon_dos' in results


def test_phonopy_pdos(generate_calc_job_node, generate_parser, generate_pdos_inputs, tmpdir):
    """Test a `phonopy` calculation where pdos will be parsed from temporary directory."""
    entry_point_calc_job = 'phonopy.phonopy'
    entry_point_parser = 'phonopy.phonopy'

    attributes = {'retrieve_temporary_list': ['phonopy.yaml', 'projected_dos.dat']}

    node = generate_calc_job_node(
        entry_point_calc_job,
        test_name='default_pdos',
        inputs=generate_pdos_inputs,
        attributes=attributes,
        retrieve_temporary=(tmpdir, attributes['retrieve_temporary_list']),
    )

    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert 'output_parameters' in results
    assert 'projected_phonon_dos' in results


@pytest.mark.parametrize(
    ('parameters', 'test_name', 'output_results', 'temporary_list'),
    (
        ({
            'BAND': 'AUTO'
        }, 'default_outputs_band', ['phonon_bands'], ['band.hdf5']),
        ({
            'TPROP': True,
            'MESH': 50,
            'WRITE_MESH': False,
        }, 'default_outputs_tprop', ['thermal_properties'], ['thermal_properties.yaml']),
        ({
            'TPROP': True,
            'MESH': 50,
            'WRITE_MESH': True,
        }, 'default_outputs_mesh', ['thermal_properties', 'qpoints_mesh'], ['thermal_properties.yaml', 'mesh.hdf5']),
        ({
            'QPOINTS': [1, 1, 1, 2, 2, 2]
        }, 'default_outputs_qpoints', ['qpoints'], ['qpoints.hdf5']),
        ({
            'IRREPS': [0, 0, 0],
        }, 'default_outputs_irreps', ['irreducible_representations'], ['irreps.yaml']),
    ),
)
def test_phonopy_outputs(
    generate_calc_job_node, generate_parser, generate_phonopy_calculation_inputs, parameters, test_name, output_results,
    temporary_list, tmpdir
):
    """Test a `phonopy` calculation with different outputs."""
    entry_point_calc_job = 'phonopy.phonopy'
    entry_point_parser = 'phonopy.phonopy'

    temporary_list += ['phonopy.yaml']
    attributes = {'retrieve_temporary_list': temporary_list}

    node = generate_calc_job_node(
        entry_point_calc_job,
        test_name=test_name,
        inputs=generate_phonopy_calculation_inputs(parameters),
        attributes=attributes,
        retrieve_temporary=(tmpdir, attributes['retrieve_temporary_list']),
    )

    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message

    assert 'output_parameters' in results
    for result in output_results:
        assert result in results


def test_phonopy_files_missing(generate_calc_job_node, generate_parser, generate_band_inputs, tmpdir):
    """Test a `phonopy` calculation where files are missing from temporary directory."""
    entry_point_calc_job = 'phonopy.phonopy'
    entry_point_parser = 'phonopy.phonopy'

    attributes = {'retrieve_temporary_list': ['phonopy.yaml']}

    node = generate_calc_job_node(
        entry_point_calc_job,
        test_name='default_files_missing',
        inputs=generate_band_inputs,
        attributes=attributes,
        retrieve_temporary=(tmpdir, attributes['retrieve_temporary_list']),
    )

    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_FILES_MISSING.status


def test_phonopy_stdout_missing(fixture_localhost, generate_calc_job_node, generate_parser):
    """Test a `phonopy` calculation where the stdout is missing."""
    entry_point_calc_job = 'phonopy.phonopy'
    entry_point_parser = 'phonopy.phonopy'

    node = generate_calc_job_node(
        entry_point_calc_job,
        fixture_localhost,
        'default_stdout_missing',
    )

    parser = generate_parser(entry_point_parser)
    _, calcfunction = parser.parse_from_node(node, store_provenance=False)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_failed, calcfunction.exit_status
    assert calcfunction.exit_status == node.process_class.exit_codes.ERROR_OUTPUT_STDOUT_MISSING.status

# -*- coding: utf-8 -*-
"""Tests for the `PhonopyPpParser` class."""

import pytest

from aiida import orm


@pytest.fixture
def generate_minimal_inputs():
    """Return only those inputs that the parser will expect to be there."""
    return {
        "parameters": orm.Dict(dict={}),
    }


def test_phonopypp_minimal(generate_calc_job_node, generate_parser, generate_minimal_inputs, tmpdir):
    """Test a `phonopy` calculation where files will be parsed from temporary directory."""
    entry_point_calc_job = "phonopy.pp"
    entry_point_parser = "phonopy.pp"

    attributes = {"retrieve_temporary_list": ["phonopy.yaml", "force_constants.hdf5"]}

    node = generate_calc_job_node(
        entry_point_calc_job,
        test_name="default_minimal",
        inputs=generate_minimal_inputs,
        attributes=attributes,
        retrieve_temporary=(tmpdir, attributes["retrieve_temporary_list"]),
    )

    parser = generate_parser(entry_point_parser)
    results, calcfunction = parser.parse_from_node(node, store_provenance=False, retrieved_temporary_folder=tmpdir)

    assert calcfunction.is_finished, calcfunction.exception
    assert calcfunction.is_finished_ok, calcfunction.exit_message
    assert "parameters" in results
    assert "force_constants" in results

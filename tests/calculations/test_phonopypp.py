# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `PhonopyPpCalculation` class."""
import pytest

from aiida import orm
from aiida.common import datastructures

from aiida_phonopy.utils.resources import get_default_options


@pytest.fixture
def generate_inputs(
    fixture_localhost,
    fixture_sandbox,
    fixture_code,
    generate_remote_data,
    generate_structure,
    generate_workchain_force_sets,
):
    """Return a minimal input for the `PhonopyPpCalculation`."""

    def _generate_inputs(parameters=None, remote_data=False):
        import numpy as np
        from aiida_phonopy.calculations.functions.get_force_sets import generate_get_force_sets

        # First we need a pre ran `ForceSetsWorkChain`
        process = generate_workchain_force_sets()
        process.setup()  # we run the phonopy preprocess

        # Generating displacement_dataset
        displacement_dataset = process.outputs["displacement_dataset"]

        # Generating force_sets
        forces_dict = {}
        supercells = process.ctx.supercells
        base_force = np.array([1, 0, 0])

        for key, value in supercells.items():  # key: e.g. "supercell_001", "phonon_supercell_001"
            num = key.split("_")[-1]  # e.g. "001"
            if num == "forces":
                num = 0

            fake_forces = []
            for atom in range(len(value.sites)):
                fake_forces.append(base_force)
                base_force += base_force
            fake_forces = np.array(fake_forces)

            forces = orm.ArrayData()
            forces.set_array("forces", fake_forces)

            forces_dict[f"forces_{num}"] = forces  # fake ArrayData

        get_force_sets = generate_get_force_sets()
        force_sets = get_force_sets(**forces_dict)["force_sets"]

        # Generating parameters
        if parameters is None:
            parameters = orm.Dict(dict={})

        # Generating structure
        structure = generate_structure()

        # Building input dict to return
        ret_dic = {
            "code": fixture_code("phonopy"),
            "structure": structure,
            "supercell_matrix": orm.List(list=[1, 1, 1]),
            "force_sets": force_sets,
            "displacement_dataset": displacement_dataset,
            "parameters": parameters,
            "metadata": {"options": get_default_options()},
        }

        if remote_data:
            ret_dic.update(
                {"parent_folder": generate_remote_data(fixture_localhost, fixture_sandbox.abspath, "phonopy")}
            )

        return ret_dic

    return _generate_inputs


def test_phonopy_default(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test a default `PhonopyPpCalculation`."""
    entry_point_name = "phonopy.pp"
    inputs = generate_inputs()

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    raw_inputs = ["phonopy_cells.yaml", "FORCE_SETS"]
    # retrieve_temporary_list = ["phonopy.yaml", "force_constants.hdf5"]

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert len(calc_info.local_copy_list) == 0
    assert len(calc_info.retrieve_list) == 0
    # assert len(calc_info.retrieve_temporary_list)    == len(retrieve_temporary_list)
    # assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(raw_inputs)
    assert len(calc_info.codes_info) == 1


def test_phonopy_cmdsline(fixture_sandbox, generate_calc_job, generate_inputs):
    """Test a `PhonopyPpCalculation` with user-defined cmdline from `parameters`."""
    entry_point_name = "phonopy.pp"
    parameters = orm.Dict(dict={"band": "auto", "writedm": False})
    inputs = generate_inputs(parameters=parameters)

    cmd_readcell = ["-c", "phonopy_cells.yaml"]
    cmd_writefc = ["--writefc", "--writefc-format=hdf5"]
    cmd_readfc = ["--readfc", "--readfc-format=hdf5"]

    calc_info = generate_calc_job(fixture_sandbox, entry_point_name, inputs)

    assert sorted(calc_info.codes_info[0].cmdline_params) == sorted(cmd_readcell + cmd_writefc)
    assert sorted(calc_info.codes_info[1].cmdline_params) == sorted(["phonopy.conf"] + cmd_readcell + cmd_readfc)
    assert len(calc_info.codes_info) == 2

    raw_inputs = ["phonopy.conf", "phonopy_cells.yaml", "FORCE_SETS"]
    # retrieve_temporary_list = ["phonopy.yaml", "force_constants.hdf5"]

    # Check the attributes of the returned `CalcInfo`
    assert isinstance(calc_info, datastructures.CalcInfo)
    assert len(calc_info.local_copy_list) == 0
    assert len(calc_info.retrieve_list) == 0
    # assert len(calc_info.retrieve_temporary_list)    == len(retrieve_temporary_list)
    # assert sorted(calc_info.retrieve_temporary_list) == sorted(retrieve_temporary_list)

    # Checks on the files written to the sandbox folder as raw input
    assert sorted(fixture_sandbox.get_content_list()) == sorted(raw_inputs)

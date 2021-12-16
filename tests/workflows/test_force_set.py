# -*- coding: utf-8 -*-
# pylint: disable=no-member,redefined-outer-name
"""Tests for the `ForceSetsWorkChain` class."""
import pytest


@pytest.fixture
def generate_displacements():
    """Return a dictionary for `displacements` namespace in ForceSetsWorkChain."""

    def _generate_displacements(type="default"):
        """Return a dictionary for `displacements` namespace in ForceSetsWorkChain."""
        from aiida.orm import Dict

        if type == "default":
            return Dict(dict={"distance": 0.01})
        elif type == "random":
            return Dict(dict={"random_seed": 1, "number_of_snapshots": 5})
        elif type == "raise_error":
            return Dict(dict={"random_seed": 1, "number_of_snapshots": 5, "dataset": {}})
            # TODO add generate_dataset which generates suitable dataset Dict for phonopy

    return _generate_displacements


def test_setup(generate_workchain_force_sets):
    """Test `ForceSetsWorkChain.setup`."""
    process = generate_workchain_force_sets()
    process.setup()

    for key in ["cells", "primitive_matrix", "displacement_dataset"]:
        assert key in process.outputs

    for key in ["primitive", "supercell", "supercell_1"]:
        assert key in process.outputs["cells"]

    assert "primitive" not in process.ctx.supercells
    assert "supercell" not in process.ctx.supercells
    assert "supercell_1" in process.ctx.supercells


def test_setup_with_mapping(generate_workchain_force_sets):
    """Test `ForceSetsWorkChain.setup` when mapping is needed."""
    process = generate_workchain_force_sets(structure_id="silicon-with-names")
    process.setup()

    # assert 'supercells' not in process.ctx.supercells # need to modify
    for key in ["cells", "phonopy_cells", "primitive_matrix", "displacement_dataset"]:
        assert key in process.outputs

    for key in ["primitive", "supercell", "supercell_1", "supercell_2"]:
        assert key in process.outputs["cells"]
        assert key in process.outputs["phonopy_cells"]

    assert "primitive" not in process.ctx.supercells
    assert "supercell" not in process.ctx.supercells
    assert "supercell_1" in process.ctx.supercells
    assert "supercell_2" in process.ctx.supercells


def test_setup_with_subract_forces(generate_workchain_force_sets):
    """Test `ForceSetsWorkChain.setup`."""
    from aiida.orm import Bool

    extra = {"subtract_residual_forces": Bool(True)}
    process = generate_workchain_force_sets(append_inputs=extra)
    process.setup()

    assert "supercell" in process.ctx.supercells

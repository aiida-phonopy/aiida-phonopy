# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name
"""Tests for the `PhonopyPpCalculation` class."""
import numpy as np

from aiida import orm
from aiida_phonopy.calculations.functions.get_force_sets import generate_get_force_sets


def test_get_force_sets(generate_workchain_force_sets):
    """Return only those inputs that the parser will expect to be there."""
    # First we need a pre ran `ForceSetsWorkChain`
    process = generate_workchain_force_sets()
    process.setup()  # we run the phonopy preprocess

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
    output = get_force_sets(**forces_dict)

    assert "force_sets" in output
    assert len(output["force_sets"].get_array("force_sets")) == 1
    assert isinstance(output["force_sets"], orm.ArrayData)

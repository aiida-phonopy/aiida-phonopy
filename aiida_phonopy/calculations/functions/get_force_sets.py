# -*- coding: utf-8 -*-
"""Calcfunctions to gather computed forces into force sets."""

import numpy as np
from aiida.engine import calcfunction
from aiida import orm


def generate_get_force_sets(force_label="forces", force_index=0):
    """Return a `get_force_sets` calcfunction, when the array name for
    the forces is different from the default `forces`.

    :param force_label: label for the forces arrayname
    :param force_index: index for the forces arrayname; e.g. useful when dealing with `TrajectoryData`
    :return: `get_force_sets` calcfunction
    """

    @calcfunction
    def get_force_sets(**forces_dict):
        """Return force sets from supercell forces.

        :param forces_dict: the dictionary must contain keys starting with ``forces`` and ``energy``
        :return: dictionary with force sets and optionally energies, pristine supercell forces and energy
        :raises: KeyError if the input is not in the correct format
        """

        for key in forces_dict.keys():
            if not (key.startswith("energy") or key.startswith("forces")):
                raise KeyError(f"{key} is not an accepted key for the `forces_dict` input.")

        forces_0 = forces_dict.pop("forces_0", None)

        # Setting force sets array
        force_sets = [0 for key in forces_dict if key.startswith("forces")]

        # Setting energies array
        num_energies = len([0 for key in forces_dict if key.startswith("energy")])
        if num_energies > 0:
            energies = np.zeros(num_energies)
        else:
            energies = None

        # Filling arrays
        for key, value in forces_dict.items():
            index = int(key.split("_")[-1])  # e.g. "001" --> 1
            if key.startswith("forces"):
                if force_index == 0:
                    force_sets[index - 1] = value.get_array(force_label)
                else:
                    force_sets[index - 1] = value.get_array(force_label)[force_label]
            elif key.startswith("forces"):
                energies[index - 1] = value.value

        # Finilizing force sets array
        force_sets = np.array(force_sets)
        if forces_0 is not None:
            force_sets = force_sets - forces_0.get_array(force_label)

        # Finalizing the output
        output_dict = {}

        force_sets_data = orm.ArrayData()
        force_sets_data.set_array("force_sets", force_sets)
        output_dict["force_sets"] = force_sets_data

        if energies is not None:
            energies_out = orm.ArrayData()
            energies_out.set_array("energies", energies)
            output_dict["energies"] = energies_out

        return output_dict

    return get_force_sets

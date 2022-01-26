# -*- coding: utf-8 -*-
"""File raw parsers."""

import h5py
import yaml

import numpy as np

from aiida import orm


def parse_yaml(f):
    """Parse a `yaml` file and returns orm.Dict object."""
    yaml_dict = yaml.safe_load(f)

    return orm.Dict(dict=yaml_dict)


def parse_FORCE_CONSTANTS(f):
    """Parse `force_constants.hdf5` output file."""
    with h5py.File(f, "r") as f5:
        force_constants = f5["force_constants"][:]
    fc_array = orm.ArrayData()
    fc_array.set_array("force_constants", force_constants)
    fc_array.label = "force_constants"

    return fc_array


def parse_total_dos(f):
    """Parse `total_dos.dat` output file."""
    data = np.loadtxt(f)
    total_dos = {"frequency_points": data[:, 0], "total_dos": data[:, 1]}
    dos = orm.XyData()
    dos.set_x(total_dos["frequency_points"], "Frequency", "THz")
    dos.set_y(total_dos["total_dos"], "DOS", "1/THz")
    dos.label = "Total DOS"

    return dos


def parse_projected_dos(f):
    """Parse `projected_dos.dat` output file."""
    data = np.loadtxt(f)
    projected_dos = {"frequency_points": data[:, 0], "projected_dos": data[:, 1:].T}
    pdos = orm.XyData()
    pdos_list = [pd for pd in projected_dos["projected_dos"]]
    pdos.set_x(projected_dos["frequency_points"], "Frequency", "THz")
    pdos.set_y(
        pdos_list,
        [
            "Projected DOS",
        ]
        * len(pdos_list),
        [
            "1/THz",
        ]
        * len(pdos_list),
    )
    pdos.label = "Projected DOS"

    return pdos

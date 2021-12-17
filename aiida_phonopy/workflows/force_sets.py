# -*- coding: utf-8 -*-
"""Turn-key solution for automatic frozen phonons calculations."""

from abc import ABCMeta, abstractmethod
from aiida.engine import WorkChain
from aiida import orm
from aiida_phonopy.calculations.functions.preprocess import phonopy_preprocess
from aiida_phonopy.calculations.functions.get_force_sets import generate_get_force_sets


def validate_supercell_matrix(value, _):
    """Validate the `validate_supercell_matrix` input."""
    if not len(value) == 3:
        return "need exactly 3 diagonal elements or 3x3 arrays."

    for row in value:
        if isinstance(row, list):
            if not len(row):
                return "supercell matrix need to have 3x1 or 3x3 shape."
            for element in row:
                if not isinstance(row, int):
                    return (
                        f"type `{type(element)}` of {element} is not an accepted "
                        "type in supercell matrix; only `int` is valid."
                    )
        else:
            if not isinstance(row, int):
                return (
                    f"type `{type(row)}` of {row} is not an accepted type in" "supercell matrix; only `int` is valid."
                )


def validate_positive_integer(value, _):
    """Validate that `value` is positive."""
    if not value.value > 0:
        return f"{value} is not positive."


class ForceSetsWorkChain(WorkChain, metaclass=ABCMeta):
    """
    Workflow to compute automatically the force set of a given structure
    using the frozen phonons approach.

    Phonopy is used to produce structures with displacements,
    while the forces are calculated with a quantum engine of choice.
    """

    _ENABLED_DISPLACEMENTS_FLAGS = {
        "distance": [float],
        "is_plusminus": ["auto", float],
        "is_diagonal": [bool],
        "is_trigonal": [bool],
        "number_of_snapshots": [int, None],
        "random_seed": [int, None],
        "temperature": [float, None],
        "cutoff_frequency": [float, None],
    }

    _DEFAULT_DISPLACEMENTS = {
        "distance": 0.01,
        "is_plusminus": "auto",
        "is_diagonal": True,
        "is_trigonal": False,
        "temperature": None,
        "cutoff_frequency": None,
    }

    _FORCE_LABEL = "forces"
    _FORCE_INDEX = 0

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)
        spec.input("structure", valid_type=orm.StructureData, help="The structure at equilibrium volume.")
        spec.input(
            "supercell_matrix",
            valid_type=orm.List,
            required=True,
            validator=validate_supercell_matrix,
            help=(
                "The matrix used to generate the supercell from the input "
                "structure in the List format. Allowed shapes are 3x1 and 3x3 lists."
            ),
        )
        spec.input(
            "displacements",
            valid_type=orm.Dict,
            default=lambda: orm.Dict(dict=cls._DEFAULT_DISPLACEMENTS),
            help=(
                "Displacements info. The following flags are allowed:\n "
                + "\n ".join(f"{flag_name}" for flag_name in cls._ENABLED_DISPLACEMENTS_FLAGS)
            ),
            validator=cls._validate_displacements,
        )
        spec.input(
            "symmetry_tolerance",
            valid_type=orm.Float,
            default=lambda: orm.Float(1e-5),
            help="Symmetry tolerance for space group analysis on the input structure.",
        )
        spec.input(
            "subtract_residual_forces",
            valid_type=orm.Bool,
            default=lambda: orm.Bool(False),
            help="If `True` it performs a pristine supercell calculation.",
        )  # to explain better?

        spec.outline(
            cls.setup,
            cls.run_forces,
            cls.inspect_forces,
            cls.run_results,
        )

        spec.output_namespace(
            "cells", valid_type=orm.StructureData, dynamic=True, help="The supercells and primitive structures."
        )
        spec.output_namespace(
            "phonopy_cells",
            valid_type=orm.StructureData,
            dynamic=True,
            required=False,
            help=(
                "The structures with the mapped different atoms when kind names are used "
                "for the input structure. Use these ones in the post processing."
            ),
        )
        spec.output("primitive_matrix", valid_type=orm.List, help="The output regarding the primitive matrix.")
        spec.output(
            "displacement_dataset", valid_type=orm.Dict, help="The set of displacements referred to the (super)cell."
        )
        spec.output_namespace(
            "supercells_forces", valid_type=orm.ArrayData, help="The forces acting on the atoms of each (super)cell."
        )
        spec.output(
            "force_sets",
            valid_type=orm.ArrayData,
            required=False,
            help="The sets of forces needed for the phonopy post-process.",
        )
        spec.output("energies", valid_type=orm.ArrayData, required=False, help="The supercells energies.")
        spec.output(
            "cells_mapping",
            valid_type=orm.List,
            required=False,
            help="The mapping between cells for the computation and the `phonopy` cells.",
        )

    @classmethod
    def _validate_displacements(cls, value, _):
        """Validate the ``displacements`` input namespace."""
        if value:
            value_dict = value.get_dict()
            enabled_dict = cls._ENABLED_DISPLACEMENTS_FLAGS
            unknown_flags = set(value_dict.keys()) - set(enabled_dict.keys())
            if unknown_flags:
                return (
                    f"Unknown flags in 'displacements': {unknown_flags}, "
                    f"allowed flags are {cls._ENABLED_DISPLACEMENTS_FLAGS.keys()}."
                )
            invalid_values = [
                value_dict[key]
                for key in value_dict.keys()
                if not (type(value_dict[key]) in enabled_dict[key] or value_dict[key] in enabled_dict[key])
            ]
            if invalid_values:
                return f"Displacement options must be of the correct type; got invalid values {invalid_values}."

    def setup(self):
        """Setup the workflow generating (super)cells and primitive cell."""
        preprocess = phonopy_preprocess(
            self.inputs.structure,
            self.inputs.symmetry_tolerance,
            self.inputs.displacements,
            self.inputs.supercell_matrix,
        )

        for key in preprocess.keys():
            self.out(key, preprocess[key])

        self.ctx.supercells = {}
        cells = preprocess["cells"]
        for item in cells.items():
            self.ctx.supercells.update({item})
        self.ctx.supercells.pop("primitive", None)
        if not self.inputs.subtract_residual_forces.value:
            self.ctx.supercells.pop("supercell")

    @abstractmethod
    def run_forces(self):
        """Run supercells using the selected quantum engine code."""

    @abstractmethod
    def inspect_forces(self):
        """Inspect and collect all runs and return an error message if one or more did not finished correctly."""

    def run_results(self):
        """Collect and expose the outputs."""
        get_force_sets = generate_get_force_sets(force_label=self._FORCE_LABEL, force_index=self._FORCE_INDEX)
        results = get_force_sets(**self.ctx.forces)
        for key, value in results.items():
            self.out(key, value)

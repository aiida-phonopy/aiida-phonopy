# -*- coding: utf-8 -*-
"""CalcJob for phonopy post-process."""

from aiida.engine import CalcJob
from aiida.common import exceptions
from aiida import orm


class BasePhonopyCalculation(CalcJob):
    """Base `CalcJob` implementation for Phonopy and Phono3py post-processing."""

    # it would be nice to group the tags with their relative
    # output file, in order to have a selective parsing
    _AVAILABLE_TAGS = {}  # { 'TAG':[allowed_type], ... }
    _BLOCKED_TAGS = []

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)
        spec.input(
            "parameters",
            valid_type=(orm.Dict, orm.List),
            help=(
                "Phonopy parameters (`setting tags`) for post processing. "
                "The following tags, along their type, are allowed:\n"
                + "\n".join(f"{tag_name}" for tag_name in cls._AVAILABLE_TAGS.keys())
            ),
            validator=cls._validate_parameters,
        )
        spec.input("structure", valid_type=orm.StructureData, required=False, help="Unit cell structure.")
        spec.input(
            "primitive",
            valid_type=orm.StructureData,
            required=False,
            help="Primitive cell structure is necessary only if NAC is applied.",
        )
        spec.input("supercell_matrix", valid_type=orm.List, required=False)
        spec.input("primitive_matrix", valid_type=orm.ArrayData, required=False)
        spec.input("force_sets", valid_type=orm.ArrayData, required=False, help="Sets of forces in supercells.")
        spec.input("symmetry_tolerance", valid_type=orm.Float, required=False, default=lambda: orm.Float(1e-5))
        spec.input(
            "nac_params",
            valid_type=orm.ArrayData,
            required=False,
            help=(
                "Non analytical constants parameters: dielectric constant (`epsilon`) "
                "and effective charges (`born_charges`) tensors."
            ),
        )
        spec.input("displacement_dataset", valid_type=orm.Dict, required=False, help="Displacements dataset.")
        spec.input(
            "dataset", valid_type=(orm.Dict, orm.ArrayData), required=False, help="Displacements dataset with forces."
        )
        spec.input(
            "parent_folder",
            valid_type=(orm.RemoteData, orm.FolderData),
            required=False,
            help="Folder of a phonopy post processing calculations.",
        )
        spec.input("metadata.options.withmpi", valid_type=bool, default=False)

    @classmethod
    def _validate_parameters(cls, value, _):
        """Validate the ``parameters`` input namespace."""

        def __validate_dict(value_dict):
            enabled_dict = cls._AVAILABLE_TAGS
            unknown_tags = set(value_dict.keys()) - set(enabled_dict.keys())
            if unknown_tags:
                return (
                    f"Unknown tags in 'parameters': {unknown_tags}, " f"allowed tags are {cls._AVAILABLE_TAGS.keys()}."
                )
            invalid_values = [
                value_dict[key] for key in value_dict.keys() if not type(value_dict[key]) in enabled_dict[key]
            ]
            if invalid_values:
                return f"Parameters tags must be of the correct type; got invalid values {invalid_values}."

        if value:
            if isinstance(value, orm.Dict):
                __validate_dict(value.get_dict())
            else:
                for val in value.get_list():
                    if isinstance(val, dict):
                        __validate_dict(val)
                    else:
                        return (
                            f"Parameters in List object need to contain dictionarties; got invalid type `{type(val)}`."
                        )


def _case_transform_dict(dictionary, dict_name, func_name, transform):
    from collections import Counter

    if not isinstance(dictionary, dict):
        raise TypeError(f"{func_name} accepts only dictionaries as argument, got {type(dictionary)}")
    new_dict = dict((transform(str(k)), v) for k, v in dictionary.items())
    if len(new_dict) != len(dictionary):
        num_items = Counter(transform(str(k)) for k in dictionary.keys())
        double_keys = ",".join([k for k, v in num_items if v > 1])
        raise exceptions.InputValidationError(
            "Inside the dictionary '{}' there are the following keys that "
            "are repeated more than once when compared case-insensitively: {}."
            "This is not allowed.".format(dict_name, double_keys)
        )
    return new_dict


def _lowercase_dict(dictionary, dict_name):
    return _case_transform_dict(dictionary, dict_name, "_lowercase_dict", str.lower)


def _uppercase_dict(dictionary, dict_name):
    return _case_transform_dict(dictionary, dict_name, "_uppercase_dict", str.upper)

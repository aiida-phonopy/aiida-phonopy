"""Utilities related to process builder or inputs dist."""
import copy
from aiida.engine import calcfunction
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.common import InputValidationError
from aiida.orm import Str, Bool, Code, load_group

KpointsData = DataFactory("array.kpoints")
Dict = DataFactory("dict")
StructureData = DataFactory("structure")
PotcarData = DataFactory("vasp.potcar")


def get_plugin_names(calculator_settings):
    """Return plugin names of calculators."""
    code_strings = []
    if "sequence" in calculator_settings.keys():
        for key in calculator_settings["sequence"]:
            code_strings.append(calculator_settings[key]["code_string"])
    else:
        code_strings.append(calculator_settings["code_string"])

    plugin_names = []
    for code_string in code_strings:
        code = Code.get_from_string(code_string)
        plugin_names.append(code.get_input_plugin_name())

    return plugin_names


def get_calcjob_inputs(
    calculator_settings, structure, calc_type=None, label=None, ctx=None
):
    """Return builder inputs of a calculation."""
    if "sequence" in calculator_settings.keys():
        key = calculator_settings["sequence"][ctx.iteration - 1]
        calculator_inputs = calculator_settings[key]
    else:
        calculator_inputs = calculator_settings

    code = Code.get_from_string(calculator_inputs["code_string"])
    plugin_name = code.get_input_plugin_name()
    if plugin_name == "vasp.vasp":
        if isinstance(calculator_inputs["options"], dict):
            options = Dict(dict=calculator_inputs["options"])
        else:
            options = calculator_inputs["options"]
        if isinstance(calculator_inputs["potential_family"], str):
            potential_family = Str(calculator_inputs["potential_family"])
        else:
            potential_family = calculator_inputs["potential_family"]
        if isinstance(calculator_inputs["potential_mapping"], dict):
            potential_mapping = Dict(dict=calculator_inputs["potential_mapping"])
        else:
            potential_mapping = calculator_inputs["potential_mapping"]
        builder_inputs = {
            "options": options,
            "parameters": _get_parameters_Dict(calculator_inputs),
            "settings": _get_vasp_settings(calculator_inputs),
            "kpoints": _get_kpoints(calculator_inputs, structure),
            "clean_workdir": Bool(False),
            "structure": structure,
            "code": code,
            "potential_family": potential_family,
            "potential_mapping": potential_mapping,
        }
        if label:
            builder_inputs["metadata"] = {"label": label}
    elif plugin_name == "quantumespresso.pw":
        family = load_group(calculator_inputs["pseudo_family_string"])
        pseudos = family.get_pseudos(structure=structure)
        metadata = {"options": calculator_inputs["options"]}
        if label:
            metadata["label"] = label
        pw = {
            "metadata": metadata,
            "parameters": _get_parameters_Dict(calculator_inputs),
            "structure": structure,
            "pseudos": pseudos,
            "code": code,
        }
        builder_inputs = {
            "kpoints": _get_kpoints(calculator_inputs, structure),
            "pw": pw,
        }
    elif plugin_name == "quantumespresso.ph":
        qpoints = KpointsData()
        qpoints.set_kpoints_mesh([1, 1, 1], offset=[0, 0, 0])
        metadata = {"options": calculator_inputs["options"]}
        if label:
            metadata["label"] = label
        ph = {
            "metadata": metadata,
            "qpoints": qpoints,
            "parameters": _get_parameters_Dict(calculator_inputs),
            "parent_folder": ctx.nac_params_calcs[0].outputs.remote_folder,
            "code": code,
        }
        builder_inputs = {"ph": ph}
    else:
        raise RuntimeError("Code could not be found.")

    return builder_inputs


def get_calculator_process(code_string=None, plugin_name=None):
    """Return WorkChain or CalcJob."""
    if plugin_name is None:
        code = Code.get_from_string(code_string)
        _plugin_name = code.get_input_plugin_name()
    else:
        _plugin_name = plugin_name
    if _plugin_name == "vasp.vasp":
        return WorkflowFactory(_plugin_name)
    elif _plugin_name in ("quantumespresso.pw", "quantumespresso.ph"):
        return WorkflowFactory(_plugin_name + ".base")
    else:
        raise RuntimeError("Code could not be found.")


def get_vasp_immigrant_inputs(folder_path, calculator_settings, label=None):
    """Return VASP immigrant inputs.

    folder_path : str
        VASP directory path.
    calculator_settings : dict
        aiida-phonopy calculator settings for forces or nac params.

    """
    code = Code.get_from_string(calculator_settings["code_string"])

    if code.get_input_plugin_name() == "vasp.vasp":
        inputs = {}
        inputs["code"] = code
        inputs["folder_path"] = Str(folder_path)
        if "settings" in calculator_settings:
            settings = copy.deepcopy(calculator_settings["settings"])
        else:
            settings = {}
        if "parser_settings" in calculator_settings:
            if "parser_settings" in settings:
                settings["parser_settings"].update(
                    calculator_settings["parser_settings"]
                )
            else:
                settings["parser_settings"] = calculator_settings["parser_settings"]
        if settings:
            inputs["settings"] = Dict(dict=settings)
        if "options" in calculator_settings:
            inputs["options"] = Dict(dict=calculator_settings["options"])
        if "metadata" in calculator_settings:
            inputs["metadata"] = calculator_settings["metadata"]
            if label:
                inputs["metadata"]["label"] = label
        elif label:
            inputs["metadata"] = {"label": label}
        if "potential_family" in calculator_settings:
            inputs["potential_family"] = Str(calculator_settings["potential_family"])
        if "potential_mapping" in calculator_settings:
            inputs["potential_mapping"] = Dict(
                dict=calculator_settings["potential_mapping"]
            )
    else:
        raise RuntimeError("Code could not be found.")

    return inputs


def _get_parameters_Dict(calculator_inputs):
    """Return parameters for inputs.parameters.

    If calculator_inputs["parameters"] is already a Dict,
    a new Dict will not be made, and just it will be returned.

    """
    if isinstance(calculator_inputs["parameters"], dict):
        return Dict(dict=calculator_inputs["parameters"])
    else:
        return calculator_inputs["parameters"]


def _get_vasp_settings(calculator_inputs):
    """Update VASP settings.

    If no update of settings and calculator_inputs["settings"] is already a Dict,
    a new Dict will not be made, and just it will be returned.

    """
    updated = False
    if "settings" in calculator_inputs.keys():
        settings = calculator_inputs["settings"]
    else:
        settings = {}
    if "parser_settings" in calculator_inputs.keys():
        settings["parser_settings"] = calculator_inputs["parser_settings"]
        updated = True
    if (
        "parser_settings" not in settings
        or "add_forces" not in settings["parser_settings"]
    ):
        settings["parser_settings"].update({"add_forces": True})
        updated = True

    assert settings

    if updated:
        return create_vasp_inputs_settings(Dict(dict=settings))
    else:
        return settings


@calcfunction
def create_vasp_inputs_settings(settings):
    """Store Dict for VaspWorkChain.inputs.settings."""
    return Dict(dict=settings.get_dict())


def _get_kpoints(calculator_inputs, structure):
    """Return KpointsData."""
    if "kpoints" in calculator_inputs.keys():
        assert isinstance(calculator_inputs["kpoints"], KpointsData)
        return calculator_inputs["kpoints"]
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    if "kpoints_density" in calculator_inputs.keys():
        kpoints.set_kpoints_mesh_from_density(calculator_inputs["kpoints_density"])
    elif "kpoints_mesh" in calculator_inputs.keys():
        if "kpoints_offset" in calculator_inputs.keys():
            kpoints_offset = calculator_inputs["kpoints_offset"]
        else:
            kpoints_offset = [0.0, 0.0, 0.0]

        kpoints.set_kpoints_mesh(
            calculator_inputs["kpoints_mesh"], offset=kpoints_offset
        )
    else:
        raise InputValidationError(
            "no kpoint definition in input. "
            "Define either kpoints_density or kpoints_mesh"
        )

    return kpoints

"""Utilities related to process builder or inputs dist."""

from aiida.engine import calcfunction
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.common import InputValidationError
from aiida.orm import Str, Bool, Code, load_group

KpointsData = DataFactory("array.kpoints")
Dict = DataFactory("dict")
StructureData = DataFactory("structure")
PotcarData = DataFactory("vasp.potcar")


def get_calcjob_inputs(
    calculator_settings, structure, calc_type=None, label=None, ctx=None
):
    """Return builder inputs of a calculation."""
    return _get_calcjob_inputs(
        calculator_settings, structure, calc_type=calc_type, label=label, ctx=ctx
    )


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


def _get_calcjob_inputs(
    calculator_settings, structure, calc_type=None, label=None, ctx=None
):
    """Return builder inputs of a calculation."""
    if calc_type is None:
        if "sequence" in calculator_settings.keys():
            key = calculator_settings["sequence"][ctx.iteration - 1]
            settings = calculator_settings[key]
        else:
            settings = calculator_settings
    else:
        settings = Dict(dict=calculator_settings[calc_type])

    code = Code.get_from_string(settings["code_string"])
    plugin_name = code.get_input_plugin_name()
    if plugin_name == "vasp.vasp":
        builder_inputs = {
            "options": _get_vasp_options(settings),
            "parameters": _get_parameters(settings),
            "settings": get_vasp_settings(settings),
            "kpoints": _get_kpoints(settings, structure),
            "clean_workdir": Bool(False),
            "structure": structure,
            "code": code,
        }
        if label:
            builder_inputs.update({"metadata": {"label": label}})
        potential_family = Str(settings["potential_family"])
        potential_mapping = Dict(dict=settings["potential_mapping"])
        builder_inputs.update(
            {
                "potential_family": potential_family,
                "potential_mapping": potential_mapping,
            }
        )
    elif plugin_name == "quantumespresso.pw":
        family = load_group(settings["pseudo_family_string"])
        pseudos = family.get_pseudos(structure=structure)
        pw = {
            "metadata": {"options": _get_options(settings), "label": label},
            "parameters": _get_parameters(settings),
            "structure": structure,
            "pseudos": pseudos,
            "code": code,
        }
        builder_inputs = {"kpoints": _get_kpoints(settings, structure), "pw": pw}
    elif plugin_name == "quantumespresso.ph":
        qpoints = KpointsData()
        qpoints.set_kpoints_mesh([1, 1, 1], offset=[0, 0, 0])
        ph = {
            "metadata": {"options": _get_options(settings), "label": label},
            "qpoints": qpoints,
            "parameters": _get_parameters(settings),
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
        if "parser_settings" in calculator_settings:
            inputs["settings"] = Dict(
                dict={"parser_settings": calculator_settings["parser_settings"]}
            )
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


def _get_options(settings_dict):
    return settings_dict["options"]


def _get_vasp_options(settings):
    return Dict(dict=settings["options"])


def _get_parameters(settings):
    parameters = settings["parameters"]
    return Dict(dict=parameters)


@calcfunction
def get_vasp_settings(settings):
    """Update VASP settings."""
    if "parser_settings" in settings.keys():
        parser_settings_dict = settings["parser_settings"]
    else:
        parser_settings_dict = {}
    if "add_forces" not in parser_settings_dict:
        parser_settings_dict.update({"add_forces": True})
    return Dict(dict={"parser_settings": parser_settings_dict})


def _get_kpoints(settings, structure):
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    if "kpoints_density" in settings.keys():
        kpoints.set_kpoints_mesh_from_density(settings["kpoints_density"])
    elif "kpoints_mesh" in settings.keys():
        if "kpoints_offset" in settings.keys():
            kpoints_offset = settings["kpoints_offset"]
        else:
            kpoints_offset = [0.0, 0.0, 0.0]

        kpoints.set_kpoints_mesh(settings["kpoints_mesh"], offset=kpoints_offset)
    else:
        raise InputValidationError(
            "no kpoint definition in input. "
            "Define either kpoints_density or kpoints_mesh"
        )

    return kpoints

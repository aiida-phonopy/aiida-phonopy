from aiida.plugins import DataFactory, WorkflowFactory, CalculationFactory
from aiida.common import InputValidationError
from aiida.orm import Str, Bool, Code

KpointsData = DataFactory("array.kpoints")
Dict = DataFactory('dict')
StructureData = DataFactory('structure')


def get_calcjob_builder(structure, calculator_settings, calc_type=None,
                        pressure=0.0, label=None):
    if calc_type:
        code = Code.get_from_string(
            calculator_settings[calc_type]['code_string'])
    else:
        code = Code.get_from_string(
            calculator_settings['code_string'])

    if code.attributes['input_plugin'] in ['vasp.vasp']:
        if calc_type is None:
            settings_dict = calculator_settings.get_dict()
        else:
            settings_dict = calculator_settings[calc_type]

        return _get_vasp_builder(structure,
                                 settings_dict,
                                 pressure=pressure,
                                 label=label)
    else:
        raise RuntimeError("Code could not be found.")


def get_immigrant_builder(calculation_folder,
                          calculator_settings,
                          calc_type=None):
    if calc_type:
        code = Code.get_from_string(
            calculator_settings[calc_type]['code_string'])
    else:
        code = Code.get_from_string(
            calculator_settings['code_string'])

    if code.attributes['input_plugin'] in ['vasp.vasp']:
        if calc_type is None:
            settings_dict = calculator_settings.get_dict()
        else:
            settings_dict = calculator_settings[calc_type]

        calc_cls = CalculationFactory('vasp.vasp')
        metadata = {'options': settings_dict['options']}
        process, builder = calc_cls.immigrant(code,
                                              calculation_folder,
                                              metadata=metadata,
                                              settings=settings_dict)
        builder.metadata['options']['parser_name'] = 'vasp.vasp'
    else:
        raise RuntimeError("Code could not be found.")

    return builder


def _get_vasp_builder(structure, settings_dict, pressure=0.0, label=None):
    """
    Generate the input paramemeters needed to run a calculation for VASP

    :param structure:  StructureData object containing the crystal structure
    :param settings:  dict object containing a dictionary with the
        INCAR parameters
    :return: Calculation process object, input dictionary
    """

    code_string = settings_dict['code_string']
    VaspWorkflow = WorkflowFactory('vasp.vasp')
    builder = VaspWorkflow.get_builder()
    if label:
        builder.metadata.label = label
    builder.code = Code.get_from_string(code_string)
    builder.structure = structure
    options = Dict(dict=settings_dict['options'])
    builder.options = options
    builder.clean_workdir = Bool(False)

    if 'parser_settings' in settings_dict:
        parser_settings_dict = settings_dict['parser_settings']
    else:
        parser_settings_dict = {}
    if 'add_forces' not in parser_settings_dict:
        parser_settings_dict.update({'add_forces': True})

    builder.settings = DataFactory('dict')(
        dict={'parser_settings': parser_settings_dict})

    incar = dict(settings_dict['parameters'])
    keys_lower = [key.lower() for key in incar]
    if 'ediff' not in keys_lower:
        incar.update({'EDIFF': 1.0E-8})
    builder.parameters = Dict(dict=incar)
    builder.potential_family = Str(settings_dict['potential_family'])
    builder.potential_mapping = Dict(
        dict=settings_dict['potential_mapping'])

    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    if 'kpoints_density' in settings_dict:
        kpoints.set_kpoints_mesh_from_density(settings_dict['kpoints_density'])
    elif 'kpoints_mesh' in settings_dict:
        if 'kpoints_offset' in settings_dict:
            kpoints_offset = settings_dict['kpoints_offset']
        else:
            kpoints_offset = [0.0, 0.0, 0.0]

        kpoints.set_kpoints_mesh(settings_dict['kpoints_mesh'],
                                 offset=kpoints_offset)
    else:
        raise InputValidationError(
            'no kpoint definition in input. '
            'Define either kpoints_density or kpoints_mesh')

    builder.kpoints = kpoints

    return builder

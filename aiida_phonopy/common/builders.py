from aiida.engine import calcfunction
from aiida.plugins import DataFactory, WorkflowFactory, CalculationFactory
from aiida.common import InputValidationError
from aiida.orm import Str, Bool, Code

KpointsData = DataFactory("array.kpoints")
Dict = DataFactory('dict')
StructureData = DataFactory('structure')
PotcarData = DataFactory('vasp.potcar')

@calcfunction
def get_force_calcjob_inputs(calculator_settings, supercell):
    return _get_calcjob_inputs(calculator_settings, supercell, 'forces')

@calcfunction
def get_phonon_force_calcjob_inputs(calculator_settings, supercell):
    return _get_calcjob_inputs(calculator_settings, supercell, 'phonon_forces')

@calcfunction
def get_nac_calcjob_inputs(calculator_settings, unitcell):
    return _get_calcjob_inputs(calculator_settings, unitcell, 'nac')


def _get_calcjob_inputs(calculator_settings, supercell, calc_type):
    code = Code.get_from_string(calculator_settings[calc_type]['code_string'])
    if code.attributes['input_plugin'] in ['vasp.vasp']:
        settings = calculator_settings[calc_type]
        builder_inputs = {'options': _get_vasp_options(settings),
                          'parameters': _get_vasp_parameters(settings),
                          'settings': _get_vasp_settings(settings),
                          'kpoints': _get_vasp_kpoints(settings, supercell)}
    else:
        raise RuntimeError("Code could not be found.")

    potential_family = Str(settings['potential_family'])
    potential_mapping = Dict(dict=settings['potential_mapping'])
    builder_inputs.update({'potential_family': potential_family,
                           'potential_mapping': potential_mapping})
    return builder_inputs


def get_calcjob_builder(structure, code_string, builder_inputs, label=None):
    code = Code.get_from_string(code_string)
    if code.attributes['input_plugin'] in ['vasp.vasp']:
        VaspWorkflow = WorkflowFactory('vasp.vasp')
        builder = VaspWorkflow.get_builder()
        # VaspCalc = CalculationFactory('vasp.vasp')
        # potential = PotcarData.get_potcars_from_structure(
        #     structure, builder_inputs['potential_family'].value,
        #     mapping=builder_inputs['potential_mapping'].get_dict())
        # builder = VaspCalc.get_builder()
        if label:
            builder.metadata.label = label
        builder.code = Code.get_from_string(code_string)
        builder.structure = structure
        builder.settings = builder_inputs['settings']
        builder.parameters = builder_inputs['parameters']
        builder.kpoints = builder_inputs['kpoints']
        builder.potential_family = builder_inputs['potential_family']
        builder.potential_mapping = builder_inputs['potential_mapping']
        builder.options = builder_inputs['options']
        builder.clean_workdir = Bool(False)
        # builder.potential = potential
        # builder.metadata['options'].update(builder_inputs['options'].get_dict())
    else:
        raise RuntimeError("Code could not be found.")

    return builder


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
        params = {'metadata': {'options': settings_dict['options']},
                  'settings': settings_dict}
        if 'potential_family' in settings_dict:
            params['potential_family'] = settings_dict['potential_family']
        if 'potential_mapping' in settings_dict:
            params['potential_mapping'] = settings_dict['potential_mapping']

        _, builder = calc_cls.immigrant(code, calculation_folder, **params)
        builder.metadata['options']['parser_name'] = 'vasp.vasp'
    else:
        raise RuntimeError("Code could not be found.")

    return builder


def _get_vasp_options(settings_dict):
    return Dict(dict=settings_dict['options'])


def _get_vasp_parameters(settings_dict):
    parameters = settings_dict['parameters']
    incar = parameters['incar']
    keys_lower = [key.lower() for key in incar]
    if 'ediff' not in keys_lower:
        incar.update({'EDIFF': 1.0E-8})
    return Dict(dict=parameters)


def _get_vasp_settings(settings_dict):
    if 'parser_settings' in settings_dict:
        parser_settings_dict = settings_dict['parser_settings']
    else:
        parser_settings_dict = {}
    if 'add_forces' not in parser_settings_dict:
        parser_settings_dict.update({'add_forces': True})
    return Dict(dict={'parser_settings': parser_settings_dict})


def _get_vasp_kpoints(settings_dict, structure):
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

    return kpoints


def _get_vasp_builder(structure, settings_dict, label=None):
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
    builder.options = _get_vasp_options(settings_dict)
    builder.clean_workdir = Bool(False)
    builder.settings = _get_vasp_settings(settings_dict)
    builder.parameters = _get_vasp_parameters(settings_dict)
    builder.potential_family = Str(settings_dict['potential_family'])
    builder.potential_mapping = Dict(
        dict=settings_dict['potential_mapping'])
    builder.kpoints = _get_vasp_kpoints(settings_dict, structure)

    return builder

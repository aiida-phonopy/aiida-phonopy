from aiida.engine import calcfunction
from aiida.plugins import DataFactory, WorkflowFactory, CalculationFactory
from aiida.common import InputValidationError
from aiida.orm import Str, Bool, Code, load_group

KpointsData = DataFactory("array.kpoints")
Dict = DataFactory('dict')
StructureData = DataFactory('structure')
PotcarData = DataFactory('vasp.potcar')


def get_calcjob_inputs(calculator_settings, structure,
                       calc_type=None, label=None, ctx=None):
    """Return builder inputs of a calculation."""
    return _get_calcjob_inputs(calculator_settings, structure,
                               calc_type=calc_type, label=label, ctx=ctx)


@calcfunction
def get_force_calcjob_inputs(calculator_settings, supercell):
    """Return builder inputs of force calculations."""
    return _get_calcjob_inputs(calculator_settings, supercell,
                               calc_type='forces')


@calcfunction
def get_phonon_force_calcjob_inputs(calculator_settings, supercell):
    """Return builder inputs of force calculations for phono3py fc2."""
    return _get_calcjob_inputs(calculator_settings, supercell,
                               calc_type='phonon_forces')


@calcfunction
def get_nac_calcjob_inputs(calculator_settings, unitcell):
    """Return builder inputs of an NAC params calculation."""
    return _get_calcjob_inputs(calculator_settings, unitcell, 'nac')


def _get_calcjob_inputs(calculator_settings, structure, calc_type=None,
                        label=None, ctx=None):
    """Return builder inputs of a calculation."""
    if calc_type is None:
        if 'sequence' in calculator_settings.keys():
            key = calculator_settings['sequence'][ctx.iteration - 1]
            settings = calculator_settings[key]
        else:
            settings = calculator_settings
    else:
        settings = Dict(dict=calculator_settings[calc_type])

    code = Code.get_from_string(settings['code_string'])
    plugin_name = code.get_input_plugin_name()
    if plugin_name == 'vasp.vasp':
        builder_inputs = {'options': _get_vasp_options(settings),
                          'parameters': _get_parameters(settings),
                          'settings': get_vasp_settings(settings),
                          'kpoints': _get_kpoints(settings, structure),
                          'clean_workdir': Bool(False),
                          'structure': structure,
                          'code': code}
        if label:
            builder_inputs.update({'metadata': {'label': label}})
        potential_family = Str(settings['potential_family'])
        potential_mapping = Dict(dict=settings['potential_mapping'])
        builder_inputs.update({'potential_family': potential_family,
                               'potential_mapping': potential_mapping})
    elif plugin_name == 'quantumespresso.pw':
        family = load_group(settings['pseudo_family_string'])
        pseudos = family.get_pseudos(structure=structure)
        pw = {'metadata': {'options': _get_options(settings),
                           'label': label},
              'parameters': _get_parameters(settings),
              'structure': structure,
              'pseudos': pseudos,
              'code': code}
        builder_inputs = {'kpoints': _get_kpoints(settings, structure),
                          'pw': pw}
    elif plugin_name == 'quantumespresso.ph':
        qpoints = KpointsData()
        qpoints.set_kpoints_mesh([1, 1, 1], offset=[0, 0, 0])
        ph = {'metadata': {'options': _get_options(settings),
                           'label': label},
              'qpoints': qpoints,
              'parameters': _get_parameters(settings),
              'parent_folder': ctx['nac_params_1'].outputs.remote_folder,
              'code': code}
        builder_inputs = {'ph': ph}
    else:
        raise RuntimeError("Code could not be found.")

    return builder_inputs


def get_calculator_process(code_string):
    """Return WorkChain or CalcJob."""
    code = Code.get_from_string(code_string)
    plugin_name = code.get_input_plugin_name()
    if plugin_name == 'vasp.vasp':
        return WorkflowFactory(plugin_name)
    elif plugin_name in ('quantumespresso.pw', 'quantumespresso.ph'):
        return WorkflowFactory(plugin_name + ".base")
    else:
        raise RuntimeError("Code could not be found.")


def get_calcjob_builder(structure, code_string, builder_inputs, label=None):
    """Return process builder.

    This method supposes createing a process builder of a force calculator
    (VASP, QE, etc.).

    """
    code = Code.get_from_string(code_string)
    if code.get_input_plugin_name() == 'vasp.vasp':
        VaspWorkflow = WorkflowFactory('vasp.vasp')
        builder = VaspWorkflow.get_builder()
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

    if code.get_input_plugin_name() == 'vasp.vasp':
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


def _get_options(settings_dict):
    return settings_dict['options']


def _get_vasp_options(settings):
    return Dict(dict=settings['options'])


def _get_parameters(settings):
    parameters = settings['parameters']
    return Dict(dict=parameters)


@calcfunction
def get_vasp_settings(settings):
    """Update VASP settings."""
    if 'parser_settings' in settings.keys():
        parser_settings_dict = settings['parser_settings']
    else:
        parser_settings_dict = {}
    if 'add_forces' not in parser_settings_dict:
        parser_settings_dict.update({'add_forces': True})
    return Dict(dict={'parser_settings': parser_settings_dict})


def _get_kpoints(settings, structure):
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    if 'kpoints_density' in settings.keys():
        kpoints.set_kpoints_mesh_from_density(settings['kpoints_density'])
    elif 'kpoints_mesh' in settings.keys():
        if 'kpoints_offset' in settings.keys():
            kpoints_offset = settings['kpoints_offset']
        else:
            kpoints_offset = [0.0, 0.0, 0.0]

        kpoints.set_kpoints_mesh(settings['kpoints_mesh'],
                                 offset=kpoints_offset)
    else:
        raise InputValidationError(
            'no kpoint definition in input. '
            'Define either kpoints_density or kpoints_mesh')

    return kpoints

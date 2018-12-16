import sys
from aiida.orm import Code, CalculationFactory, DataFactory, WorkflowFactory
from aiida.common.exceptions import InputValidationError
from aiida.orm.data.base import Str

KpointsData = DataFactory("array.kpoints")
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')


vasp_configs = {'optimize':
                {'PREC': 'Accurate',
                 'ISTART': 0,
                 'IBRION': 2,
                 'ISIF': 3,
                 'LWAVE': '.FALSE.',
                 'LCHARG': '.FALSE.',
                 'ADDGRID': '.TRUE.',
                 'LREAL': '.FALSE.'},
                'optimize_constant_volume':
                {'PREC': 'Accurate',
                 'ISTART': 0,
                 'IBRION': 2,
                 'ISIF': 4,
                 'LWAVE': '.FALSE.',
                 'LCHARG': '.FALSE.',
                 'ADDGRID': '.TRUE.',
                 'LREAL': '.FALSE.'},
                'forces':
                {'PREC': 'Accurate',
                 'ISYM': 0,
                 'ISTART': 0,
                 'IBRION': -1,
                 'NSW': 0,
                 'LWAVE': '.FALSE.',
                 'LCHARG': '.FALSE.',
                 'ADDGRID': '.TRUE.',
                 'LREAL': '.FALSE.'},
                'born_charges':
                {'PREC': 'Accurate',
                 'LEPSILON': '.TRUE.',
                 'ISTART': 0,
                 'IBRION': -1,
                 'NSW': 0,
                 'LWAVE': '.FALSE.',
                 'LCHARG': '.FALSE.',
                 'ADDGRID': '.TRUE.',
                 'LREAL': '.FALSE.'}}


def generate_vasp_params(structure, settings, calc_type=None, pressure=0.0):
    """
    Generate the input paramemeters needed to run a calculation for VASP

    :param structure:  StructureData object containing the crystal structure
    :param settings:  ParametersData object containing a dictionary with the
        INCAR parameters
    :return: Calculation process object, input dictionary
    """

    if calc_type is None:
        settings_dict = settings.get_dict()
    else:
        settings_dict = settings.get_dict()[calc_type]
    code_name = settings_dict['code']
    code = Code.get_from_string(code_name)
    VaspWorkflow = WorkflowFactory('vasp.vasp')
    builder = VaspWorkflow.get_builder()

    builder.code = code
    builder.structure = structure
    options = ParameterData(dict=settings_dict['options'])
    builder.options = options
    parser_settings_dict = settings_dict['parser_settings']
    if 'parser_settings' in settings_dict:
        parser_settings_dict.update({'add_forces': True})
    else:
        parser_settings_dict = {'add_forces': True}
    builder.settings = DataFactory('parameter')(
        dict={'parser_settings': parser_settings_dict})

    # INCAR (parameters)
    incar = dict(settings_dict['parameters'])
    keys_lower = [key.lower() for key in incar]
    if 'ediff' not in keys_lower:
        incar.update({'EDIFF': 1.0E-8})
    if 'ediffg' not in keys_lower:
        incar.update({'EDIFFG': -1.0E-6})

    builder.parameters = ParameterData(dict=incar)
    builder.potential_family = Str(settings_dict['potential_family'])
    builder.potential_mapping = ParameterData(
        dict=settings_dict['potential_mapping'])

    # KPOINTS
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)

    if 'kpoints_density_{}'.format(calc_type) in settings_dict:
        kpoints.set_kpoints_mesh_from_density(
            settings_dict['kpoints_density_{}'.format(calc_type)])

    elif 'kpoints_density' in settings_dict:
        kpoints.set_kpoints_mesh_from_density(settings_dict['kpoints_density'])

    elif 'kpoints_mesh_{}'.format(calc_type) in settings_dict:
        if 'kpoints_offset' in settings_dict:
            kpoints_offset = settings_dict['kpoints_offset']
        else:
            kpoints_offset = [0.0, 0.0, 0.0]

        kpoints.set_kpoints_mesh(
            settings_dict['kpoints_mesh_{}'.format(calc_type)],
            offset=kpoints_offset)

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


def generate_inputs(structure, es_settings, calc_type=None, pressure=0.0):
    if calc_type:
        plugin = Code.get_from_string(
            es_settings.get_dict()[calc_type]['code']).get_attr('input_plugin')
    else:
        plugin = Code.get_from_string(
            es_settings.get_dict()['code']).get_attr('input_plugin')

    if plugin in ['vasp.vasp']:
        return generate_vasp_params(structure, es_settings,
                                    calc_type=calc_type, pressure=pressure)
    else:
        raise RuntimeError("Code could not be found.")

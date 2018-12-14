import sys
from aiida.orm import Code, CalculationFactory, DataFactory, WorkflowFactory
from aiida.common.exceptions import InputValidationError
from aiida.orm.data.base import Str

KpointsData = DataFactory("array.kpoints")
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')


# Function obtained from aiida's quantumespresso plugin. Copied here
# for convenience
def get_pseudos_qe(structure, family_name):
    """
    Set the pseudo to use for all atomic kinds, picking pseudos from the
    family with name family_name.

    :note: The structure must already be set.

    :param family_name: the name of the group containing the pseudos
    """
    from collections import defaultdict
    from aiida.orm.data.upf import get_pseudos_from_structure

    # A dict {kind_name: pseudo_object}
    kind_pseudo_dict = get_pseudos_from_structure(structure, family_name)

    # We have to group the species by pseudo, I use the pseudo PK
    # pseudo_dict will just map PK->pseudo_object
    pseudo_dict = {}
    # Will contain a list of all species of the pseudo with given PK
    pseudo_species = defaultdict(list)

    for kindname, pseudo in kind_pseudo_dict.iteritems():
        pseudo_dict[pseudo.pk] = pseudo
        pseudo_species[pseudo.pk].append(kindname)

    pseudos = {}
    for pseudo_pk in pseudo_dict:
        pseudo = pseudo_dict[pseudo_pk]
        kinds = pseudo_species[pseudo_pk]
        for kind in kinds:
            pseudos[kind] = pseudo

    return pseudos


def generate_qe_params(structure, settings, pressure=0.0, type=None):
    """Generate the input parameters needed to run a calculation for PW

    :param structure: StructureData object containing the crystal structure
    :param machine: ParametersData object containing a dictionary with the
                    computational resources information
    :param settings: ParametersData object containing a dictionary with the
                     INCAR parameters
    :return: Calculation process object, input dictionary
    """

    try:
        code = settings.dict.code[type]
    except:
        code = settings.dict.code

    plugin = Code.get_from_string(code).get_attr('input_plugin')
    PwCalculation = CalculationFactory(plugin)
    inputs = PwCalculation.process().get_inputs_template()

    # code
    inputs.code = Code.get_from_string(code)

    # structure
    inputs.structure = structure

    # machine
    inputs._options.resources = settings.dict.machine['resources']
    inputs._options.max_wallclock_seconds = settings.dict.machine['max_wallclock_seconds']
    if 'queue_name' in settings.get_dict()['machine']:
        inputs._options.queue_name = settings.dict.machine['queue_name']
    if 'import_sys_environment' in settings.get_dict()['machine']:
        inputs._options.import_sys_environment = settings.dict.machine['import_sys_environment']

    # Parameters
    parameters = dict(settings.dict.parameters)

    parameters['CONTROL'] = {'calculation': 'scf'}

    if type == 'optimize':
        parameters['CONTROL'].update({'calculation': 'vc-relax',
                                      'tstress': True,
                                      'tprnfor': True,
                                      'etot_conv_thr': 1.e-8,
                                      'forc_conv_thr': 1.e-8})
        parameters['CELL'] = {'press': pressure,
                              'press_conv_thr': 1.e-3,
                              'cell_dynamics': 'bfgs',  # Quasi-Newton algorithm
                              #   'cell_dofree': 'all'
                              }  # Degrees of movement
        parameters['IONS'] = {'ion_dynamics': 'bfgs',
                              'ion_nstepe': 10}

    if type == 'forces':
        parameters['CONTROL'].update({'tstress': True,
                                      'tprnfor': True,
                                      'etot_conv_thr': 1.e-8,
                                      'forc_conv_thr': 1.e-8
                                      })

    if type == 'born_charges':  # in development (not really usable)
        parameters['CONTROL'].update({'tstress': True,
                                      'tprnfor': True,
                                      'etot_conv_thr': 1.e-8,
                                      'forc_conv_thr': 1.e-8
                                      })

        #parameters['INPUTPH'] = {'epsil': True,
        #                         'zeu': True}  # Degrees of movement

    inputs.parameters = ParameterData(dict=parameters)

    # Kpoints
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints_mesh_from_density(settings.dict.kpoints_density)

    inputs.kpoints = kpoints

    inputs.pseudo = get_pseudos_qe(structure, settings.dict.pseudos_family)

    return PwCalculation.process(), inputs


def generate_lammps_params(structure, settings, type=None, pressure=0.0):
    """
    Generate the input paramemeters needed to run a calculation for LAMMPS

    :param structure: StructureData object
    :param settings: ParametersData object containing a dictionary with the
                     LAMMPS parameters
    :return: Calculation process object, input dictionary
    """

    try:
        code = settings.dict.code[type]
    except:
        code = settings.dict.code

    plugin = Code.get_from_string(code).get_attr('input_plugin')
    LammpsCalculation = CalculationFactory(plugin)
    inputs = LammpsCalculation.process().get_inputs_template()
    inputs.code = Code.get_from_string(code)

    # machine
    inputs._options.resources = settings.dict.machine['resources']
    inputs._options.max_wallclock_seconds = settings.dict.machine[
        'max_wallclock_seconds']

    if 'queue_name' in settings.get_dict()['machine']:
        inputs._options.queue_name = settings.dict.machine['queue_name']
    if 'import_sys_environment' in settings.get_dict()['machine']:
        inputs._options.import_sys_environment = settings.dict.machine[
            'import_sys_environment']

    inputs.structure = structure
    inputs.potential = ParameterData(dict=settings.dict.potential)

    if type == 'forces':
        if 'parameters' in settings.get_dict():
            lammps_parameters = dict(settings.dict.parameters)
            inputs.parameters = ParameterData(dict=lammps_parameters)

    # if code.get_input_plugin_name() == 'lammps.optimize':
    if type == 'optimize':
        print('optimize inside')

        lammps_parameters = dict(settings.dict.parameters)
        lammps_parameters.update({'pressure': pressure})  # pressure kb
        inputs.parameters = ParameterData(dict=lammps_parameters)

    return LammpsCalculation.process(), inputs


def generate_vasp_params(structure, settings, calc_type=None, pressure=0.0):
    """
    Generate the input paramemeters needed to run a calculation for VASP

    :param structure:  StructureData object containing the crystal structure
    :param settings:  ParametersData object containing a dictionary with the
        INCAR parameters
    :return: Calculation process object, input dictionary
    """

    code_name = settings.get_dict()['code']
    if calc_type in code_name:
        code_name = settings.get_dict()['code'][calc_type]
    code = Code.get_from_string(code_name)
    VaspWorkflow = WorkflowFactory('vasp.vasp')
    builder = VaspWorkflow.get_builder()
    builder.code = code
    builder.structure = structure
    options = ParameterData(dict=settings.get_dict()['options'])
    builder.options = options
    parser_settings_dict = settings.get_dict()['parser_settings']
    if 'parser_settings' not in settings.get_dict():
        parser_settings_dict = {'add_forces': True}
    else:
        parser_settings_dict.update({'add_forces': True})
    settings_dict = {'parser_settings': parser_settings_dict}
    builder.settings = DataFactory('parameter')(dict=settings_dict)

    # INCAR (parameters)
    incar = dict(settings.get_dict()['parameters'])

    if calc_type == 'optimize':
        incar.update({
            'PREC': 'Accurate',
            'ISTART': 0,
            'IBRION': 2,
            'ISIF': 3,
            'LWAVE': '.FALSE.',
            'LCHARG': '.FALSE.',
            'ADDGRID': '.TRUE.',
            'LREAL': '.FALSE.',
            'PSTRESS': pressure})  # unit: kb -> kB

        if 'NSW' not in incar:
            incar.update({'NSW': 300})

    elif calc_type == 'optimize_constant_volume':
        incar.update({
            'PREC': 'Accurate',
            'ISTART': 0,
            'IBRION': 2,
            'ISIF': 4,
            'LWAVE': '.FALSE.',
            'LCHARG': '.FALSE.',
            'ADDGRID': '.TRUE.',
            'LREAL': '.FALSE.'})

        if 'NSW' not in incar:
            incar.update({'NSW': 300})

    elif calc_type == 'forces':
        incar.update({
            'PREC': 'Accurate',
            'ISYM': 0,
            'ISTART': 0,
            'IBRION': -1,
            'NSW': 0,
            'LWAVE': '.FALSE.',
            'LCHARG': '.FALSE.',
            'ADDGRID': '.TRUE.',
            'LREAL': '.FALSE.'})

    elif calc_type == 'born_charges':
        incar.update({
            'PREC': 'Accurate',
            'LEPSILON': '.TRUE.',
            'ISTART': 0,
            'IBRION': -1,
            'NSW': 0,
            'LWAVE': '.FALSE.',
            'LCHARG': '.FALSE.',
            'ADDGRID': '.TRUE.',
            'LREAL': '.FALSE.'})

    if 'EDIFF' not in incar:
        incar.update({'EDIFF': 1.0E-8})
    if 'EDIFFG' not in incar:
        incar.update({'EDIFFG': -1.0E-6})

    builder.parameters = ParameterData(dict=incar)
    builder.potential_family = Str(settings.get_dict()['potential_family'])
    builder.potential_mapping = ParameterData(
        dict=settings.get_dict()['potential_mapping'])

    # KPOINTS
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)

    if 'kpoints_density_{}'.format(calc_type) in settings.get_dict():
        kpoints.set_kpoints_mesh_from_density(
            settings.get_dict()['kpoints_density_{}'.format(calc_type)])

    elif 'kpoints_density' in settings.get_dict():
        kpoints.set_kpoints_mesh_from_density(settings.dict.kpoints_density)

    elif 'kpoints_mesh_{}'.format(calc_type) in settings.get_dict():
        if 'kpoints_offset' in settings.get_dict():
            kpoints_offset = settings.dict.kpoints_offset
        else:
            kpoints_offset = [0.0, 0.0, 0.0]

        kpoints.set_kpoints_mesh(
            settings.get_dict()['kpoints_mesh_{}'.format(calc_type)],
            offset=kpoints_offset)

    elif 'kpoints_mesh' in settings.get_dict():
        if 'kpoints_offset' in settings.get_dict():
            kpoints_offset = settings.dict.kpoints_offset
        else:
            kpoints_offset = [0.0, 0.0, 0.0]

        kpoints.set_kpoints_mesh(settings.dict.kpoints_mesh,
                                 offset=kpoints_offset)
    else:
        raise InputValidationError(
            'no kpoint definition in input. '
            'Define either kpoints_density or kpoints_mesh')

    builder.kpoints = kpoints

    return builder, None


def generate_inputs(structure, es_settings, calc_type=None, pressure=0.0):

    try:
        plugin = Code.get_from_string(
            es_settings.get_dict()['code'][calc_type]).get_attr('input_plugin')
    except:
        try:
            plugin = Code.get_from_string(
                es_settings.get_dict()['code']).get_attr('input_plugin')
        except InputValidationError:
            raise InputValidationError(
                'No code provided for {} calculation type'.format(calc_type))

    if plugin in ['vasp.vasp']:
        return generate_vasp_params(structure, es_settings,
                                    calc_type=calc_type, pressure=pressure)

    elif plugin in ['quantumespresso.pw']:
        return generate_qe_params(structure, es_settings,
                                  calc_type=calc_type, pressure=pressure)

    elif plugin in ['lammps.force', 'lammps.optimize', 'lammps.md']:
        return generate_lammps_params(structure, es_settings,
                                      calc_type=calc_type, pressure=pressure)
    else:
        sys.exit(0)

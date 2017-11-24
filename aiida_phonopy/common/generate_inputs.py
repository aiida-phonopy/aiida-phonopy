# This file is used to generate inputs for the 3 different supported codes: VASP, QE and LAMMPS. This is implemented
# in 3 functions: generate_vasp_params(), generate_qe_params() and generate_lammps_params().
# generate_inputs() function at the end of the file decides which function to use according to the plugin

from aiida.work.run import run, submit, async
from aiida.orm import Code, CalculationFactory, DataFactory
# from aiida.orm.data.upf import UpfData


KpointsData = DataFactory("array.kpoints")
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')

# Function obtained from aiida's quantumespresso plugin. Copied here for convinence
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

    """
    Generate the input paramemeters needed to run a calculation for PW (Quantum Espresso)

    :param structure:  StructureData object containing the crystal structure
    :param machine:  ParametersData object containing a dictionary with the computational resources information
    :param settings:  ParametersData object containing a dictionary with the INCAR parameters
    :return: Calculation process object, input dictionary
    """

    if type is None:
        code = settings.dict.code
    else:
        code = settings.dict.code[type]

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

        parameters['CONTROL'].update({'tstress': True,
                                      'tprnfor': True})

    if type == 'born_charges':
        parameters['CONTROL'].update({'tstress': True,
                                      'tprnfor': True
                                      })
        parameters['INPUTPH'] = {'epsil': True,
                                 'zeu': True}  # Degrees of movement

    inputs.parameters = ParameterData(dict=parameters)

    # Kpoints
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints_mesh_from_density(settings.dict.kpoints_density)

    inputs.kpoints = kpoints

    inputs.pseudo = get_pseudos_qe(structure, settings.dict.pseudos_family)

    return PwCalculation.process(), inputs


def generate_lammps_params(structure, settings, pressure=0.0, type=None):
    """
    Generate the input paramemeters needed to run a calculation for LAMMPS

    :param structure: StructureData object
    :param settings: ParametersData object containing a dictionary with the LAMMPS parameters
    :return: Calculation process object, input dictionary
    """

    if type is None:
        code = settings.dict.code
    else:
        code = settings.dict.code[type]

    plugin = Code.get_from_string(code).get_attr('input_plugin')
    LammpsCalculation = CalculationFactory(plugin)
    inputs = LammpsCalculation.process().get_inputs_template()
    inputs.code = Code.get_from_string(code)

    inputs._options.resources = settings.dict.machine['resources']
    inputs._options.max_wallclock_seconds = settings.dict.machine['max_wallclock_seconds']

    inputs.structure = structure
    inputs.potential = ParameterData(dict=settings.dict.potential)

    # if code.get_input_plugin_name() == 'lammps.optimize':
    if type == 'optimize':
        lammps_parameters = dict(settings.dict.parameters)
        lammps_parameters.update({'pressure': pressure})  # pressure kb
        inputs.parameters = ParameterData(dict=lammps_parameters)

    return LammpsCalculation.process(), inputs


def get_pseudos_vasp(structure, family_name, folder_path=None):
    """
    Set the pseudo to use for all atomic kinds, picking pseudos from the
    family with name family_name.

    :note: The structure must already be set.

    :param family_name: the name of the group containing the pseudos
    """
    import numpy as np

    PawData = DataFactory('vasp.paw')

    unique_symbols = np.unique([site.kind_name for site in structure.sites]).tolist()

    paw_cls = PawData()
    if folder_path is not None:
        paw_cls.import_family(folder_path,
                              familyname=family_name,
                              family_desc='This is a test family',
                              # store=True,
                              stop_if_existing=False
                              )

    pseudos = {}
    for symbol in unique_symbols:
        pseudos[symbol] = paw_cls.load_paw(family=family_name,
                                           symbol=symbol)[0]

    return pseudos


def generate_vasp_params(structure, settings, type=None, pressure=0.0):
    """
    Generate the input paramemeters needed to run a calculation for VASP

    :param structure:  StructureData object containing the crystal structure
    :param settings:  ParametersData object containing a dictionary with the INCAR parameters
    :return: Calculation process object, input dictionary
    """


    if type is None:
        code = settings.dict.code
    else:
        code = settings.dict.code[type]

    plugin = Code.get_from_string(code).get_attr('input_plugin')

    VaspCalculation = CalculationFactory(plugin)

    class VaspCalculation2(VaspCalculation):
            _default_parser = 'vasp.pymatgen'

    inputs = VaspCalculation2.process().get_inputs_template()

    # code
    inputs.code = Code.get_from_string(code)

    # structure
    inputs.structure = structure

    # machine
    inputs._options.resources = settings.dict.machine['resources']
    inputs._options.max_wallclock_seconds = settings.dict.machine['max_wallclock_seconds']

    # INCAR (parameters)
    incar = dict(settings.dict.parameters)

    if type == 'optimize':
        incar.update({
            'PREC': 'Accurate',
            'ISTART': 0,
            'IBRION': 2,
            'ISIF': 3,
            'NSW': 100,
            'LWAVE': '.FALSE.',
            'LCHARG': '.FALSE.',
            'EDIFF': 0,
            'EDIFFG': -1e-08,
            'ADDGRID': '.TRUE.',
            'LREAL': '.FALSE.',
            'PSTRESS': pressure})  # unit: kb -> kB

    elif type == 'optimize_constant_volume':
        incar.update({
            'PREC': 'Accurate',
            'ISTART': 0,
            'IBRION': 2,
            'ISIF': 4,
            'NSW': 100,
            'LWAVE': '.FALSE.',
            'LCHARG': '.FALSE.',
            'EDIFF': 1e-08,
            'EDIFFG': -1e-08,
            'ADDGRID': '.TRUE.',
            'LREAL': '.FALSE.'})

    elif type == 'forces':
        incar.update({
            'PREC': 'Accurate',
            'ISYM': 0,
            'ISTART': 0,
            'IBRION': -1,
            'NSW': 1,
            'LWAVE': '.FALSE.',
            'LCHARG': '.FALSE.',
            'EDIFF': 1e-08,
            'ADDGRID': '.TRUE.',
            'LREAL': '.FALSE.'})

    elif type == 'born_charges':
        incar.update({
            'PREC': 'Accurate',
            'LEPSILON': '.TRUE.',
            'ISTART': 0,
            'IBRION': 1,
            'NSW': 0,
            'LWAVE': '.FALSE.',
            'LCHARG': '.FALSE.',
            'EDIFF': 1e-08,
            'ADDGRID': '.TRUE.',
            'LREAL': '.FALSE.'})

    inputs.parameters = ParameterData(dict=incar)

    # POTCAR (pseudo potentials)
    inputs.paw = get_pseudos_vasp(structure, settings.dict.pseudos_family,
                                  folder_path=settings.dict.family_folder)


    # Kpoints
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    kpoints.set_kpoints_mesh_from_density(settings.dict.kpoints_density)
    inputs.kpoints = kpoints

    return VaspCalculation.process(), inputs

def generate_inputs(structure,  es_settings, type=None, pressure=0.0, machine=None):

    try:
        plugin = Code.get_from_string(es_settings.dict.code[type]).get_attr('input_plugin')
    except:
        plugin = Code.get_from_string(es_settings.dict.code).get_attr('input_plugin')

    if plugin in ['vasp.vasp']:
        return generate_vasp_params(structure, es_settings, type=type, pressure=pressure)

    elif plugin in ['quantumespresso.pw']:
        return generate_qe_params(structure, es_settings, type=type, pressure=pressure)

    elif plugin in ['lammps.force', 'lammps.optimize', 'lammps.md']:
        return generate_lammps_params(structure, es_settings, type=type, pressure=pressure)
    else:
        print ('No supported plugin')
        exit()

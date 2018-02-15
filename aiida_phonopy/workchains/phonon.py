# Works run by the daemon (using submit)

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction

from aiida.orm import Code, CalculationFactory, load_node, DataFactory, WorkflowFactory
from aiida.work.run import run, submit

from aiida.orm.data.base import Str, Float, Bool
from aiida.work.workchain import _If, _While

import numpy as np
from aiida_phonopy.common.generate_inputs import generate_inputs

# Should be improved by some kind of WorkChainFactory
# For now all workchains should be copied to aiida/workflows
# from aiida.workflows.wc_optimize import OptimizeStructure

ForceConstantsData = DataFactory('phonopy.force_constants')
ForceSetsData = DataFactory('phonopy.force_sets')
BandStructureData = DataFactory('phonopy.band_structure')
PhononDosData = DataFactory('phonopy.phonon_dos')
NacData = DataFactory('phonopy.nac')

ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

OptimizeStructure = WorkflowFactory('phonopy.optimize')

__testing__ = False

def generate_phonopy_params(structure, ph_settings, force_sets=None, force_constants=None, nac_data=None, bands=None):
    """
    Generate inputs parameters needed to do a remote phonopy calculation

    :param structure: StructureData Object that constains the crystal structure unit cell
    :param ph_settings: ParametersData object containing a dictionary with the phonopy input data
    :param force_sets: ForceSetssData object containing the atomic forces and displacement information
    :return: Calculation process object, input dictionary
    """

    try:
        code = Code.get_from_string(ph_settings.dict.code['fc2'])
    except :
        code = Code.get_from_string(ph_settings.dict.code)

    plugin = code.get_attr('input_plugin')
    PhonopyCalculation = CalculationFactory(plugin)

    # The inputs
    inputs = PhonopyCalculation.process().get_inputs_template()

    # code
    inputs.code = code

    # structure
    inputs.structure = structure

    # parameters
    inputs.parameters = ph_settings

    # resources
    inputs._options.resources = ph_settings.dict.machine['resources']
    inputs._options.max_wallclock_seconds = ph_settings.dict.machine['max_wallclock_seconds']

    # data_sets
    inputs.data_sets = force_sets

    # data_sets
    if force_sets is not None:
        inputs.data_sets = force_sets

    # force constants
    if force_constants is not None:
        inputs.force_constants = force_constants

    # non-analytical corrections
    if nac_data is not None:
        inputs.nac_data = nac_data

    # bands
    if bands is not None:
        inputs.bands = bands

    if force_constants is None and force_sets is None:
        Exception('Either force sets or force constants must be set!')

    return PhonopyCalculation.process(), inputs


def wf_like_calculation(work_function):
    """
    This function defines decorator to emulate the output stored in self.ctx of remote calculations
    in local workfunctions
    :param work_function: @workfunction decorated function
    :return: corresponding workcalculation function
    """
    def work_calculation(*args, **kwargs):
        return work_function(*args, **kwargs).values()[0].get_inputs_dict().values()[0]
    return work_calculation


def phonopy_bulk_from_structure(structure):
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
    bulk = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)
    return bulk


@workfunction
def create_supercells_with_displacements_using_phonopy(structure, ph_settings):
    """
    Use phonopy to create the supercells with displacements to calculate the force constants by using
    finite displacements methodology

    :param structure: StructureData object
    :param phonopy_input: ParametersData object containing a dictionary with the data needed for phonopy
    :return: A set of StructureData Objects containing the supercells with displacements
    """
    from phonopy import Phonopy

    # Generate phonopy phonon object
    phonon = Phonopy(phonopy_bulk_from_structure(structure),
                     supercell_matrix=ph_settings.dict.supercell,
                     primitive_matrix=ph_settings.dict.primitive,
                     symprec=ph_settings.dict.symmetry_precision)

    phonon.generate_displacements(distance=ph_settings.dict.distance)

    cells_with_disp = phonon.get_supercells_with_displacements()

    # Transform cells to StructureData and set them ready to return
    data_sets = phonon.get_displacement_dataset()
    data_sets_object = ForceSetsData(data_sets=data_sets)

    disp_cells = {'data_sets': data_sets_object}
    for i, phonopy_supercell in enumerate(cells_with_disp):
        supercell = StructureData(cell=phonopy_supercell.get_cell())
        for symbol, position in zip(phonopy_supercell.get_chemical_symbols(),
                                    phonopy_supercell.get_positions()):
            supercell.append_atom(position=position, symbols=symbol)
        disp_cells["structure_{}".format(i)] = supercell

    return disp_cells


@workfunction
def create_forces_set(**kwargs):
    """
    Build data_sets from forces of supercells with displacments

    :param forces_X: ArrayData objects that contain the atomic forces for each supercell with displacement, respectively (X is integer)
    :param data_sets: ForceSetsData object that contains the displacements info (This info should match with forces_X)
    :return: ForceSetsData object that contains the atomic forces and displacements info (datasets dict in phonopy)

    """
    data_sets = kwargs.pop('data_sets')
    force_sets = ForceSetsData(data_sets=data_sets.get_data_sets())

    forces = []
    for i in range(data_sets.get_number_of_displacements()):
        forces.append(kwargs.pop('forces_{}'.format(i)).get_array('forces')[-1])

    force_sets.set_forces(forces)

    return {'force_sets': force_sets}

@workfunction
def get_force_constants_from_phonopy(structure, ph_settings, force_sets):
    """
    Calculate the force constants locally using phonopy

    :param structure:
    :param phonopy_input: ParameterData object that contains phonopy settings
    :param force_sets: ForceSetsData object that contains the atomic forces and displacements info (datasets dict in phonopy)
    :return: ForceConstantsData object containing the 2nd order force constants calculated with phonopy
    """
    from phonopy import Phonopy

    # Generate phonopy phonon object

    phonon = Phonopy(phonopy_bulk_from_structure(structure),
                     ph_settings.dict.supercell,
                     primitive_matrix=ph_settings.dict.primitive,
                     symprec=ph_settings.dict.symmetry_precision)

    phonon.generate_displacements(distance=ph_settings.dict.distance)

    # Build data_sets from forces of supercells with displacments
    phonon.set_displacement_dataset(force_sets.get_force_sets())
    phonon.produce_force_constants()

    force_constants = ForceConstantsData(data=phonon.get_force_constants())

    return {'force_constants': force_constants}


@workfunction
def get_nac_from_data(**kwargs):
    """
    Worfunction to extract nac ArrayData object from calc
    """

    nac_data = NacData(structure=kwargs.pop('structure'),
                       born_charges=kwargs.pop('born_charges').get_array('born_charges'),
                       epsilon=kwargs.pop('epsilon').get_array('epsilon')[-1])

    return {'nac_data': nac_data}


def get_path_using_seekpath(structure, band_resolution=30):
    import seekpath

    phonopy_structure = phonopy_bulk_from_structure(structure)
    cell = phonopy_structure.get_cell()
    scaled_positions = phonopy_structure.get_scaled_positions()
    numbers = phonopy_structure.get_atomic_numbers()

    structure = (cell, scaled_positions, numbers)
    path_data = seekpath.get_path(structure)

    labels = path_data['point_coords']

    band_ranges = []
    for set in path_data['path']:
        band_ranges.append([labels[set[0]], labels[set[1]]])

    bands =[]
    for q_start, q_end in band_ranges:
        band = []
        for i in range(band_resolution+1):
            band.append(np.array(q_start) + (np.array(q_end) - np.array(q_start)) / band_resolution * i)
        bands.append(band)

    band_structure = BandStructureData(bands=bands,
                                       labels=path_data['path'],
                                       unitcell=phonopy_structure.get_cell())
    return band_structure


@workfunction
def get_primitive(structure, ph_settings):

    from phonopy import Phonopy

    phonon = Phonopy(phonopy_bulk_from_structure(structure),
                     supercell_matrix=ph_settings.dict.supercell,
                     primitive_matrix=ph_settings.dict.primitive,
                     symprec=ph_settings.dict.symmetry_precision)

    primitive_phonopy = phonon.get_primitive()

    primitive_cell = primitive_phonopy.get_cell()
    symbols = primitive_phonopy.get_chemical_symbols()
    positions = primitive_phonopy.get_positions()

    primitive_structure = StructureData(cell=primitive_cell)
    for symbol, position in zip(symbols, positions):
        primitive_structure.append_atom(position=position, symbols=symbol)

    return {'primitive_structure': primitive_structure}


@wf_like_calculation
@workfunction
def get_properties_from_phonopy(**kwargs):
    """
    Calculate DOS and thermal properties using phonopy (locally)
    :param structure: StructureData Object
    :param ph_settings: Parametersdata object containing a dictionary with the data needed to run phonopy:
            supercells matrix, primitive matrix and q-points mesh.
    :param force_constants: (optional)ForceConstantsData object containing the 2nd order force constants
    :param force_sets: (optional) ForceSetsData object containing the phonopy force sets
    :param nac: (optional) ArrayData object from a single point calculation data containing dielectric tensor and Born charges
    :return: phonon band structure, force constants, thermal properties and DOS
    """

    structure = kwargs.pop('structure')
    ph_settings = kwargs.pop('ph_settings')
    bands = kwargs.pop('bands')

    from phonopy import Phonopy

    phonon = Phonopy(phonopy_bulk_from_structure(structure),
                     supercell_matrix=ph_settings.dict.supercell,
                     primitive_matrix=ph_settings.dict.primitive,
                     symprec=ph_settings.dict.symmetry_precision)

    if 'force_constants' in kwargs:
        force_constants = kwargs.pop('force_constants')
        phonon.set_force_constants(force_constants.get_data())

    else:
        force_sets = kwargs.pop('force_sets')
        phonon.set_displacement_dataset(force_sets.get_force_sets())
        phonon.produce_force_constants()
        force_constants = ForceConstantsData(data=phonon.get_force_constants())

    if 'nac_data' in kwargs:
        print ('use born charges')
        nac_data = kwargs.pop('nac_data')
        primitive = phonon.get_primitive()
        nac_parameters = nac_data.get_born_parameters_phonopy(primitive_cell=primitive.get_cell())
        phonon.set_nac_params(nac_parameters)

    # Normalization factor primitive to unit cell
    normalization_factor = phonon.unitcell.get_number_of_atoms()/phonon.primitive.get_number_of_atoms()

    # DOS
    phonon.set_mesh(ph_settings.dict.mesh, is_eigenvectors=True, is_mesh_symmetry=False)
    phonon.set_total_DOS(tetrahedron_method=True)
    phonon.set_partial_DOS(tetrahedron_method=True)

    total_dos = phonon.get_total_DOS()
    partial_dos = phonon.get_partial_DOS()
    dos = PhononDosData(frequencies=total_dos[0],
                        dos=total_dos[1]*normalization_factor,
                        partial_dos=np.array(partial_dos[1])*normalization_factor,
                        atom_labels=np.array(phonon.primitive.get_chemical_symbols()))

    # THERMAL PROPERTIES (per primtive cell)
    phonon.set_thermal_properties()
    t, free_energy, entropy, cv = phonon.get_thermal_properties()

    # Stores thermal properties (per unit cell) data in DB as a workflow result
    thermal_properties = ArrayData()
    thermal_properties.set_array('temperature', t)
    thermal_properties.set_array('free_energy', free_energy * normalization_factor)
    thermal_properties.set_array('entropy', entropy * normalization_factor)
    thermal_properties.set_array('heat_capacity', cv * normalization_factor)

    # BAND STRUCTURE
    phonon.set_band_structure(bands.get_bands())
    band_structure = BandStructureData(bands=bands.get_bands(),
                                       labels=bands.get_labels(),
                                       unitcell=bands.get_unitcell())

    band_structure.set_band_structure_phonopy(phonon.get_band_structure())
    return {'thermal_properties': thermal_properties,
            'dos': dos,
            'band_structure': band_structure,
            'force_constants': force_constants}


class PhononPhonopy(WorkChain):
    """
    Workchain to do a phonon calculation using phonopy

    :param structure: StructureData object that contains the crystal structure unit cell
    :param ph_settings: ParametersData object that contains a dictionary with the data needed to run phonopy:
                                  'supercell': [[2,0,0],
                                                [0,2,0],
                                                [0,0,2]],
                                  'primitive': [[1.0, 0.0, 0.0],
                                                [0.0, 1.0, 0.0],
                                                [0.0, 0.0, 1.0]],
                                  'distance': 0.01,
                                  'mesh': [40, 40, 40],
                                  # 'code': 'phonopy@boston'  # include this to run phonopy remotely otherwise run phonopy localy

    :param es_settings: ParametersData object that contains a dictionary with the setting needed to calculate the electronic structure.
                        The structure of this dictionary strongly depends on the software (VASP, QE, LAMMPS, ...)
    :param optimize: Set true to perform a crystal structure optimization before the phonon calculation (default: True)
    :param pressure: Set the external pressure (stress tensor) at which the optimization is performed in KBar (default: 0)
    """
    @classmethod
    def define(cls, spec):
        super(PhononPhonopy, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("ph_settings", valid_type=ParameterData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional arguments
        spec.input("optimize", valid_type=Bool, required=False, default=Bool(True))
        spec.input("pressure", valid_type=Float, required=False, default=Float(0.0))
        spec.input("use_nac", valid_type=Bool, required=False, default=Bool(False))

        spec.outline(_If(cls.use_optimize)(cls.optimize),
                     cls.create_displacement_calculations,
                     cls.calculate_phonon_properties,
                     cls.collect_data)

    def use_optimize(self):
        print('start phonon (pk={})'.format(self.pid))
        return self.inputs.optimize

    def optimize(self):
        print ('start optimize')
        future = submit(OptimizeStructure,
                        structure=self.inputs.structure,
                        es_settings=self.inputs.es_settings,
                        pressure=self.inputs.pressure,
                        )
        if __testing__:
            self.ctx._content['optimize'] = load_node(9357)
            return

        print ('optimize workchain: {}'.format(future.pid))

        return ToContext(optimized=future)

    def create_displacement_calculations(self):
        self.report('create displacements')

        if 'optimized' in self.ctx:
            self.ctx.final_structure = self.ctx.optimized.out.optimized_structure
            self.out('optimized_data', self.ctx.optimized.out.optimized_structure_data)
        else:
            self.ctx.final_structure = self.inputs.structure

        self.ctx.primitive_structure = get_primitive(self.ctx.final_structure, self.inputs.ph_settings)['primitive_structure']

        supercells = create_supercells_with_displacements_using_phonopy(self.ctx.final_structure,
                                                                        self.inputs.ph_settings)

        self.ctx.data_sets = supercells.pop('data_sets')
        self.ctx.number_of_displacements = len(supercells)

        calcs = {}

        # Load data from nodes
        if __testing__:
            from aiida.orm import load_node
            nodes = [96661, 96664]  # VASP
            labels = ['structure_1', 'structure_0']
            for pk, label in zip(nodes, labels):
                future = load_node(pk)
                self.ctx._content[label] = future

            #self.ctx._content['single_point'] = load_node(96667)
            return

        # Forces
        for label, supercell in supercells.iteritems():

            JobCalculation, calculation_input = generate_inputs(supercell,
                                                                # self.inputs.machine,
                                                                self.inputs.es_settings,
                                                                # pressure=self.input.pressure,
                                                                type='forces')

            calculation_input._label = label
            future = submit(JobCalculation, **calculation_input)
            self.report('{} pk = {}'.format(label, future.pid))

            calcs[label] = future

        # Born charges
        if bool(self.inputs.use_nac):
            self.report('calculate born charges')
            JobCalculation, calculation_input = generate_inputs(self.ctx.primitive_structure,
                                                                self.inputs.es_settings,
                                                                # pressure=self.input.pressure,
                                                                type='born_charges')
            future = submit(JobCalculation, **calculation_input)
            print ('single_point: {}'.format(future.pid))
            calcs['single_point'] = future

        return ToContext(**calcs)

    def calculate_phonon_properties(self):

        wf_inputs = {}
        for i in range(self.ctx.number_of_displacements):
            # This has to be changed to make uniform plugin interface
            try:
                wf_inputs['forces_{}'.format(i)] = self.ctx.get('structure_{}'.format(i)).out.output_trajectory
            except:
                wf_inputs['forces_{}'.format(i)] = self.ctx.get('structure_{}'.format(i)).out.output_array
        wf_inputs['data_sets'] = self.ctx.data_sets

        self.ctx.force_sets = create_forces_set(**wf_inputs)['force_sets']

        bands = get_path_using_seekpath(self.ctx.primitive_structure, band_resolution=30)
        phonopy_inputs = {'structure': self.ctx.final_structure,
                          'ph_settings': self.inputs.ph_settings,
                          'force_sets' : self.ctx.force_sets,
                          'bands': bands}

        if 'single_point' in self.ctx:
            nac_data = get_nac_from_data(born_charges=self.ctx.single_point.out.born_charges,
                                         epsilon=self.ctx.single_point.out.output_array,
                                         structure=self.ctx.single_point.inp.structure)
            phonopy_inputs.update(nac_data)
            self.out('nac_data', nac_data['nac_data'])

        if 'code' in self.inputs.ph_settings.get_dict():
            self.report('remote phonopy calculation')

            JobCalculation, calculation_input = generate_phonopy_params(**phonopy_inputs)
            future = submit(JobCalculation, **calculation_input)
            print ('phonopy calculation: {}'.format(future.pid))

            return ToContext(phonon_properties=future)
        else:
            self.report('local phonopy calculation')

            self.ctx.phonon_properties = get_properties_from_phonopy(**phonopy_inputs)

        return

    def collect_data(self):

        self.out('force_constants', self.ctx.phonon_properties.out.force_constants)
        self.out('thermal_properties', self.ctx.phonon_properties.out.thermal_properties)
        self.out('dos', self.ctx.phonon_properties.out.dos)
        self.out('band_structure', self.ctx.phonon_properties.out.band_structure)
        self.out('final_structure', self.ctx.final_structure)
        self.out('force_sets', self.ctx.force_sets)

        self.report('finish phonon')

        return

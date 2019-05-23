from aiida.engine import WorkChain, ToContext
from aiida.engine import workfunction

from aiida.plugins import Code, CalculationFactory, load_node, DataFactory, WorkflowFactory
from aiida.engine import run, submit

from aiida.orm import Str, Float, Bool, Int
from aiida.engine import _If, _While

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

ParameterData = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

OptimizeStructure = WorkflowFactory('phonopy.optimize')

__testing__ = False


@workfunction
def create_supercells_with_displacements_using_phono3py(structure, ph_settings, cutoff):
    """
    Use phono3py to create the supercells with displacements to calculate the force constants by using
    finite displacements methodology

    :param structure: StructureData object
    :param phonopy_input: ParametersData object containing a dictionary with the data needed for phonopy
    :param cutoff: FloatData object containing the value of the cutoff for 3rd order FC in Angstroms (if 0 no cutoff is applied)
    :return: A set of StructureData Objects containing the supercells with displacements
    """
    from phono3py.phonon3 import Phono3py

    from aiida_phonopy.common.utils import phonopy_atoms_from_structure

    # Generate phonopy phonon object
    phono3py = Phono3py(phonopy_atoms_from_structure(structure),
                        supercell_matrix=ph_settings.dict.supercell,
                        primitive_matrix=ph_settings.dict.primitive,
                        symprec=ph_settings.dict.symmetry_tolerance,
                        log_level=1)

    if float(cutoff) == 0:
        cutoff = None
    else:
        cutoff = float(cutoff)

    phono3py.generate_displacements(distance=ph_settings.dict.distance,
                                    cutoff_pair_distance=cutoff)

    cells_with_disp = phono3py.get_supercells_with_displacements()

    # Transform cells to StructureData and set them ready to return
    data_sets = phono3py.get_displacement_dataset()
    data_sets_object = ForceSetsData(data_sets3=data_sets)

    disp_cells = {'data_sets': data_sets_object}
    for i, phonopy_supercell in enumerate(cells_with_disp):
        if phonopy_supercell is None:
            print ('structure_{} cutoff skip'.format(i))
            continue
        supercell = StructureData(cell=phonopy_supercell.get_cell())
        for symbol, position in zip(phonopy_supercell.get_chemical_symbols(),
                                    phonopy_supercell.get_positions()):
            supercell.append_atom(position=position, symbols=symbol)
        disp_cells['structure_{}'.format(i)] = supercell

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
    force_sets = ForceSetsData(data_sets3=data_sets.get_data_sets3())

    forces = []
    for i in range(data_sets.get_number_of_displacements()):
        if 'forces_{}'.format(i) in kwargs:
            forces.append(kwargs.pop('forces_{}'.format(i)).get_array('forces')[-1])

    force_sets.set_forces(forces)

    return {'force_sets': force_sets}

@workfunction
def get_force_constants3(data_sets, structure, ph_settings):

    from phono3py.phonon3 import Phono3py
    from aiida_phonopy.common.utils import phonopy_atoms_from_structure

    # Generate phonopy phonon object
    phono3py = Phono3py(phonopy_atoms_from_structure(structure),
                        supercell_matrix=ph_settings.dict.supercell,
                        primitive_matrix=ph_settings.dict.primitive,
                        symprec=ph_settings.dict.symmetry_tolerance,
                        log_level=1)

    phono3py.produce_fc3(data_sets.get_forces3(),
                         displacement_dataset=data_sets.get_data_sets3(),
                         is_translational_symmetry=True,
                         is_permutation_symmetry=True,
                         is_permutation_symmetry_fc2=True)
    fc3 = phono3py.get_fc3()
    fc2 = phono3py.get_fc2()

    force_constants_2 = ForceConstantsData(data=fc2)
    force_constants_3 = ForceConstantsData(data=fc3)

    return {'force_constants_2order': force_constants_2,
            'force_constants_3order': force_constants_3}


def get_recover_calc(wc, label):

    if type(wc) is Bool:
        return None

    try:
        pk_list = wc.out.calcs_pk
    except:
        pk_list = []

    out_list = [load_node(pk) for pk in pk_list]

    for calc in wc.get_outputs() + out_list:
        if calc.label == label:
            if calc.get_state() == 'FINISHED':
                return calc

    return None


def cut_supercells(supercells, data_sets):

    if type(data_sets) is Bool:
        return supercells

    data_sets_dict = data_sets.get_data_sets3()

    if data_sets_dict is None:
        data_sets_dict = data_sets.get_data_sets()
        for i in range(data_sets.get_number_of_displacements()):
            del supercells['structure_{}'.format(i)]
            print ('removed: structure_{}'.format(i))
    else:
        calculation_list = []
        i = len(data_sets_dict['first_atoms'])
        for first_atoms in data_sets_dict['first_atoms']:
            for second_atoms in first_atoms['second_atoms']:
                if 'included' in second_atoms:
                    if not second_atoms['included']:
                        calculation_list.append('structure_{}'.format(i))
                i += 1

        supercells = {key:value for key, value in supercells.items() if key in calculation_list}
    return supercells


def get_forces_from_sets(data_sets, index):

    if type(data_sets) is Bool:
        return None

    forces = data_sets.get_forces3()
    data_sets_dict = data_sets.get_data_sets3()

    if data_sets_dict is None:
        harmonic_displacements = data_sets.get_number_of_displacements()
    else:
        harmonic_displacements = len(data_sets_dict['first_atoms'])

    if index < harmonic_displacements:
        return forces[index]

    if data_sets_dict is not None:
        i = i_ref = len(data_sets_dict['first_atoms'])
        #i_ref = len(data_sets['first_atoms'])
        for first_atoms in data_sets_dict['first_atoms']:
            for second_atoms in first_atoms['second_atoms']:
                if 'included' in second_atoms:
                    if second_atoms['included']:
                        if index == i_ref:
                            return forces[i]
                        i += 1
                else:
                    if index == i_ref:
                        return forces[i]
                    i += 1
                i_ref += 1

    return None



class PhononPhono3py(WorkChain):
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
    :param optimize: (optional) Set true to perform a crystal structure optimization before the phonon calculation (default: True)
    :param pressure: (optional) Set the external pressure (stress tensor) at which the optimization is performed in KBar (default: 0)
    :param use_nac: (optional) Bool type object that determines if non-analytical correction term is calculated or not
    :param calculate_fc: (optional) Bool type object that determines if non-analytical correction term is calculated or not
    :param chunks: (optional) Int type object that defines the  maximum number of maximum calculations allowed to run simultaneously (default:100)
    :param cutoff: (optional) Float type object that defined the distance cutoff used in the generation of supercells with displacements in phono3py.
    :param data_sets: (optional) ForceSetsData type object that contains the forces data from a previous calculation to be reused

    """
    @classmethod
    def define(cls, spec):
        super(PhononPhono3py, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("ph_settings", valid_type=ParameterData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional arguments
        spec.input("optimize", valid_type=Bool, required=False, default=Bool(True))
        spec.input("pressure", valid_type=Float, required=False, default=Float(0.0))
        spec.input("use_nac", valid_type=Bool, required=False, default=Bool(False))  # false by default
        spec.input("calculate_fc", valid_type=Bool, required=False, default=Bool(False))  # false by default
        spec.input("chunks", valid_type=Int, required=False, default=Int(100))
        spec.input("cutoff", valid_type=Float, required=False, default=Float(0))
        spec.input("recover", required=False, default=Bool(False)) # temporal patch
        spec.input("data_sets", required=False, default=Bool(False))

        spec.outline(_If(cls.use_optimize)(cls.optimize),
                     # cls.create_displacement_calculations,
                     _While(cls.continue_submitting)(cls.create_displacement_calculations_chunk),
                     cls.collect_data,
                     _If(cls.calculate_fc)(cls.calculate_force_constants))
        # spec.outline(cls.calculate_force_constants)  # testing

    def use_optimize(self):
        print('start phonon3 (pk={})'.format(self.pid))
        return self.inputs.optimize

    def calculate_fc(self):
        return self.inputs.calculate_fc

    def continue_submitting(self):

        if 'i_disp' in self.ctx:
            if self.ctx.i_disp < 1:
                return False
            self.ctx.i_disp -= 1
        return True

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

        from aiida_phonopy.workchains.phonon import get_primitive

        print ('create displacements')
        self.report('create displacements')

        if 'optimized' in self.ctx:
            self.ctx.final_structure = self.ctx.optimized.out.optimized_structure
            self.out('optimized_data', self.ctx.optimized.out.optimized_structure_data)
        else:
            self.ctx.final_structure = self.inputs.structure

        self.ctx.primitive_structure = get_primitive(self.ctx.final_structure,
                                                     self.inputs.ph_settings)['primitive_structure']

        supercells = create_supercells_with_displacements_using_phono3py(self.ctx.final_structure,
                                                                         self.inputs.ph_settings,
                                                                         self.inputs.cutoff)

        self.ctx.data_sets = supercells.pop('data_sets')
        self.ctx.number_of_displacements = len(supercells)

        if __testing__:
            f = open('labels', 'r')
            lines = f.readlines()
            f.close()

            from aiida.orm import load_node
            nodes = [int(line.split()[3]) for line in lines]
            print (nodes)
            labels = [line.split()[0] for line in lines]
            print (labels)
            for pk, label in zip(nodes, labels):
                future = load_node(pk)
                self.ctx._content[label] = future
                print ('{} pk = {}'.format(label, pk))

            return

        calcs = {}
        for label, supercell in supercells.iteritems():

            JobCalculation, calculation_input = generate_inputs(supercell,
                                                                # self.inputs.machine,
                                                                self.inputs.es_settings,
                                                                # pressure=self.input.pressure,
                                                                type='forces')

            calculation_input._label = label
            future = submit(JobCalculation, **calculation_input)
            print ('{} pk = {}'.format(label, future.pid))
            # self.report('{} pk = {}'.format(label, future.pid))

            calcs[label] = future

        # Born charges (for primitive cell)
        if bool(self.inputs.use_nac):
            self.report('calculate born charges')
            JobCalculation, calculation_input = generate_inputs(self.ctx.primitive_structure,
                                                                # self.inputs.machine,
                                                                self.inputs.es_settings,
                                                                # pressure=self.input.pressure,
                                                                type='born_charges')
            future = submit(JobCalculation, **calculation_input)
            print ('single_point: {}'.format(future.pid))
            calcs['single_point'] = future

        return ToContext(**calcs)

    def create_displacement_calculations_chunk(self):

        from aiida_phonopy.workchains.phonon import get_primitive

        if 'optimized' in self.ctx:
            self.ctx.final_structure = self.ctx.optimized.out.optimized_structure
            self.out('optimized_data', self.ctx.optimized.out.optimized_structure_data)
        else:
            self.ctx.final_structure = self.inputs.structure

        self.ctx.primitive_structure = get_primitive(self.ctx.final_structure,
                                                     self.inputs.ph_settings)['primitive_structure']

        supercells = create_supercells_with_displacements_using_phono3py(self.ctx.final_structure,
                                                                         self.inputs.ph_settings,
                                                                         self.inputs.cutoff)


        self.ctx.data_sets = supercells.pop('data_sets')
        self.ctx.number_of_displacements = len(supercells)

        supercells = cut_supercells(supercells, self.inputs.data_sets)

        calcs = {}

        n_disp = len(supercells)
        if 'i_disp' in self.ctx:
            list = range(self.ctx.i_disp * int(self.inputs.chunks),
                         (self.ctx.i_disp + 1) * int(self.inputs.chunks))
        else:
            self.ctx.i_disp = n_disp / int(self.inputs.chunks)
            list = range(self.ctx.i_disp * int(self.inputs.chunks), n_disp)
            print ('create displacements')
            self.report('create displacements')
            print ('total displacements: {}'.format(n_disp))

            # Born charges (for primitive cell)
            if bool(self.inputs.use_nac):
                self.report('calculate born charges')
                JobCalculation, calculation_input = generate_inputs(self.ctx.primitive_structure,
                                                                    # self.inputs.machine,
                                                                    self.inputs.es_settings,
                                                                    # pressure=self.input.pressure,
                                                                    type='born_charges')
                future = submit(JobCalculation, **calculation_input)
                print ('single_point: {}'.format(future.pid))
                calcs['single_point'] = future

        supercell_list = np.array(supercells.items())[list]

        for label, supercell in supercell_list:

            recover_calc = get_recover_calc(self.inputs.recover, label)

            # Recover calculations (temporal beta patch)
            if recover_calc is not None:
                self.ctx._content[label] = recover_calc
                print ('recovered: {}'.format(label))
                continue
            # Recover calculations end (temporal beta patch)

            JobCalculation, calculation_input = generate_inputs(supercell,
                                                                # self.inputs.machine,
                                                                self.inputs.es_settings,
                                                                # pressure=self.input.pressure,
                                                                type='forces')

            calculation_input['_label'] = label
            future = submit(JobCalculation, **calculation_input)
            print ('{} pk = {}'.format(label, future.pid))

            calcs[label] = future

        return ToContext(**calcs)

    def collect_data(self):

        from aiida_phonopy.workchains.phonon import get_nac_from_data
        self.report('collect data and create force_sets')

        wf_inputs = {}
        for i in range(self.ctx.data_sets.get_number_of_displacements()):

            forces = get_forces_from_sets(self.inputs.data_sets, i)
            if forces is not None:
                print ('read_{}'.format(i))
                array_data = ArrayData()
                array_data.set_array('forces', np.array([forces]))
                wf_inputs['forces_{}'.format(i)] = array_data
                continue

            calc = self.ctx.get('structure_{}'.format(i))
            #print ('collect structure_{}'.format(i))
            if calc is not None:
                print ('calculated_{} OK'.format(i))

                # This has to be changed to make uniform plugin interface
                try:
                    wf_inputs['forces_{}'.format(i)] = calc.out.output_trajectory
                except:
                    wf_inputs['forces_{}'.format(i)] = calc.out.output_array

        wf_inputs['data_sets'] = self.ctx.data_sets

        self.ctx.force_sets = create_forces_set(**wf_inputs)['force_sets']

        if 'single_point' in self.ctx:
            nac_data = get_nac_from_data(born_charges=self.ctx.single_point.out.born_charges,
                                         epsilon=self.ctx.single_point.out.output_array,
                                         structure=self.ctx.primitive_structure)

            self.out('nac_data', nac_data['nac_data'])

        self.out('force_sets', self.ctx.force_sets)
        self.out('final_structure', self.ctx.final_structure)

        self.report('phonon3py calculation finished ')

    def calculate_force_constants(self):

        force_constants = get_force_constants3(self.ctx.force_sets,
                                               self.ctx.final_structure,
                                               self.inputs.ph_settings)

        self.out('force_constants_2order', force_constants['force_constants_2order'])
        self.out('force_constants_3order', force_constants['force_constants_3order'])

        return

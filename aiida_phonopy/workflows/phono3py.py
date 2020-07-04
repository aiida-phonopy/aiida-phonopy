import numpy as np
from phonopy import Phonopy
from aiida.engine import WorkChain
from aiida.plugins import WorkflowFactory, DataFactory
from aiida.orm import Float, Bool, Str, Code
from aiida.engine import if_
from aiida_phonopy.common.generate_inputs import (get_calcjob_builder,
                                                  get_immigrant_builder)
from aiida_phonopy.common.utils import (
    get_force_sets, get_force_constants, get_nac_params, get_phonon,
    get_phonon_setting_info, check_imported_supercell_structure,
    from_node_id_to_aiida_node_id, get_data_from_node_id)


PhonopyWorkChain = WorkflowFactory('phonopy.phonopy')
Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

@workfunction
def create_supercells_with_displacements_using_phono3py(structure, ph_settings, cutoff):
    """
    Use phono3py to create the supercells with displacements to calculate the force constants by using
    finite displacements methodology

    :param structure: StructureData object
    :param phonopy_input: Dict object containing a dictionary with the data needed for phonopy
    :param cutoff: FloatData object containing the value of the cutoff for 3rd order FC in Angstroms (if 0 no cutoff is applied)
    :return: A set of StructureData Objects containing the supercells with displacements
    """
    from phono3py import Phono3py

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


class Phono3pyWorkChain(WorkChain):
    """Phono3py workchain


    """
    @classmethod
    def define(cls, spec):
        super(Phono3pyWorkChain, cls).define(spec)
        spec.expose_inputs(PhonopyWorkChain,
                           exclude=['immigrant_calculation_folders',
                                    'calculation_nodes', 'dry_run'])

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

        from aiida_phonopy.workflows.phonopy import get_primitive

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

        from aiida_phonopy.workflows.phonopy import get_primitive

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

        from aiida_phonopy.workflows.phonopy import get_nac_from_data
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

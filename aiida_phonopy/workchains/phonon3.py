# Works run by the daemon (using submit)

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction

from aiida.orm import Code, CalculationFactory, load_node, DataFactory, WorkflowFactory
from aiida.work.run import run, submit, async

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


@workfunction
def create_supercells_with_displacements_using_phono3py(structure, ph_settings):
    """
    Use phono3py to create the supercells with displacements to calculate the force constants by using
    finite displacements methodology

    :param structure: StructureData object
    :param phonopy_input: ParametersData object containing a dictionary with the data needed for phonopy
    :return: A set of StructureData Objects containing the supercells with displacements
    """
    from phono3py.phonon3 import Phono3py

    from aiida_phonopy.workchains.phonon import phonopy_bulk_from_structure

    # Generate phonopy phonon object
    phono3py = Phono3py(phonopy_bulk_from_structure(structure),
                        supercell_matrix=ph_settings.dict.supercell,
                        primitive_matrix=ph_settings.dict.primitive,
                        symprec=ph_settings.dict.symmetry_precision,
                        log_level=1)

    phono3py.generate_displacements(distance=ph_settings.dict.distance)

    cells_with_disp = phono3py.get_supercells_with_displacements()

    # Transform cells to StructureData and set them ready to return
    data_sets = phono3py.get_displacement_dataset()
    print data_sets
    data_sets_object = ForceSetsData(data_sets3=data_sets)

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


def get_force_constants3(data_sets, structure, ph_settings):

    print ('create 3rd order force constants')
    from phono3py.phonon3 import Phono3py

    from aiida_phonopy.workchains.phonon import phonopy_bulk_from_structure

    # Generate phonopy phonon object
    phono3py = Phono3py(phonopy_bulk_from_structure(structure),
                        supercell_matrix=ph_settings.dict.supercell,
                        primitive_matrix=ph_settings.dict.primitive,
                        symprec=ph_settings.dict.symmetry_precision,
                        log_level=1)

    phono3py.generate_displacements(distance=ph_settings.dict.distance)


    phono3py.produce_fc3(data_sets.get_forces3(),
                         displacement_dataset=data_sets.get_data_sets3(),
                         is_translational_symmetry=True,
                         is_permutation_symmetry=True,
                         is_permutation_symmetry_fc2=True)
    fc3 = phono3py.get_fc3()
    fc2 = phono3py.get_fc2()

    return fc3, fc2




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
    :param optimize: Set true to perform a crystal structure optimization before the phonon calculation (default: True)
    :param pressure: Set the external pressure (stress tensor) at which the optimization is performed in KBar (default: 0)
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
        spec.input("use_nac", valid_type=Bool, required=False, default=Bool(True))

        spec.outline(_If(cls.use_optimize)(cls.optimize),
                     cls.create_displacement_calculations, cls.collect_data)

    def use_optimize(self):
        print('start phonon3 (pk={})'.format(self.pid))
        return self.inputs.optimize

    def optimize(self):
        print ('start optimize')
        future = submit(OptimizeStructure,
                        structure=self.inputs.structure,
                        es_settings=self.inputs.es_settings,
                        pressure=self.inputs.pressure,
                        )
        # For testing
        testing = False
        if testing:
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
                                                                         self.inputs.ph_settings)

        self.ctx.data_sets = supercells.pop('data_sets')
        self.ctx.number_of_displacements = len(supercells)

        calcs = {}

        # Load data from nodes
        testing = False
        if testing:
            from aiida.orm import load_node
            nodes = [9378, 9381]  # VASP
            labels = ['structure_1', 'structure_0']
            for pk, label in zip(nodes, labels):
                future = load_node(pk)
                self.ctx._content[label] = future

            self.ctx._content['single_point'] = load_node(9385)
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
            # print label, future.pid
            self.report('{} pk = {}'.format(label, future.pid))

            calcs[label] = future

        # Born charges
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

    def collect_data(self):

        from aiida_phonopy.workchains.phonon import get_nac_from_data
        print ('calculate force constants')
        self.report('calculate force constants')

        wf_inputs = {}
        for i in range(self.ctx.number_of_displacements):
            # This has to be changed to make uniform plugin interface
            try:
                wf_inputs['forces_{}'.format(i)] = self.ctx.get('structure_{}'.format(i)).out.output_trajectory
            except:
                wf_inputs['forces_{}'.format(i)] = self.ctx.get('structure_{}'.format(i)).out.output_array

        wf_inputs['data_sets'] = self.ctx.data_sets

        self.ctx.force_sets = create_forces_set(**wf_inputs)['force_sets']

        if 'single_point' in self.ctx:
            nac_data = get_nac_from_data(born_charges=self.ctx.single_point.out.born_charges,
                                         epsilon=self.ctx.single_point.out.output_array,
                                         structure=self.ctx.primitive_structure)

            self.out('nac_data', nac_data['nac_data'])

        self.out('force_sets', self.ctx.force_sets)

        self.report('phonon3py calculation finished ')

        # test

        force_constants = get_force_constants3(self.ctx.force_sets,
                                               self.ctx.final_structure,
                                               self.inputs.ph_settings)

        print (force_constants)

        return

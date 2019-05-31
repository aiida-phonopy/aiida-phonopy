from aiida.engine import WorkChain, ToContext

from aiida.common.extendeddicts import AttributeDict
from aiida.engine import workfunction
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.orm import Float, Bool, Int, Str, load_node, Code
from aiida.engine import if_
import numpy as np
from aiida_phonopy.common.generate_inputs import generate_inputs
from aiida_phonopy.common.utils import (phonopy_atoms_from_structure,
                                        phonopy_atoms_to_structure,
                                        get_forces_from_uuid,
                                        get_born_epsilon_from_uuid,
                                        get_force_sets,
                                        get_force_constants,
                                        get_nac_params,
                                        get_path_using_seekpath,
                                        get_phonon)

# Should be improved by some kind of WorkChainFactory
# For now all workchains should be copied to aiida/workflows
# from aiida.workflows.wc_optimize import OptimizeStructure

Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')
OptimizeStructure = WorkflowFactory('phonopy.optimize')
BandStructureData = DataFactory('phonopy.band_structure')
PhononDosData = DataFactory('phonopy.phonon_dos')


class PhononPhonopy(WorkChain):
    """ Workchain to do a phonon calculation using phonopy

    :param structure: StructureData object that contains the crystal structure
        unit cell
    :param cell_matrices: ArrayData object having supercell matrix and
        optionally primitive matrix. (default: primitive_matrix = unit matrix)
            'supercell_matrix': [[2, 0, 0],
                                 [0, 2, 0],
                                 [0, 0, 2]],
            'primitive_matrix': [[1.0, 0.0, 0.0],
                                 [0.0, 1.0, 0.0],
                                 [0.0, 0.0, 1.0]]
    :param calculator_settings: Dict object that contains a
        dictionary with the setting needed to calculate the electronic
        structure:
            {'forces': force_config,
             'nac': nac_config}
        where force_config and nac_config are used for the supercell force
        calculation and Born effective charges and dielectric constant
        calculation in primitive cell, respectively.
    :ph_settings: Dict object. Needed to run phonon calculation.
        {'mesh': [20, 20, 20]}. Optional. (default: {mesh: [20, 20, 20]})
    :param optimize: Set true to perform a crystal structure optimization
        before the phonon calculation (default: True)
    :param pressure: Set the external pressure (stress tensor) at which the
        optimization is performed in KBar (default: 0)
    :is_nac: Bool object. Whether running non-analytical term correction.
        Optional. (default: False)
    :run_phonopy: Bool object. Whether running phonon calculation or not.
        Optional. (default: False)
    :remote_phonopy: Bool object. Whether running phonon calculation or not.
        Optional. (default: False)
    :code_string: Str object. Needed to run phonopy remotely. Optional.
    :options: Dict object. Needed to run phonopy remotely. Optional.
    :distance: Float object. Displacement distance. Optional. (default: 0.01)
    :symmetry_tolerance: Float object. Symmetry tolerance. Optional.
        (default: 1e-5)

    """
    @classmethod
    def define(cls, spec):
        super(PhononPhonopy, cls).define(spec)
        spec.input('structure', valid_type=StructureData, required=True)
        spec.input('cell_matrices', valid_type=ArrayData, required=True)
        # Optional arguments
        spec.input('calculation_uuids',
                   valid_type=Dict, required=False)
        spec.input('calculator_settings',
                   valid_type=Dict, required=False)
        spec.input('displacement_dataset',
                   valid_type=Dict, required=False)
        spec.input('code_string', valid_type=Str, required=False)
        spec.input('options', valid_type=Dict, required=False)
        spec.input('phonon_settings', valid_type=Dict, required=False)
        spec.input('is_nac',
                   valid_type=Bool, required=False, default=Bool(False))
        spec.input('optimize',
                   valid_type=Bool, required=False, default=Bool(True))
        spec.input('pressure',
                   valid_type=Float, required=False, default=Float(0.0))
        spec.input('distance',
                   valid_type=Float, required=False, default=Float(0.01))
        spec.input('symmetry_tolerance',
                   valid_type=Float, required=False, default=Float(1e-5))
        spec.input('run_phonopy',
                   valid_type=Bool, required=False, default=Bool(False))
        spec.input('remote_phonopy',
                   valid_type=Bool, required=False, default=Bool(False))

        spec.outline(cls.check_settings,
                     if_(cls.use_optimize)(cls.optimize),
                     cls.set_unitcell,
                     cls.create_displacements,
                     if_(cls.import_nodes)(
                         cls.read_force_and_nac_calculations,
                     ).else_(
                         cls.run_force_and_nac_calculations,
                     ),
                     cls.create_force_sets,
                     cls.create_force_constants,
                     if_(cls.is_nac)(cls.create_nac_params),
                     cls.set_default_outputs,
                     cls.prepare_phonopy_inputs,
                     if_(cls.run_phonopy)(
                         if_(cls.remote_phonopy)(
                             cls.run_phonopy_remote,
                             cls.collect_data,
                         ).else_(
                             cls.run_phonopy_in_workchain,
                         )))

        spec.output('force_constants', valid_type=ArrayData, required=False)
        spec.output('final_structure', valid_type=StructureData, required=True)
        spec.output('force_sets', valid_type=ArrayData, required=False)
        spec.output('displacement_dataset', valid_type=Dict, required=False)
        spec.output('nac_params', valid_type=ArrayData, required=False)
        spec.output('thermal_properties', valid_type=ArrayData, required=False)
        spec.output('band_structure', valid_type=BandStructureData, required=False)
        spec.output('dos', valid_type=PhononDosData, required=False)

    def use_optimize(self):
        return self.inputs.optimize

    def remote_phonopy(self):
        return self.inputs.remote_phonopy

    def run_phonopy(self):
        return self.inputs.run_phonopy

    def is_nac(self):
        return self.inputs.is_nac

    def import_nodes(self):
        return ('calculation_uuids' in self.inputs and
                'displacement_dataset' in self.inputs)

    def check_settings(self):
        self.report('check phonon_settings')

        if 'phonon_settings' in self.inputs:
            phonon_settings_dict = self.inputs.phonon_settings.get_dict()
        else:
            phonon_settings_dict = {}

        cell_keys = self.inputs.cell_matrices.get_arraynames()
        if 'supercell_matrix' not in cell_keys:
            raise RuntimeError(
                "supercell_matrix is not found in cell_matrices.")
        else:
            dim = self.inputs.cell_matrices.get_array('supercell_matrix')
            if len(dim.ravel()) == 3:
                smat = np.diag(dim)
            else:
                smat = np.array(dim)
            if not np.issubdtype(smat.dtype, np.integer):
                raise TypeError("supercell_matrix is not integer matrix.")
            else:
                phonon_settings_dict['supercell_matrix'] = smat.tolist()
        if 'primitive_matrix' in cell_keys:
            pmat = self.inputs.cell_matrices.get_array('primitive_matrix')
            pmat = pmat.astype(float).tolist()
            phonon_settings_dict['primitive_matrix'] = pmat

        if 'mesh' not in phonon_settings_dict:
            phonon_settings_dict['mesh'] = [20, 20, 20]
        phonon_settings_dict['distance'] = float(self.inputs.distance)
        tolerance = float(self.inputs.symmetry_tolerance)
        phonon_settings_dict['symmetry_tolerance'] = tolerance

        if self.inputs.run_phonopy and self.inputs.remote_phonopy:
            if ('code_string' not in self.inputs or
                'options' not in self.inputs):
                raise RuntimeError(
                    "code_string and options have to be specified.")

        self.ctx.phonon_settings = Dict(dict=phonon_settings_dict)
        self.ctx.phonon_settings.label = 'phonon_settings'

        self.report("Phonon settings: %s" % phonon_settings_dict)

    def optimize(self):
        self.report('start optimize')
        future = self.submit(
            OptimizeStructure,
            structure=self.inputs.structure,
            calculator_settings=self.inputs.calculator_settings,
            pressure=self.inputs.pressure)
        self.report('optimize workchain: {}'.format(future.pk))

        return ToContext(optimized=future)

    def set_unitcell(self):
        self.report('set unit cell')

        if 'optimized' in self.ctx:
            self.ctx.final_structure = self.ctx.optimized.out.optimized_structure
            self.out('optimized_data',
                     self.ctx.optimized.out.optimized_structure_data)
        else:
            self.ctx.final_structure = self.inputs.structure

    def create_displacements(self):
        """Create supercells with displacements and primitive cell

        Use phonopy to create the supercells with displacements to
        calculate the force constants by using finite displacements
        methodology

        """

        self.report('create cells for phonon calculation')

        from phonopy import Phonopy

        phonon_settings_dict = self.ctx.phonon_settings.get_dict()
        if 'primitive_matrix' in phonon_settings_dict:
            phonon = Phonopy(
                phonopy_atoms_from_structure(self.ctx.final_structure),
                phonon_settings_dict['supercell_matrix'],
                primitive_matrix=phonon_settings_dict['primitive_matrix'],
                symprec=phonon_settings_dict['symmetry_tolerance'])
        else:
            phonon = Phonopy(
                phonopy_atoms_from_structure(self.ctx.final_structure),
                phonon_settings_dict['supercell_matrix'],
                symprec=phonon_settings_dict['symmetry_tolerance'])

        if 'displacement_dataset' in self.inputs:
            phonon.dataset = self.inputs.displacement_dataset.get_dict()
        else:
            phonon.generate_displacements(
                distance=phonon_settings_dict['distance'])

        self.ctx.disp_dataset = {
            'dataset': Dict(dict=phonon.dataset),
            'supercells': []}
        self.ctx.disp_dataset['dataset'].label = 'displacement_dataset'

        self.ctx.primitive_structure = phonopy_atoms_to_structure(
            phonon.primitive)

        cells_with_disp = phonon.supercells_with_displacements
        for i, phonopy_supercell in enumerate(cells_with_disp):
            supercell = phonopy_atoms_to_structure(phonopy_supercell)
            self.ctx.disp_dataset['supercells'].append(supercell)

    def run_force_and_nac_calculations(self):
        self.report('run force calculations')

        # Forces
        for i, supercell in enumerate(self.ctx.disp_dataset['supercells']):
            label = "supercell_%03d" % (i + 1)
            builder = generate_inputs(supercell,
                                      self.inputs.calculator_settings,
                                      calc_type='forces',
                                      label=label)
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

            # append_ can't be used because the order can change.
            # self.to_context(supercell_calcs=append_(future))

        # Born charges and dielectric constant
        if self.inputs.is_nac:
            self.report('calculate born charges and dielectric constant')
            builder = generate_inputs(self.ctx.primitive_structure,
                                      self.inputs.calculator_settings,
                                      calc_type='nac',
                                      label='born_and_epsilon')
            future = self.submit(builder)
            self.report('born_and_epsilon: {}'.format(future.pk))
            self.to_context(**{'born_and_epsilon': future})

    def read_force_and_nac_calculations(self):
        self.report('read calculation nodes')

        uuids = self.inputs.calculation_uuids.get_dict()['uuids']
        if self.inputs.is_nac:
            uuid = uuids.pop()
            label = 'born_and_epsilon'
            n = load_node(uuid)
            self.ctx[label] = {'output_born_charges': n.outputs.output_born_charges,
                               'output_dielectrics': n.outputs.output_dielectrics}
        for i, uuid in enumerate(uuids):
            label = "supercell_%03d" % (i + 1)
            n = load_node(uuid)
            self.ctx[label] = {'output_forces': n.outputs.output_forces}

    def create_force_sets(self):
        """Build datasets from forces of supercells with displacments"""

        self.report('create force sets')

        # VASP specific
        forces_dict = {}
        for i in range(len(self.ctx.disp_dataset['supercells'])):
            label = "supercell_%03d" % (i + 1)
            calc = self.ctx[label]
            if type(calc) is dict:
                out_dict = calc
            else:
                out_dict = calc.outputs
            if ('output_forces' in out_dict and
                'final' in out_dict['output_forces'].get_arraynames()):
                label = "forces_%03d" % (i + 1)
                forces_dict[label] = out_dict['output_forces']

        if len(forces_dict) != len(self.ctx.disp_dataset['supercells']):
            raise RuntimeError("Forces could not be retrieved.")

        self.ctx.force_sets = get_force_sets(**forces_dict)

    def create_nac_params(self):
        self.report('create nac data')

        # VASP specific
        # Call workfunction to make links
        if type(self.ctx.born_and_epsilon) is dict:
            out_dict = self.ctx.born_and_epsilon
        else:
            out_dict = self.ctx.born_and_epsilon.outputs
        self.ctx.nac_params = get_nac_params(
            born_charges=out_dict['output_born_charges'],
            epsilon=out_dict['output_dielectrics'])

    def create_force_constants(self):
        self.report('create force constants')

        self.ctx.force_constants = get_force_constants(
            self.ctx.final_structure,
            self.ctx.phonon_settings,
            self.ctx.force_sets,
            self.ctx.disp_dataset['dataset'])

    def set_default_outputs(self):
        self.out('force_constants', self.ctx.force_constants)
        self.out('final_structure', self.ctx.final_structure)
        self.out('force_sets', self.ctx.force_sets)
        self.out('displacement_dataset', self.ctx.disp_dataset['dataset'])
        if 'nac_params' in self.ctx:
            self.out('nac_params', self.ctx.nac_params)

    def prepare_phonopy_inputs(self):
        self.ctx.band_paths = get_path_using_seekpath(
            self.ctx.primitive_structure,
            band_resolution=Int(30))

    def run_phonopy_remote(self):
        """Run phonopy at remote computer"""

        self.report('remote phonopy calculation')

        code_string = self.inputs.code_string.value
        builder = Code.get_from_string(code_string).get_builder()
        builder.structure = self.ctx.final_structure
        builder.parameters = self.ctx.phonon_settings
        builder.metadata.options.update(self.inputs.options)
        builder.metadata.options.parser_name = 'phonopy'
        builder.metadata.options.input_filename = 'phonopy.conf'
        builder.metadata.options.output_filename = 'phonopy.stdout'
        builder.metadata.label = self.inputs.metadata.label

        if 'force_constants' in self.ctx:
            builder.force_constants = self.ctx.force_constants
        else:
            builder.force_sets = self.ctx.force_sets
            builder.displacement_dataset = self.ctx.disp_dataset['dataset']
        if 'nac_params' in self.ctx:
            builder.nac_params = self.ctx.nac_params
        if 'band_paths' in self.ctx:
            builder.bands = self.ctx.band_paths

        future = self.submit(builder)

        self.report('phonopy calculation: {}'.format(future.pk))
        self.to_context(**{'phonon_properties': future})
        # return ToContext(phonon_properties=future)

    def collect_data(self):
        self.report('collect data')
        self.out('thermal_properties',
                 self.ctx.phonon_properties.outputs.thermal_properties)
        self.out('dos', self.ctx.phonon_properties.outputs.dos)
        self.out('band_structure',
                 self.ctx.phonon_properties.outputs.band_structure)

        self.report('finish phonon')

    def run_phonopy_in_workchain(self):
        self.report('phonopy calculation in workchain')

        params = {}
        if 'nac_params' in self.ctx:
            params['nac_params'] = self.ctx.nac_params
        result = get_phonon(self.ctx.final_structure,
                            self.ctx.phonon_settings,
                            self.ctx.force_constants,
                            self.ctx.band_paths,
                            **params)
        self.out('thermal_properties', result['thermal_properties'])
        self.out('dos', result['dos'])
        self.out('band_structure', result['band_structure'])

        self.report('finish phonon')

from aiida.engine import WorkChain
from aiida.plugins import DataFactory, CalculationFactory
from aiida.orm import Float, Bool, Str, Code
from aiida.engine import if_, while_
from aiida_phonopy.common.builders import (
    get_calcjob_builder, get_force_calcjob_inputs, get_immigrant_builder)
from aiida_phonopy.common.utils import (
    get_force_constants, get_nac_params, get_phonon,
    generate_phonopy_cells, compare_structures,
    from_node_id_to_aiida_node_id, get_data_from_node_id,
    get_vasp_force_sets_dict, collect_vasp_forces_and_energies)
from aiida_phonopy.workflows.nac import NacParamsWorkChain


# Should be improved by some kind of WorkChainFactory
# For now all workchains should be copied to aiida/workflows

Dict = DataFactory('dict')
ArrayData = DataFactory('array')
XyData = DataFactory('array.xy')
StructureData = DataFactory('structure')
BandsData = DataFactory('array.bands')
PhonopyCalculation = CalculationFactory('phonopy.phonopy')


class PhonopyWorkChain(WorkChain):
    """Phonopy workchain

    inputs
    ------
    structure : StructureData
        Unit cell structure.
    calculator_settings : Dict
        Settings to run force and nac calculations. For example,
            {'forces': force_config,
             'nac': nac_config}
        At least 'forces' key is necessary. 'nac' is optional.
        force_config is used for supercell force calculation. nac_config
        are used for Born effective charges and dielectric constant calculation
        in primitive cell. The primitive cell is chosen by phonopy
        automatically.
    phonon_settings : Dict
        Setting to run phonon calculation. Keys are:
        supercell_matrix : list or list of list
            Multiplicity to create supercell from unit cell. Three integer
            values (list) or 3x3 integer values (list of list).
        mesh : list of float, optional
            List of three integer values or float to represent distance between
            neighboring q-points. Default is 100.0.
        distance : float, optional
            Atomic displacement distance. Default is 0.01.
        is_nac : bool, optional
            Whether running non-analytical term correction or not. Default is
            False.
        displacement_dataset : dict
            Atomic displacement dataset that phonopy can understand.
        options : dict
            AiiDA calculation options for phonon calculation used when both of
            run_phonopy and remote_phonopy are True.
    subtract_residual_forces : Bool, optional
        Run a perfect supercell force calculation and subtract the residual
        forces from forces in supercells with displacements. Default is False.
    run_phonopy : Bool, optional
        Whether running phonon calculation or not. Default is False.
    remote_phonopy : Bool, optional
        Whether running phonon calculation or not at remote. Default is False.
    code_string : Str, optional
        Code string of phonopy needed when both of run_phonopy and
        remote_phonopy are True.
    symmetry_tolerance : Float, optional
        Symmetry tolerance. Default is 1e-5.
    immigrant_calculation_folders : Dict, optional
        'force' key has to exist and 'nac' is necessary when
        phonon_settings['is_nac'] is True. The value of the 'force' key is
        the list of strings of remote directories. The value of 'nac' is the
        string of remote directory.
    calculation_nodes : Dict, optional
        This works similarly as immigrant_calculation_folders but contains
        PK or UUID instead of string of remote folder.

    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.input('structure', valid_type=StructureData)
        spec.input('phonon_settings', valid_type=Dict)
        spec.input('symmetry_tolerance', valid_type=Float,
                   default=lambda: Float(1e-5))
        spec.input('dry_run', valid_type=Bool, default=lambda: Bool(False))
        spec.input('subtract_residual_forces', valid_type=Bool,
                   default=lambda: Bool(False))
        spec.input('run_phonopy', valid_type=Bool, default=lambda: Bool(False))
        spec.input('remote_phonopy', valid_type=Bool,
                   default=lambda: Bool(False))
        spec.input('displacement_dataset', valid_type=Dict, required=False)
        spec.input('immigrant_calculation_folders', valid_type=Dict,
                   required=False)
        spec.input('calculation_nodes', valid_type=Dict, required=False)
        spec.input('calculator_settings', valid_type=Dict, required=False)
        spec.input('code_string', valid_type=Str, required=False)

        spec.outline(
            cls.initialize,
            if_(cls.import_calculations)(
                if_(cls.import_calculations_from_files)(
                    while_(cls.continue_import)(
                        cls.read_force_calculations_from_files,
                    ),
                    cls.read_nac_calculations_from_files,
                ),
                if_(cls.import_calculations_from_nodes)(
                    cls.read_calculation_data_from_nodes,
                ),
                cls.check_imported_structures,
            ).else_(
                cls.run_force_and_nac_calculations,
            ),
            if_(cls.dry_run)(
                cls.postprocess_of_dry_run,
            ).else_(
                cls.create_force_sets,
                if_(cls.is_nac)(cls.create_nac_params),
                if_(cls.run_phonopy)(
                    if_(cls.remote_phonopy)(
                        cls.run_phonopy_remote,
                        cls.collect_data,
                    ).else_(
                        cls.create_force_constants,
                        cls.run_phonopy_in_workchain,
                    ),
                ),
            ),
            cls.finalize,
        )
        spec.output('force_constants', valid_type=ArrayData, required=False)
        spec.output('primitive', valid_type=StructureData, required=False)
        spec.output('supercell', valid_type=StructureData, required=False)
        spec.output('force_sets', valid_type=ArrayData, required=False)
        spec.output('supercell_forces', valid_type=ArrayData, required=False)
        spec.output('supercell_energy', valid_type=Float, required=False)
        spec.output('nac_params', valid_type=ArrayData, required=False)
        spec.output('thermal_properties', valid_type=XyData, required=False)
        spec.output('band_structure', valid_type=BandsData, required=False)
        spec.output('dos', valid_type=XyData, required=False)
        spec.output('pdos', valid_type=XyData, required=False)
        spec.output('phonon_setting_info', valid_type=Dict)

    def dry_run(self):
        return self.inputs.dry_run

    def remote_phonopy(self):
        return self.inputs.remote_phonopy

    def run_phonopy(self):
        return self.inputs.run_phonopy

    def is_nac(self):
        if 'is_nac' in self.inputs.phonon_settings.dict:
            return self.inputs.phonon_settings['is_nac']
        else:
            False

    def import_calculations_from_files(self):
        return 'immigrant_calculation_folders' in self.inputs

    def import_calculations_from_nodes(self):
        return 'calculation_nodes' in self.inputs

    def import_calculations(self):
        if 'immigrant_calculation_folders' in self.inputs:
            self.ctx.num_imported = 0
            self.ctx.num_supercell_forces = len(
                self.inputs.immigrant_calculation_folders['forces'])
            return True
        if 'calculation_nodes' in self.inputs:
            return True
        return False

    def continue_import(self):
        return self.ctx.num_imported < self.ctx.num_supercell_forces

    def initialize(self):
        """Set default settings and create supercells and primitive cell"""

        self.report('initialize')

        if self.inputs.run_phonopy and self.inputs.remote_phonopy:
            if 'code_string' not in self.inputs:
                raise RuntimeError("code_string has to be specified.")

        if 'supercell_matrix' not in self.inputs.phonon_settings.dict:
            raise RuntimeError(
                "supercell_matrix was not found in phonon_settings.")

        kwargs = {}
        if 'displacement_dataset' in self.inputs:
            kwargs['dataset'] = self.inputs.displacement_dataset
        return_vals = generate_phonopy_cells(self.inputs.phonon_settings,
                                             self.inputs.structure,
                                             self.inputs.symmetry_tolerance,
                                             **kwargs)

        for key in ('phonon_setting_info', 'primitive', 'supercell'):
            self.ctx[key] = return_vals[key]
            self.out(key, self.ctx[key])
        self.ctx.supercells = {}
        for key in return_vals:
            if "supercell_" in key:
                self.ctx.supercells[key] = return_vals[key]
        if self.inputs.subtract_residual_forces:
            digits = len(str(len(self.ctx.supercells)))
            label = "supercell_%s" % "0".zfill(digits)
            self.ctx.supercells[label] = return_vals['supercell']

    def run_force_and_nac_calculations(self):
        self._run_force_calculations()
        if self.is_nac():
            builder = NacParamsWorkChain.get_builder()
            builder.structure = self.ctx.primitive
            builder.calculator_settings = Dict(dict=self.inputs.calculator_settings['nac'])
            future = self.submit(builder)
            self.report('born_and_epsilon: {}'.format(future.pk))
            self.to_context(**{'born_and_epsilon_calc': future})

    def _run_force_calculations(self):
        """Force calculation"""
        self.report('run force calculations')
        builder_inputs = get_force_calcjob_inputs(
            self.inputs.calculator_settings, self.ctx.supercell)
        for key in self.ctx.supercells:
            builder = get_calcjob_builder(
                self.ctx.supercells[key],
                self.inputs.calculator_settings['forces']['code_string'],
                builder_inputs,
                label=key)
            future = self.submit(builder)
            label = "force_calc_%s" % key.split('_')[-1]
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

    def read_force_calculations_from_files(self):
        self.report('import supercell force calculation data in files.')
        num_batch = 50
        self.report('%d calculations per batch.' % num_batch)

        initial_count = self.ctx.num_imported
        for i in range(initial_count, initial_count + num_batch):
            calc_folders_Dict = self.inputs.immigrant_calculation_folders
            digits = len(str(self.ctx.num_supercell_forces))

            force_folder = calc_folders_Dict['forces'][i]
            label = "force_calc_%s" % str(i + 1).zfill(digits)
            builder = get_immigrant_builder(force_folder,
                                            self.inputs.calculator_settings,
                                            calc_type='forces')
            builder.metadata.label = label
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

            self.ctx.num_imported += 1
            if not self.continue_import():
                break

    def read_nac_calculations_from_files(self):
        if self.is_nac():  # NAC the last one
            self.report('import NAC calculation data in files')
            calc_folders_Dict = self.inputs.immigrant_calculation_folders
            label = 'born_and_epsilon_calc'
            builder = get_immigrant_builder(calc_folders_Dict['nac'][0],
                                            self.inputs.calculator_settings,
                                            calc_type='nac')
            builder.metadata.label = label
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

    def read_calculation_data_from_nodes(self):
        self.report('import calculation data from nodes')

        calc_nodes_Dict = self.inputs.calculation_nodes

        digits = len(str(len(calc_nodes_Dict['force'])))
        for i, node_id in enumerate(calc_nodes_Dict['force']):
            label = "force_calc_%s" % str(i + 1).zfill(digits)
            aiida_node_id = from_node_id_to_aiida_node_id(node_id)
            # self.ctx[label]['forces'] -> ArrayData()('final')
            self.ctx[label] = get_data_from_node_id(aiida_node_id)

        if self.is_nac():
            label = 'born_and_epsilon_cal'
            node_id = calc_nodes_Dict['nac'][0]
            aiida_node_id = from_node_id_to_aiida_node_id(node_id)
            # self.ctx[label]['born_charges'] -> ArrayData()('born_charges')
            # self.ctx[label]['dielectrics'] -> ArrayData()('epsilon')
            self.ctx[label] = get_data_from_node_id(aiida_node_id)

    def check_imported_structures(self):
        self.report('check imported supercell structures')

        msg = ("Immigrant failed because of inconsistency of supercell"
               "structure")

        for key in self.ctx.supercells:
            num = key.split('_')[-1]
            calc = self.ctx["force_calc_%s" % num]
            if type(calc) is dict:
                calc_dict = calc
            else:
                calc_dict = calc.inputs
            supercell_ref = self.ctx.supercells["supercell_%s" % num]
            supercell_calc = calc_dict['structure']
            if not compare_structures(supercell_ref,
                                      supercell_calc,
                                      self.inputs.symmetry_tolerance):
                raise RuntimeError(msg)

    def postprocess_of_dry_run(self):
        self.report('Finish here because of dry-run setting')

    def create_force_sets(self):
        """Build datasets from forces of supercells with displacments"""

        self.report('create force sets')

        forces_dict = collect_vasp_forces_and_energies(
            self.ctx, self.ctx.supercells)
        for key, val in get_vasp_force_sets_dict(**forces_dict).items():
            self.ctx[key] = val
            self.out(key, self.ctx[key])

    def create_nac_params(self):
        self.report('create nac data')

        calc = self.ctx.born_and_epsilon_calc
        if type(calc) is dict:
            calc_dict = calc
            structure = calc['structure']
        else:
            calc_dict = calc.outputs
            structure = calc.inputs.structure

        if 'born_charges' not in calc_dict:
            raise RuntimeError(
                "Born effective charges could not be found "
                "in the calculation. Please check the calculation setting.")
        if 'dielectrics' not in calc_dict:
            raise RuntimeError(
                "Dielectric constant could not be found "
                "in the calculation. Please check the calculation setting.")

        kwargs = {}
        if self.import_calculations():
            kwargs['primitive'] = self.ctx.primitive
        self.ctx.nac_params = get_nac_params(
            calc_dict['born_charges'],
            calc_dict['dielectrics'],
            structure,
            self.inputs.symmetry_tolerance,
            **kwargs)
        self.out('nac_params', self.ctx.nac_params)

    def run_phonopy_remote(self):
        """Run phonopy at remote computer"""

        self.report('remote phonopy calculation')

        code_string = self.inputs.code_string.value
        builder = Code.get_from_string(code_string).get_builder()
        builder.structure = self.inputs.structure
        builder.settings = self.ctx.phonon_setting_info
        builder.metadata.label = self.inputs.metadata.label
        builder.metadata.options.update(self.inputs.phonon_settings['options'])
        builder.force_sets = self.ctx.force_sets
        if 'nac_params' in self.ctx:
            builder.nac_params = self.ctx.nac_params
            builder.primitive = self.ctx.primitive
        future = self.submit(builder)

        self.report('phonopy calculation: {}'.format(future.pk))
        self.to_context(**{'phonon_properties': future})
        # return ToContext(phonon_properties=future)

    def collect_data(self):
        self.report('collect data')
        ph_props = ('thermal_properties',
                    'dos',
                    'pdos',
                    'band_structure',
                    'force_constants')

        for prop in ph_props:
            if prop in self.ctx.phonon_properties.outputs:
                self.out(prop, self.ctx.phonon_properties.outputs[prop])

        self.report('finish phonon')

    def create_force_constants(self):
        self.report('create force constants')

        self.ctx.force_constants = get_force_constants(
            self.inputs.structure,
            self.ctx.phonon_setting_info,
            self.ctx.force_sets)
        self.out('force_constants', self.ctx.force_constants)

    def run_phonopy_in_workchain(self):
        self.report('phonopy calculation in workchain')

        params = {}
        if 'nac_params' in self.ctx:
            params['nac_params'] = self.ctx.nac_params
        result = get_phonon(self.inputs.structure,
                            self.ctx.phonon_setting_info,
                            self.ctx.force_constants,
                            **params)
        self.out('thermal_properties', result['thermal_properties'])
        self.out('dos', result['dos'])
        self.out('band_structure', result['band_structure'])

        self.report('finish phonon')

    def finalize(self):
        self.report('phonopy calculation has been done.')

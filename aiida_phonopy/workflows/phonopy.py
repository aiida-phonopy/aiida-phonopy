"""WorkChan to run phonon calculation by phonopy and force calculators."""

from aiida.engine import WorkChain, if_, while_
from aiida.plugins import DataFactory
from aiida.orm import Code
from aiida_phonopy.common.utils import (
    get_force_constants, get_phonon_properties,
    setup_phonopy_calculation,
    from_node_id_to_aiida_node_id, get_data_from_node_id,
    get_force_sets, collect_forces_and_energies)
from aiida_phonopy.workflows.nac_params import NacParamsWorkChain
from aiida_phonopy.workflows.forces import ForcesWorkChain


# Should be improved by some kind of WorkChainFactory
# For now all workchains should be copied to aiida/workflows

Float = DataFactory('float')
Bool = DataFactory('bool')
Str = DataFactory('str')
Dict = DataFactory('dict')
ArrayData = DataFactory('array')
XyData = DataFactory('array.xy')
StructureData = DataFactory('structure')
BandsData = DataFactory('array.bands')


class PhonopyWorkChain(WorkChain):
    """Phonopy workchain.

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
        fc_calculator : str
            With this being 'alm', ALM is used to calculate force constants in
            the remote phonopy calculation.
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
        """Define inputs, outputs, and outline."""
        super().define(spec)
        spec.input('structure', valid_type=StructureData, required=True)
        spec.input('phonon_settings', valid_type=Dict, required=True)
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
                    cls.initialize_immigrant,
                    while_(cls.continue_import)(
                        cls.read_force_calculations_from_files,
                    ),
                    if_(cls.is_nac)(
                        cls.read_nac_calculations_from_files,
                    ),
                ),
                if_(cls.import_calculations_from_nodes)(
                    cls.read_calculation_data_from_nodes,
                ),
            ).else_(
                cls.run_force_and_nac_calculations,
            ),
            if_(cls.dry_run)(
                cls.postprocess_of_dry_run,
            ).else_(
                cls.create_force_sets,
                if_(cls.is_nac)(cls.attach_nac_params),
                if_(cls.run_phonopy)(
                    if_(cls.remote_phonopy)(
                        cls.run_phonopy_remote,
                        cls.collect_data,
                    ).else_(
                        cls.create_force_constants,
                        cls.run_phonopy_locally,
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
        spec.exit_code(
            1001, 'ERROR_NO_PHONOPY_CODE',
            message=("Phonopy Code not found though expected to run phonopy "
                     "remotely."))
        spec.exit_code(
            1002, 'ERROR_NO_SUPERCELL_MATRIX',
            message=("supercell_matrix was not found."))
        spec.exit_code(
            1003, 'ERROR_INCONSISTENT_IMMIGRANT_FORCES_FOLDERS',
            message=("Number of supercell folders is different from number "
                     "of expected supercells."))

    def dry_run(self):
        """Return boolen for outline."""
        return self.inputs.dry_run

    def remote_phonopy(self):
        """Return boolen for outline."""
        return self.inputs.remote_phonopy

    def run_phonopy(self):
        """Return boolen for outline."""
        return self.inputs.run_phonopy

    def is_nac(self):
        """Return boolen for outline."""
        if 'is_nac' in self.inputs.phonon_settings.keys():
            return self.inputs.phonon_settings['is_nac']
        else:
            False

    def import_calculations_from_files(self):
        """Return boolen for outline."""
        return 'immigrant_calculation_folders' in self.inputs

    def import_calculations_from_nodes(self):
        """Return boolen for outline."""
        return 'calculation_nodes' in self.inputs

    def import_calculations(self):
        """Return boolen for outline."""
        if 'immigrant_calculation_folders' in self.inputs:
            return True
        if 'calculation_nodes' in self.inputs:
            return True
        return False

    def continue_import(self):
        """Return boolen for outline."""
        return self.ctx.num_imported < self.ctx.num_supercell_forces

    def initialize(self):
        """Set default settings and create supercells and primitive cell.

        self.ctx.supercells contains supercells as a dict.
        The keys are like 'spercell_000', 'supercell_001', ...,
        where the number of digits depends on the number of supercells.
        'spercell_000' is only available when
            self.inputs.subtract_residual_forces = True.

        """
        self.report('initialization')

        if self.inputs.run_phonopy and self.inputs.remote_phonopy:
            if 'code_string' not in self.inputs:
                return self.exit_codes.ERROR_NO_PHONOPY_CODE

        if 'supercell_matrix' not in self.inputs.phonon_settings.keys():
            return self.exit_codes.ERROR_NO_SUPERCELL_MATRIX

        kwargs = {}
        if 'displacement_dataset' in self.inputs:
            kwargs['dataset'] = self.inputs.displacement_dataset
        return_vals = setup_phonopy_calculation(self.inputs.phonon_settings,
                                                self.inputs.structure,
                                                self.inputs.symmetry_tolerance,
                                                self.inputs.run_phonopy,
                                                **kwargs)

        for key in ('phonon_setting_info', 'primitive', 'supercell'):
            self.ctx[key] = return_vals[key]
            self.out(key, self.ctx[key])
        self.ctx.supercells = {}
        if self.inputs.subtract_residual_forces:
            digits = len(str(len(self.ctx.supercells)))
            key = "supercell_%s" % "0".zfill(digits)
            self.ctx.supercells[key] = return_vals['supercell']
        for key in return_vals:
            if "supercell_" in key:
                self.ctx.supercells[key] = return_vals[key]

    def initialize_immigrant(self):
        """Initialize immigrant numbers."""
        self.ctx.num_imported = 0
        self.ctx.num_supercell_forces = len(
            self.inputs.immigrant_calculation_folders['forces'])
        self.ctx.supercell_keys_done = []

        if len(self.ctx.supercells) != self.ctx.num_supercell_forces:
            return self.exit_codes.ERROR_INCONSISTENT_IMMIGRANT_FOLDERS

    def run_force_and_nac_calculations(self):
        """Run supercell force and NAC params calculations.

        Supercell force calculations and NAC params calculation
        are submitted in this method to make them run in parallel.

        """
        self._run_force_calculations()
        if self.is_nac():
            self._run_nac_params_calculation()

    def _run_force_calculations(self):
        """Run supercell force calculations."""
        self.report('run force calculations')

        for key, supercell in self.ctx.supercells.items():
            label = "force_calc_%s" % key.split('_')[-1]
            builder = ForcesWorkChain.get_builder()
            builder.metadata.label = label
            builder.structure = supercell
            calculator_settings = self.inputs.calculator_settings['forces']
            builder.calculator_settings = Dict(dict=calculator_settings)
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

    def _run_nac_params_calculation(self):
        """Run nac params calculation."""
        self.report('run nac params calculation')

        builder = NacParamsWorkChain.get_builder()
        builder.metadata.label = 'nac_params'
        builder.structure = self.ctx.primitive
        calculator_settings = self.inputs.calculator_settings['nac']
        builder.calculator_settings = Dict(dict=calculator_settings)
        future = self.submit(builder)
        self.report('nac_params: {}'.format(future.pk))
        self.to_context(**{'nac_params_calc': future})

    def read_force_calculations_from_files(self):
        """Import supercell forces using immigrant."""
        self.report('import supercell force calculation data in files.')
        num_batch = 50
        self.report('%d calculations per batch.' % num_batch)

        calc_folders_Dict = self.inputs.immigrant_calculation_folders
        local_count = 0
        for key, supercell in self.ctx.supercells.items():
            if key in self.ctx.supercell_keys_done:
                continue
            else:
                self.ctx.supercell_keys_done.append(key)
            label = "force_calc_%s" % key.split('_')[-1]
            number = int(key.split('_')[-1])
            if not self.inputs.subtract_residual_forces:
                number -= 1
            force_folder = calc_folders_Dict['forces'][number]
            builder = ForcesWorkChain.get_builder()
            builder.metadata.label = label
            builder.structure = supercell
            calculator_settings = self.inputs.calculator_settings['forces']
            builder.calculator_settings = Dict(dict=calculator_settings)
            builder.immigrant_calculation_folder = Str(force_folder)
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})
            self.ctx.num_imported += 1
            if not self.continue_import():
                break
            local_count += 1
            if local_count == num_batch:
                break

    def read_nac_calculations_from_files(self):
        """Import NAC params using immigrant."""
        self.report('import NAC calculation data in files')
        label = 'nac_params_calc'
        calc_folders_Dict = self.inputs.immigrant_calculation_folders
        nac_folder = calc_folders_Dict['nac'][0]
        builder = NacParamsWorkChain.get_builder()
        builder.metadata.label = label
        builder.structure = self.ctx.primitive
        calculator_settings = self.inputs.calculator_settings['nac']
        builder.calculator_settings = Dict(dict=calculator_settings)
        builder.immigrant_calculation_folder = Str(nac_folder)
        future = self.submit(builder)
        self.report('{} pk = {}'.format(label, future.pk))
        self.to_context(**{label: future})

    def read_calculation_data_from_nodes(self):
        """Read forces and NAC params from data nodes."""
        self.report('import calculation data from nodes')

        calc_nodes_Dict = self.inputs.calculation_nodes

        digits = len(str(len(calc_nodes_Dict['forces'])))
        for i, node_id in enumerate(calc_nodes_Dict['forces']):
            label = "force_calc_%s" % str(i + 1).zfill(digits)
            aiida_node_id = from_node_id_to_aiida_node_id(node_id)
            # self.ctx[label]['forces'] -> ArrayData()('final')
            self.ctx[label] = get_data_from_node_id(aiida_node_id)

        if self.is_nac():
            label = 'nac_params_calc'
            node_id = calc_nodes_Dict['nac'][0]
            aiida_node_id = from_node_id_to_aiida_node_id(node_id)
            self.ctx[label] = get_data_from_node_id(aiida_node_id)

    def postprocess_of_dry_run(self):
        """Show message."""
        self.report('Finish here because of dry-run setting')

    def create_force_sets(self):
        """Build datasets from forces of supercells with displacments."""
        self.report('create force sets')

        forces_dict = collect_forces_and_energies(
            self.ctx, self.ctx.supercells)
        for key, val in get_force_sets(**forces_dict).items():
            self.ctx[key] = val
            self.out(key, self.ctx[key])

    def attach_nac_params(self):
        """Attach nac_params ArrayData to outputs."""
        self.report('create nac params')

        self.ctx.nac_params = self.ctx.nac_params_calc.outputs.nac_params
        self.out('nac_params', self.ctx.nac_params)

    def run_phonopy_remote(self):
        """Run phonopy at remote computer."""
        self.report('remote phonopy calculation')

        code_string = self.inputs.code_string.value
        builder = Code.get_from_string(code_string).get_builder()
        builder.structure = self.inputs.structure
        builder.settings = self.ctx.phonon_setting_info
        builder.symmetry_tolerance = self.inputs.symmetry_tolerance
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
        """Collect phonon data from remove phonopy calculation."""
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
        """Create force constants for run_phonopy_locally."""
        self.report('create force constants')

        self.ctx.force_constants = get_force_constants(
            self.inputs.structure,
            self.ctx.phonon_setting_info,
            self.ctx.force_sets,
            self.inputs.symmetry_tolerance)
        self.out('force_constants', self.ctx.force_constants)

    def run_phonopy_locally(self):
        """Run phonopy calculation locally."""
        self.report('phonopy calculation in workchain')

        params = {}
        if 'nac_params' in self.ctx:
            params['nac_params'] = self.ctx.nac_params
        result = get_phonon_properties(self.inputs.structure,
                                       self.ctx.phonon_setting_info,
                                       self.ctx.force_constants,
                                       self.inputs.symmetry_tolerance,
                                       **params)
        self.out('thermal_properties', result['thermal_properties'])
        self.out('dos', result['dos'])
        self.out('band_structure', result['band_structure'])

        self.report('finish phonon')

    def finalize(self):
        """Show final message."""
        self.report('phonopy calculation has been done.')

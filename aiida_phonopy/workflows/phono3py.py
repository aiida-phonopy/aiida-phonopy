from aiida.engine import WorkChain
from aiida.plugins import WorkflowFactory, DataFactory
from aiida.orm import Float, Bool
from aiida.engine import if_
from aiida_phonopy.common.builders import (
    get_force_calcjob_inputs, get_phonon_force_calcjob_inputs,
    get_nac_calcjob_inputs, get_calcjob_builder, get_immigrant_builder)
from aiida_phonopy.common.utils import (
    generate_phono3py_cells, get_nac_params,
    get_vasp_force_sets_dict, collect_vasp_forces_and_energies,
    compare_structures)


PhonopyWorkChain = WorkflowFactory('phonopy.phonopy')
Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')


class Phono3pyWorkChain(WorkChain):
    """Phono3py workchain


    """
    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(PhonopyWorkChain,
                           exclude=['calculation_nodes',
                                    'run_phonopy',
                                    'remote_phonopy'])
        spec.input('run_phono3py', valid_type=Bool, required=False,
                   default=lambda: Bool(False))
        spec.input('remote_phono3py', valid_type=Bool, required=False,
                   default=lambda: Bool(False))

        spec.outline(
            cls.initialize,
            if_(cls.import_calculations_from_files)(
                cls.read_force_and_nac_calculations_from_files,
                cls.check_imported_structures,
            ).else_(
                cls.run_force_and_nac_calculations,
            ),
            if_(cls.dry_run)(
                cls.postprocess_of_dry_run,
            ).else_(
                cls.create_force_sets,
                if_(cls.is_nac)(cls.create_nac_params),
                if_(cls.run_phono3py)(
                    if_(cls.remote_phono3py)(
                        cls.run_phono3py_remote,
                        cls.collect_data,
                    ).else_(
                        cls.create_force_constants,
                        cls.run_phono3py_in_workchain,
                    )
                )
            )
        )
        spec.output('fc3', valid_type=ArrayData, required=False)
        spec.output('fc2', valid_type=ArrayData, required=False)
        spec.output('primitive', valid_type=StructureData, required=False)
        spec.output('supercell', valid_type=StructureData, required=False)
        spec.output('phonon_supercell', valid_type=StructureData,
                    required=False)
        spec.output('force_sets', valid_type=ArrayData, required=False)
        spec.output('supercell_forces', valid_type=ArrayData, required=False)
        spec.output('supercell_energy', valid_type=Float, required=False)
        spec.output('phonon_force_sets', valid_type=ArrayData, required=False)
        spec.output('phonon_supercell_forces', valid_type=ArrayData,
                    required=False)
        spec.output('phonon_supercell_energy', valid_type=Float,
                    required=False)
        spec.output('nac_params', valid_type=ArrayData, required=False)
        spec.output('phonon_setting_info', valid_type=Dict, required=True)
        spec.exit_code(700, 'ERROR_NO_FORCES_CALCULATOR_SETTING',
                       message='Force calculator setting not found.')
        spec.exit_code(701, 'ERROR_NO_NAC_CALCULATOR_SETTING',
                       message='NAC calculator setting not found.')
        spec.exit_code(
            710, 'ERROR_NO_FORCES_FOLDERS',
            message='Forces calculation folders for immigrant not found.')
        spec.exit_code(
            711, 'ERROR_NO_NAC_FOLDER',
            message='NAC calculation folder for immigrant not found.')

    def is_nac(self):
        if 'is_nac' in self.inputs.phonon_settings.attributes:
            return self.inputs.phonon_settings['is_nac']
        else:
            False

    def dry_run(self):
        return self.inputs.dry_run

    def remote_phono3py(self):
        return self.inputs.remote_phono3py

    def run_phono3py(self):
        return self.inputs.run_phono3py

    def import_calculations_from_files(self):
        return 'immigrant_calculation_folders' in self.inputs

    def initialize(self):
        """Set default settings and create supercells and primitive cell"""

        self.report('initialize')

        if 'forces' not in self.inputs.calculator_settings.attributes:
            self.report("'forces' item is necessary in "
                        "calculator_settings input.")
            return self.exit_code.ERROR_NO_FORCES_CALCULATOR_SETTING

        if self.is_nac():
            if 'nac' not in self.inputs.calculator_settings.attributes:
                self.report("'nac' item is necessary in "
                            "calculator_settings input.")
                return self.exit_code.ERROR_NO_NAC_CALCULATOR_SETTING

        if self.inputs.run_phono3py and self.inputs.remote_phono3py:
            if ('code_string' not in self.inputs or
                'options' not in self.inputs):
                raise RuntimeError(
                    "code_string and options have to be specified.")

        if 'supercell_matrix' not in self.inputs.phonon_settings.attributes:
            raise RuntimeError(
                "supercell_matrix was not found in phonon_settings.")

        kwargs = {}
        if 'displacement_dataset' in self.inputs:
            kwargs['dataset'] = self.inputs.displacement_dataset
        return_vals = generate_phono3py_cells(self.inputs.phonon_settings,
                                              self.inputs.structure,
                                              self.inputs.symmetry_tolerance,
                                              **kwargs)

        for key in ('phonon_setting_info', 'primitive', 'supercell'):
            self.ctx[key] = return_vals[key]
            self.out(key, self.ctx[key])
        self.ctx.supercells = {}
        self.ctx.phonon_supercells = {}
        for key in return_vals:
            if "supercell_" in key and "phonon_" not in key:
                self.ctx.supercells[key] = return_vals[key]
            if "phonon_supercell_" in key:
                self.ctx.phonon_supercells[key] = return_vals[key]
        self.ctx.primitive = return_vals['primitive']
        self.ctx.supercell = return_vals['supercell']
        if 'phonon_supercell' in return_vals:
            self.ctx.phonon_supercell = return_vals['phonon_supercell']

        if self.inputs.subtract_residual_forces:
            digits = len(str(len(self.ctx.supercells)))
            label = "supercell_%s" % "0".zfill(digits)
            self.ctx.supercells[label] = return_vals['supercell']
            digits = len(str(len(self.ctx.phonon_supercells)))
            label = "phonon_supercell_%s" % "0".zfill(digits)
            self.ctx.phonon_supercells[label] = return_vals['phonon_supercell']

    def postprocess_of_dry_run(self):
        self.report('Finish here because of dry-run setting')

    def run_force_and_nac_calculations(self):
        self._run_force_calculations()
        if 'phonon_supercell' in self.ctx:
            self._run_phonon_force_calculations()
        if self.is_nac():
            self._run_nac_calculation()

    def _run_force_calculations(self):
        """FC3 force calculation"""
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

    def _run_phonon_force_calculations(self):
        """FC2 force calculation"""
        self.report('run phonon force calculations')
        calc_settings = self.inputs.calculator_settings
        builder_inputs = get_phonon_force_calcjob_inputs(
            calc_settings, self.ctx.phonon_supercell)
        for key in self.ctx.phonon_supercells:
            builder = get_calcjob_builder(
                self.ctx.phonon_supercells[key],
                calc_settings['phonon_forces']['code_string'],
                builder_inputs,
                label=key)
            future = self.submit(builder)
            label = "phonon_force_calc_%s" % key.split('_')[-1]
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

    def _run_nac_calculation(self):
        """Born charges and dielectric constant calculation"""
        self.report('calculate born charges and dielectric constant')
        builder_inputs = get_nac_calcjob_inputs(
            self.inputs.calculator_settings, self.ctx.primitive)
        builder = get_calcjob_builder(
            self.ctx.primitive,
            self.inputs.calculator_settings['nac']['code_string'],
            builder_inputs,
            label='born_and_epsilon')
        future = self.submit(builder)
        self.report('born_and_epsilon: {}'.format(future.pk))
        self.to_context(**{'born_and_epsilon_calc': future})

    def read_force_and_nac_calculations_from_files(self):
        self.report('import calculation data in files')

        calc_folders_Dict = self.inputs.immigrant_calculation_folders

        if 'forces' not in calc_folders_Dict.attributes:
            return self.exit_code.ERROR_NO_FORCES_FOLDERS

        if self.is_nac():  # NAC the last one
            if 'nac' not in calc_folders_Dict.attributes:
                return self.exit_code.ERROR_NO_NAC_FOLDER

        digits = len(str(len(calc_folders_Dict['forces'])))
        for i, force_folder in enumerate(calc_folders_Dict['forces']):
            label = "force_calc_%s" % str(i + 1).zfill(digits)
            builder = get_immigrant_builder(force_folder,
                                            self.inputs.calculator_settings,
                                            calc_type='forces')
            builder.metadata.label = label
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

        if 'phonon_forces' in calc_folders_Dict.attributes:
            folders = calc_folders_Dict['phonon_forces']
            digits = len(str(len(folders)))
            for i, force_folder in enumerate(folders):
                label = "force_calc_%s" % str(i + 1).zfill(digits)
                builder = get_immigrant_builder(
                    force_folder,
                    self.inputs.calculator_settings,
                    calc_type='phonon_forces')
                builder.metadata.label = label
                future = self.submit(builder)
                self.report('{} pk = {}'.format(label, future.pk))
                self.to_context(**{label: future})

        if self.is_nac():  # NAC the last one
            label = 'born_and_epsilon_calc'
            builder = get_immigrant_builder(calc_folders_Dict['nac'][0],
                                            self.inputs.calculator_settings,
                                            calc_type='nac')
            builder.metadata.label = label
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

    def check_imported_structures(self):
        self.report('check imported supercell structures')

        msg = ("Immigrant failed because of inconsistency of supercell"
               "structure")

        for key in self.ctx.supercells:
            num = key.split('_')[-1]
            calc = self.ctx["force_calc_%s" % num]
            if type(calc) is dict:  # dict when using get_data_from_node_id
                calc_dict = calc
            else:
                calc_dict = calc.inputs
            supercell_ref = self.ctx.supercells["supercell_%s" % num]
            supercell_calc = calc_dict['structure']
            if not compare_structures(
                    supercell_ref,
                    supercell_calc,
                    self.inputs.symmetry_tolerance):
                raise RuntimeError(msg)

        for key in self.ctx.phonon_supercells:
            num = key.split('_')[-1]
            calc = self.ctx["phonon_force_calc_%s" % num]
            if type(calc) is dict:    # dict when using get_data_from_node_id
                calc_dict = calc
            else:
                calc_dict = calc.inputs
            supercell_ref = self.ctx.phonon_supercells[
                "phonon_supercell_%s" % num]
            supercell_calc = calc_dict['structure']
            if not compare_structures(supercell_ref,
                                      supercell_calc,
                                      self.inputs.symmetry_tolerance):
                raise RuntimeError(msg)

    def create_force_sets(self):
        """Build datasets from forces of supercells with displacments"""

        self.report('create force sets')

        forces_dict = collect_vasp_forces_and_energies(
            self.ctx, self.ctx.supercells)
        for key, val in get_vasp_force_sets_dict(**forces_dict).items():
            self.ctx[key] = val
            self.out(key, self.ctx[key])

        forces_dict = collect_vasp_forces_and_energies(
            self.ctx, self.ctx.phonon_supercells, prefix="phonon_force_calc")
        for key, val in get_vasp_force_sets_dict(**forces_dict).items():
            self.ctx["phonon_%s" % key] = val
            self.out("phonon_%s" % key, self.ctx["phonon_%s" % key])

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
        if self.import_calculations_from_files():
            kwargs['primitive'] = self.ctx.primitive
        self.ctx.nac_params = get_nac_params(
            calc_dict['born_charges'],
            calc_dict['dielectrics'],
            structure,
            self.inputs.symmetry_tolerance,
            **kwargs)
        self.out('nac_params', self.ctx.nac_params)

    def run_phono3py_remote(self):
        self.report('remote phonopy calculation')

    def collect_data(self):
        self.report('collect data')

    def create_force_constants(self):
        self.report('create force constants')

    def run_phono3py_in_workchain(self):
        self.report('phonopy calculation in workchain')

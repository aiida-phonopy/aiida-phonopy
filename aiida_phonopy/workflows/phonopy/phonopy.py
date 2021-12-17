"""PhonopyWorkChain."""

from aiida.engine import if_, while_
from aiida.plugins import DataFactory
from aiida_phonopy.workflows.phonopy.base import BasePhonopyWorkChain
from aiida_phonopy.workflows.nac_params import NacParamsWorkChain
from aiida_phonopy.workflows.forces import ForcesWorkChain

Dict = DataFactory("dict")
Str = DataFactory("str")


class PhonopyWorkChain(BasePhonopyWorkChain):
    """Phonopy workchain.

    inputs
    ------
    See most of inputs at BasePhonopyWorkChain.
    calculator_settings : Dict
        Deprecated.
        Settings to run force and nac calculations. For example,
            {'forces': force_config,
             'nac': nac_config}
        At least 'forces' key is necessary. 'nac' is optional.
        force_config is used for supercell force calculation. nac_config
        are used for Born effective charges and dielectric constant calculation
        in primitive cell. The primitive cell is chosen by phonopy
        automatically.
    immigrant_calculation_folders : Dict, optional
        'force' key has to exist and 'nac' is necessary when
        phonon_settings['is_nac'] is True. The value of the 'force' key is
        the list of strings of remote directories. The value of 'nac' is the
        string of remote directory.

    """

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)
        spec.input("calculator_settings", valid_type=Dict, required=False)
        spec.input("immigrant_calculation_folders", valid_type=Dict, required=False)

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
            ).else_(
                cls.run_force_and_nac_calculations,
            ),
            if_(cls.dry_run)(cls.postprocess_of_dry_run,).else_(
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

    def import_calculations_from_files(self):
        """Return boolen for outline."""
        return "immigrant_calculation_folders" in self.inputs

    def import_calculations(self):
        """Return boolen for outline."""
        if "immigrant_calculation_folders" in self.inputs:
            return True
        else:
            return False

    def continue_import(self):
        """Return boolen for outline."""
        return self.ctx.num_imported < self.ctx.num_supercell_forces

    def initialize_immigrant(self):
        """Initialize immigrant numbers."""
        self.ctx.num_imported = 0
        self.ctx.num_supercell_forces = len(
            self.inputs.immigrant_calculation_folders["forces"]
        )
        self.ctx.supercell_keys_done = []

        if len(self.ctx.supercells) != self.ctx.num_supercell_forces:
            return self.exit_codes.ERROR_INCONSISTENT_IMMIGRANT_FOLDERS

    def read_force_calculations_from_files(self):
        """Import supercell forces using immigrant."""
        self.report("import supercell force calculation data in files.")
        num_batch = 50
        self.report("%d calculations per batch." % num_batch)

        calc_folders_Dict = self.inputs.immigrant_calculation_folders
        local_count = 0
        for key, supercell in self.ctx.supercells.items():
            if key in self.ctx.supercell_keys_done:
                continue
            else:
                self.ctx.supercell_keys_done.append(key)
            label = "force_calc_%s" % key.split("_")[-1]
            number = int(key.split("_")[-1])
            if not self.inputs.subtract_residual_forces:
                number -= 1
            force_folder = calc_folders_Dict["forces"][number]
            builder = ForcesWorkChain.get_builder()
            builder.metadata.label = label
            builder.structure = supercell
            calculator_settings = self.inputs.calculator_settings, Str("forces")
            builder.calculator_inputs = calculator_settings
            builder.immigrant_calculation_folder = Str(force_folder)
            future = self.submit(builder)
            self.report("{} pk = {}".format(label, future.pk))
            self.to_context(**{label: future})
            self.ctx.num_imported += 1
            if not self.continue_import():
                break
            local_count += 1
            if local_count == num_batch:
                break

    def read_nac_calculations_from_files(self):
        """Import NAC params using immigrant."""
        self.report("import NAC calculation data in files")
        label = "nac_params_calc"
        calc_folders_Dict = self.inputs.immigrant_calculation_folders
        nac_folder = calc_folders_Dict["nac"][0]
        builder = NacParamsWorkChain.get_builder()
        builder.metadata.label = label
        builder.structure = self.ctx.primitive
        calculator_settings = self.inputs.calculator_settings, Str("nac")
        builder.calculator_inputs = calculator_settings
        builder.immigrant_calculation_folder = Str(nac_folder)
        future = self.submit(builder)
        self.report("{} pk = {}".format(label, future.pk))
        self.to_context(**{label: future})

    def _run_force_calculations(self):
        """Run supercell force calculations."""
        self.report("run force calculations")

        for key, supercell in self.ctx.supercells.items():
            label = "force_calc_%s" % key.split("_")[-1]
            builder = ForcesWorkChain.get_builder()
            builder.metadata.label = label
            builder.structure = supercell
            if "force" in self.inputs.calculator_inputs:
                calculator_inputs = self.inputs.calculator_inputs.force
            else:
                calculator_inputs = self.inputs.calculator_settings["forces"]
                self.logger.warning(
                    "Use calculator_inputs.force instead of "
                    "calculator_settings['forces']."
                )
            builder.calculator_inputs = calculator_inputs
            future = self.submit(builder)
            self.report("{} pk = {}".format(label, future.pk))
            self.to_context(**{label: future})

    def _run_nac_params_calculation(self):
        """Run nac params calculation."""
        self.report("run nac params calculation")

        builder = NacParamsWorkChain.get_builder()
        builder.metadata.label = "nac_params"
        builder.structure = self.ctx.primitive
        if "nac" in self.inputs.calculator_inputs:
            calculator_inputs = self.inputs.calculator_inputs.nac
        else:
            calculator_inputs = self.inputs.calculator_settings["nac"]
            self.logger.warning(
                "Use calculator_inputs.nac instead of calculator_settings['nac']."
            )

        builder.calculator_inputs = calculator_inputs
        future = self.submit(builder)
        self.report("nac_params: {}".format(future.pk))
        self.to_context(**{"nac_params_calc": future})

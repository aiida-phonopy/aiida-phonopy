import numpy as np
from phonopy import Phonopy
from aiida.engine import WorkChain
from aiida.plugins import WorkflowFactory, DataFactory
from aiida.orm import Float, Bool, Str, Code
from aiida.engine import if_
from aiida_phonopy.common.generate_inputs import get_calcjob_builder
from aiida_phonopy.common.utils import (
    generate_phono3py_cells,
    get_vasp_force_sets_dict, collect_vasp_forces_and_energies)


PhonopyWorkChain = WorkflowFactory('phonopy.phonopy')
Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')


class Phono3pyWorkChain(WorkChain):
    """Phono3py workchain


    """
    @classmethod
    def define(cls, spec):
        super(Phono3pyWorkChain, cls).define(spec)
        spec.expose_inputs(PhonopyWorkChain,
                           exclude=['immigrant_calculation_folders',
                                    'calculation_nodes',
                                    'run_phonopy',
                                    'remote_phonopy'])
        spec.input('run_phono3py', valid_type=Bool, required=False,
                   default=lambda: Bool(False))
        spec.input('remote_phono3py', valid_type=Bool, required=False,
                   default=lambda: Bool(False))

        spec.outline(
            cls.initialize,
            cls.run_force_and_nac_calculations,
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
        spec.output('forces_fc3', valid_type=ArrayData, required=False)
        spec.output('forces_fc2', valid_type=ArrayData, required=False)
        spec.output('nac_params', valid_type=ArrayData, required=False)
        spec.output('phonon_setting_info', valid_type=Dict, required=True)

    def dry_run(self):
        return self.inputs.dry_run

    def remote_phono3py(self):
        return self.inputs.remote_phonopy

    def run_phono3py(self):
        return self.inputs.run_phonopy

    def initialize(self):
        """Set default settings and create supercells and primitive cell"""

        self.report('initialize')

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
            if "supercell_" in key:
                self.ctx.supercells[key] = return_vals[key]
            if "phonon_supercell_" in key:
                self.ctx.phonon_supercells[key] = return_vals[key]
        self.ctx.primitive = return_vals['primitive']
        self.ctx.supercell = return_vals['supercell']
        if 'phonon_supercell' in return_vals:
            self.ctx.phonon_supercell = return_vals['phonon_supercell']

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

    def run_phono3py_remote(self):
        self.report('remote phonopy calculation')

    def collect_data(self):
        self.report('collect data')

    def create_force_constants(self):
        self.report('create force constants')

    def run_phono3py_in_workchain(self):
        self.report('phonopy calculation in workchain')

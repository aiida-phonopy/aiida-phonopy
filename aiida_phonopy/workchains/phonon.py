import numpy as np
from aiida.engine import WorkChain, ToContext
from aiida.plugins import DataFactory
from aiida.orm import Float, Bool, Str, Code
from aiida.engine import if_
from aiida_phonopy.common.generate_inputs import (get_calcjob_builder,
                                                  get_immigrant_builder)
from aiida_phonopy.common.utils import (get_force_sets,
                                        get_force_constants,
                                        get_nac_params,
                                        get_phonon,
                                        get_phonon_setting_info,
                                        check_imported_supercell_structure)


# Should be improved by some kind of WorkChainFactory
# For now all workchains should be copied to aiida/workflows

Dict = DataFactory('dict')
ArrayData = DataFactory('array')
XyData = DataFactory('array.xy')
StructureData = DataFactory('structure')
BandsData = DataFactory('array.bands')


class PhononPhonopy(WorkChain):
    """ Workchain to do a phonon calculation using phonopy

    :param structure: StructureData object that contains the crystal structure
        unit cell
    :param calculator_settings: Dict object that contains a
        dictionary with the setting needed to calculate the electronic
        structure:
            {'forces': force_config,
             'nac': nac_config}
        where force_config and nac_config are used for the supercell force
        calculation and Born effective charges and dielectric constant
        calculation in primitive cell, respectively.
    :phonon_settings: Dict object. Needed to run phonon calculation.
        {'mesh': [20, 20, 20]}. Optional. (default: {mesh: [20, 20, 20]})
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
        spec.input('phonon_settings', valid_type=Dict, required=True)
        spec.input('immigrant_calculation_folders',
                   valid_type=Dict, required=False)
        spec.input('calculator_settings',
                   valid_type=Dict, required=False)
        spec.input('code_string', valid_type=Str, required=False)
        spec.input('options', valid_type=Dict, required=False)
        spec.input('symmetry_tolerance',
                   valid_type=Float, required=False, default=Float(1e-5))
        spec.input('dry_run',
                   valid_type=Bool, required=False, default=Bool(False))
        spec.input('run_phonopy',
                   valid_type=Bool, required=False, default=Bool(False))
        spec.input('remote_phonopy',
                   valid_type=Bool, required=False, default=Bool(False))

        spec.outline(cls.initialize_supercell_phonon_calculation,
                     if_(cls.dry_run)(
                         if_(cls.import_calculations)(
                             cls.read_force_and_nac_calculations,
                             cls.check_imported_supercell_structures,
                         )
                     ).else_(
                         if_(cls.import_calculations)(
                             cls.read_force_and_nac_calculations,
                             cls.check_imported_supercell_structures,
                         ).else_(
                             cls.run_force_and_nac_calculations,
                         ),
                         cls.create_force_sets,
                         if_(cls.is_nac)(cls.create_nac_params),
                         if_(cls.run_phonopy)(
                             if_(cls.remote_phonopy)(
                                 cls.run_phonopy_remote,
                                 cls.collect_data,
                             ).else_(
                                 cls.create_force_constants,
                                 cls.run_phonopy_in_workchain,
                             ))))

        spec.output('force_constants', valid_type=ArrayData, required=False)
        spec.output('primitive', valid_type=StructureData, required=False)
        spec.output('supercell', valid_type=StructureData, required=False)
        spec.output('force_sets', valid_type=ArrayData, required=False)
        spec.output('nac_params', valid_type=ArrayData, required=False)
        spec.output('thermal_properties', valid_type=XyData, required=False)
        spec.output('band_structure', valid_type=BandsData, required=False)
        spec.output('dos', valid_type=XyData, required=False)
        spec.output('pdos', valid_type=XyData, required=False)
        spec.output('phonon_setting_info', valid_type=Dict, required=True)

    def dry_run(self):
        return self.inputs.dry_run

    def remote_phonopy(self):
        return self.inputs.remote_phonopy

    def run_phonopy(self):
        return self.inputs.run_phonopy

    def is_nac(self):
        if 'is_nac' in self.inputs.phonon_settings.attributes:
            return self.inputs.phonon_settings['is_nac']
        else:
            False

    def import_calculations(self):
        return 'immigrant_calculation_folders' in self.inputs

    def initialize_supercell_phonon_calculation(self):
        """Set default settings and create supercells and primitive cell"""

        self.report('initialize_supercell_phonon_calculation')

        if self.inputs.run_phonopy and self.inputs.remote_phonopy:
            if ('code_string' not in self.inputs or
                'options' not in self.inputs):
                raise RuntimeError(
                    "code_string and options have to be specified.")

        if 'supercell_matrix' not in self.inputs.phonon_settings.attributes:
            raise RuntimeError(
                "supercell_matrix was not found in phonon_settings.")

        return_vals = get_phonon_setting_info(
            self.inputs.phonon_settings,
            self.inputs.structure,
            self.inputs.symmetry_tolerance)
        self.ctx.phonon_setting_info = return_vals['phonon_setting_info']
        self.out('phonon_setting_info', self.ctx.phonon_setting_info)

        self.ctx.supercells = {}
        for i in range(len(return_vals) - 3):
            label = "supercell_%03d" % (i + 1)
            self.ctx.supercells[label] = return_vals[label]
        self.ctx.primitive = return_vals['primitive']
        self.ctx.supercell = return_vals['supercell']
        self.out('primitive', self.ctx.primitive)
        self.out('supercell', self.ctx.supercell)

    def run_force_and_nac_calculations(self):
        self.report('run force calculations')

        # Forces
        for i in range(len(self.ctx.supercells)):
            label = "supercell_%03d" % (i + 1)
            builder = get_calcjob_builder(self.ctx.supercells[label],
                                          self.inputs.calculator_settings,
                                          calc_type='forces',
                                          label=label)
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

        # Born charges and dielectric constant
        if self.ctx.phonon_setting_info['is_nac']:
            self.report('calculate born charges and dielectric constant')
            builder = get_calcjob_builder(self.ctx.primitive,
                                          self.inputs.calculator_settings,
                                          calc_type='nac',
                                          label='born_and_epsilon')
            future = self.submit(builder)
            self.report('born_and_epsilon: {}'.format(future.pk))
            self.to_context(**{'born_and_epsilon': future})

    def read_force_and_nac_calculations(self):
        self.report('import calculation data in files')

        calc_folders_Dict = self.inputs.immigrant_calculation_folders
        calc_folders = calc_folders_Dict['calculation_folders']
        for i, calc_folder in enumerate(calc_folders[:-1]):
            label = "supercell_%03d" % (i + 1)
            builder = get_immigrant_builder(calc_folder,
                                            self.inputs.calculator_settings,
                                            calc_type='forces')
            builder.metadata.label = label
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

        if self.ctx.phonon_setting_info['is_nac']:  # NAC the last one
            calc_folder = calc_folders[-1]
            label = 'born_and_epsilon'
            builder = get_immigrant_builder(calc_folder,
                                            self.inputs.calculator_settings,
                                            calc_type='nac')
            builder.metadata.label = label
            future = self.submit(builder)
            self.report('{} pk = {}'.format(label, future.pk))
            self.to_context(**{label: future})

    def check_imported_supercell_structures(self):
        self.report('check imported supercell structures')

        msg = ("Immigrant failed because of inconsistency of supercell"
               "structure")

        for i in range(len(self.ctx.supercells)):
            label = "supercell_%03d" % (i + 1)
            supercell_ref = self.ctx.supercells[label]
            supercell_calc = self.ctx[label].inputs.structure
            if not check_imported_supercell_structure(
                    supercell_ref,
                    supercell_calc,
                    self.inputs.symmetry_tolerance):
                raise RuntimeError(msg)

    def create_force_sets(self):
        """Build datasets from forces of supercells with displacments"""

        self.report('create force sets')

        # VASP specific
        forces_dict = {}

        for i in range(len(self.ctx.supercells)):
            label = "supercell_%03d" % (i + 1)
            calc = self.ctx[label]
            if type(calc) is dict:
                out_dict = calc
            else:
                out_dict = calc.outputs
            if ('forces' in out_dict and
                'final' in out_dict['forces'].get_arraynames()):
                label = "forces_%03d" % (i + 1)
                forces_dict[label] = out_dict['forces']
            else:
                msg = ("Forces could not be found in calculation %03d."
                       % (i + 1))
                self.report(msg)

        if len(forces_dict) != len(self.ctx.supercells):
            raise RuntimeError("Forces could not be retrieved.")

        self.ctx.force_sets = get_force_sets(**forces_dict)
        self.out('force_sets', self.ctx.force_sets)

    def create_nac_params(self):
        self.report('create nac data')

        # VASP specific
        # Call workfunction to make links
        out_dict = self.ctx.born_and_epsilon.outputs

        if 'born_charges' not in out_dict:
            raise RuntimeError(
                "Born effective charges could not be found "
                "in the calculation. Please check the calculation setting.")
        if 'dielectrics' not in out_dict:
            raise RuntimeError(
                "Dielectric constant could not be found "
                "in the calculation. Please check the calculation setting.")

        params = {'symmetry_tolerance':
                  Float(self.ctx.phonon_setting_info['symmetry_tolerance'])}
        if self.import_calculations():
            params['primitive'] = self.ctx.primitive
        self.ctx.nac_params = get_nac_params(
            out_dict['born_charges'],
            out_dict['dielectrics'],
            self.ctx.born_and_epsilon.inputs.structure,
            **params)
        self.out('nac_params', self.ctx.nac_params)

    def run_phonopy_remote(self):
        """Run phonopy at remote computer"""

        self.report('remote phonopy calculation')

        code_string = self.inputs.code_string.value
        builder = Code.get_from_string(code_string).get_builder()
        builder.structure = self.inputs.structure
        builder.settings = self.ctx.phonon_setting_info
        builder.metadata.options.update(self.inputs.options)
        builder.metadata.label = self.inputs.metadata.label
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
        self.out('thermal_properties',
                 self.ctx.phonon_properties.outputs.thermal_properties)
        self.out('dos', self.ctx.phonon_properties.outputs.dos)
        self.out('pdos', self.ctx.phonon_properties.outputs.pdos)
        self.out('band_structure',
                 self.ctx.phonon_properties.outputs.band_structure)

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

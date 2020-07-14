from aiida.plugins import DataFactory
from aiida_phonopy.calcs.base import BasePhonopyCalculation
from aiida.orm import Str
from aiida.common import InputValidationError
from aiida_phonopy.common.file_generators import (
    get_FORCE_SETS_txt, get_phonopy_options, get_phonopy_yaml_txt)


BandsData = DataFactory('array.bands')
ArrayData = DataFactory('array')
XyData = DataFactory('array.xy')
Dict = DataFactory('dict')


class PhonopyCalculation(BasePhonopyCalculation):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    _OUTPUT_PROJECTED_DOS = 'projected_dos.dat'
    _OUTPUT_TOTAL_DOS = 'total_dos.dat'
    _OUTPUT_THERMAL_PROPERTIES = 'thermal_properties.yaml'
    _OUTPUT_BAND_STRUCTURE = 'band.yaml'
    _INOUT_FORCE_CONSTANTS = 'FORCE_CONSTANTS'
    _INPUT_CELL = 'phonopy_cells.yaml'
    _INPUT_FORCE_SETS = 'FORCE_SETS'

    @classmethod
    def define(cls, spec):
        super().define(spec)

        spec.input('projected_dos_filename', valid_type=Str,
                   default=lambda: Str(cls._OUTPUT_PROJECTED_DOS))
        spec.input('total_dos_filename', valid_type=Str,
                   default=lambda: Str(cls._OUTPUT_TOTAL_DOS))
        spec.input('thermal_properties_filename', valid_type=Str,
                   default=lambda: Str(cls._OUTPUT_THERMAL_PROPERTIES))
        spec.input('band_structure_filename', valid_type=Str,
                   default=lambda: Str(cls._OUTPUT_BAND_STRUCTURE))
        spec.input('force_constants_filename', valid_type=Str,
                   default=lambda: Str(cls._INOUT_FORCE_CONSTANTS))
        # parser_name has to be set to invoke parsing.
        spec.inputs['metadata']['options']['parser_name'].default = 'phonopy'
        spec.inputs['metadata']['options']['output_filename'].default = 'phonopy.out'

        spec.output('force_constants', valid_type=ArrayData, required=False,
                    help='Calculated force constants')
        spec.output('dos', valid_type=XyData, required=False,
                    help='Calculated total DOS')
        spec.output('pdos', valid_type=XyData, required=False,
                    help='Calculated projected DOS')
        spec.output('thermal_properties', valid_type=XyData, required=False,
                    help='Calculated thermal properties')
        spec.output('band_structure', valid_type=BandsData, required=False,
                    help='Calculated phonon band structure')

    def prepare_for_submission(self, folder):
        return super().prepare_for_submission(folder)

    def _create_additional_files(self, folder):
        self.logger.info("create_additional_files")

        self._create_phonopy_yaml(folder)
        self._create_FORCE_SETS(folder)
        mesh_opts, fc_opts = get_phonopy_options(
            self.inputs.settings.get_dict())

        self._internal_retrieve_list = [self._INOUT_FORCE_CONSTANTS, ]
        self._additional_cmd_params = [['--writefc', ] + fc_opts, ]
        self._calculation_cmd = [['-c', self._INPUT_CELL], ]

        # First run with --writefc, and with --readfc for remaining runs
        if not self.inputs.fc_only:

            self._calculation_cmd = [
                ['-c', self._INPUT_CELL, '--pdos=auto'] + mesh_opts,
                ['-c', self._INPUT_CELL, '-t', '--dos'] + mesh_opts,
                ['-c', self._INPUT_CELL, '--band=auto', '--band-points=101',
                 '--band-const-interval']
            ]
            N = len(self._calculation_cmd) - 1
            self._additional_cmd_params += [['--readfc', ] for i in range(N)]
            self._internal_retrieve_list += [self._OUTPUT_TOTAL_DOS,
                                             self._OUTPUT_PROJECTED_DOS,
                                             self._OUTPUT_THERMAL_PROPERTIES,
                                             self._OUTPUT_BAND_STRUCTURE]

    def _create_phonopy_yaml(self, folder):
        phpy_yaml_txt = get_phonopy_yaml_txt(
            self.inputs.structure,
            self.inputs.settings.get_dict())
        with folder.open(self._INPUT_CELL, 'w', encoding='utf8') as handle:
            handle.write(phpy_yaml_txt)

    def _create_FORCE_SETS(self, folder):
        if 'force_sets' in self.inputs:
            force_sets = self.inputs.force_sets
        else:
            force_sets = None
        if 'displacement_dataset' in self.inputs.settings.attributes:
            dataset = self.inputs.settings['displacement_dataset']
        elif 'dataset' in self.inputs.settings.attributes:
            dataset = self.inputs.settings['dataset']
        elif ('dataset' in self.inputs and
              'displacements' in self.inputs.dataset.get_arraynames()):
            dataset = {'displacements':
                       self.inputs.dataset.get_array('displacements')}
            if 'forces' in self.inputs.dataset.get_arraynames():
                dataset['forces'] = self.inputs.dataset.get_array('forces')
                if force_sets is not None:
                    force_sets = None
        else:
            dataset = None

        # can work both for type-I and type-II
        force_sets_txt = get_FORCE_SETS_txt(dataset, force_sets=force_sets)
        if force_sets_txt is None:
            msg = ("Displacements or forces were not found.")
            raise InputValidationError(msg)

        with folder.open(
                self._INPUT_FORCE_SETS, 'w', encoding='utf8') as handle:
            handle.write(force_sets_txt)

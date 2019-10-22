from aiida.plugins import DataFactory
from aiida_phonopy.calcs.base import BasePhonopyCalculation
from aiida.engine import CalcJob, ExitCode
from aiida.orm import Str
from aiida.common import InputValidationError
from aiida_phonopy.common.raw_parsers import get_FORCE_SETS_txt
import six


BandsData = DataFactory('array.bands')
ArrayData = DataFactory('array')
XyData = DataFactory('array.xy')
Dict = DataFactory('dict')


class PhonopyCalculation(BasePhonopyCalculation, CalcJob):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    _OUTPUT_PROJECTED_DOS = 'projected_dos.dat'
    _OUTPUT_TOTAL_DOS = 'total_dos.dat'
    _OUTPUT_THERMAL_PROPERTIES = 'thermal_properties.yaml'
    _OUTPUT_BAND_STRUCTURE = 'band.yaml'
    _INOUT_FORCE_CONSTANTS = 'FORCE_CONSTANTS'

    @classmethod
    def define(cls, spec):
        super(PhonopyCalculation, cls).define(spec)
        spec.input('projected_dos_filename',
                   valid_type=Str, default=Str(cls._OUTPUT_PROJECTED_DOS))
        spec.input('total_dos_filename',
                   valid_type=Str, default=Str(cls._OUTPUT_TOTAL_DOS))
        spec.input('thermal_properties_filename',
                   valid_type=Str, default=Str(cls._OUTPUT_THERMAL_PROPERTIES))
        spec.input('band_structure_filename',
                   valid_type=Str, default=Str(cls._OUTPUT_BAND_STRUCTURE))
        spec.input('force_constants_filename',
                   valid_type=Str, default=Str(cls._INOUT_FORCE_CONSTANTS))
        spec.input('metadata.options.parser_name',
                   valid_type=six.string_types, default='phonopy')
        spec.input('metadata.options.input_filename',
                   valid_type=six.string_types, default='phonopy.conf')
        spec.input('metadata.options.output_filename',
                   valid_type=six.string_types, default='phonopy.stdout')

        super(PhonopyCalculation, cls)._baseclass_use_methods(spec)

        spec.output('force_constants', valid_type=ArrayData,
                    required=False,
                    help='Calculated force constants')
        spec.output('dos', valid_type=XyData,
                    required=False,
                    help='Calculated total DOS')
        spec.output('pdos', valid_type=XyData,
                    required=False,
                    help='Calculated projected DOS')
        spec.output('thermal_properties', valid_type=XyData,
                    required=False,
                    help='Calculated thermal properties')
        spec.output('band_structure', valid_type=BandsData,
                    required=False,
                    help='Calculated phonon band structure')
        spec.exit_code(100, 'ERROR_NO_RETRIEVED_FOLDER',
                       message=('The retrieved folder data node could not '
                                'be accessed.'))
        spec.exit_code(200, 'ERROR_MISSING_FILE',
                       message='An important file is missing.')
        spec.exit_code(300, 'ERROR_PARSING_FILE_FAILED',
                       message='Parsing a file has failed.')

    def _create_additional_files(self, folder):
        self.logger.info("create_additional_files")

        self._internal_retrieve_list = [self._OUTPUT_TOTAL_DOS,
                                        self._OUTPUT_PROJECTED_DOS,
                                        self._OUTPUT_THERMAL_PROPERTIES,
                                        self._OUTPUT_BAND_STRUCTURE]
        self._calculation_cmd = [
            ['--pdos=auto'],
            ['-t', '--dos'],
            ['--band=auto', '--band-points=101', '--band-const-interval']
        ]

        if 'force_sets' in self.inputs:
            force_sets = self.inputs.force_sets
        else:
            force_sets = None
        if 'displacement_dataset' in self.inputs.settings.attributes:
            displacement_dataset = self.inputs.settings[
                'displacement_dataset']
        else:
            displacement_dataset = None

        if force_sets is not None and displacement_dataset is not None:
            force_sets_txt = get_FORCE_SETS_txt(
                force_sets, displacement_dataset)
            force_sets_filename = folder.get_abs_path(self._INPUT_FORCE_SETS)
            with open(force_sets_filename, 'w') as infile:
                infile.write(force_sets_txt)
            # First run with --writefc, and with --readfc for remaining runs
            self._additional_cmd_params = [
                ['--readfc'] for i in range(len(self._calculation_cmd) - 1)]
            self._additional_cmd_params.insert(0, ['--writefc'])
            self._internal_retrieve_list.append(self._INOUT_FORCE_CONSTANTS)
        else:
            msg = ("no force_sets nor force_constants are specified for "
                   "this calculation")
            raise InputValidationError(msg)

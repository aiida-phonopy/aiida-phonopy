from aiida.plugins import DataFactory
from aiida_phonopy.calculations.phonopy.base import BasePhonopyCalculation
from aiida.engine import CalcJob
from aiida.orm import Str

from aiida.common import InputValidationError
from aiida_phonopy.common.raw_parsers import (get_FORCE_CONSTANTS_txt,
                                              get_FORCE_SETS_txt)


PhononDosData = DataFactory('phonopy.phonon_dos')
BandStructureData = DataFactory('phonopy.band_structure')
ArrayData = DataFactory('array')
Dict = DataFactory('dict')


class PhonopyCalculation(BasePhonopyCalculation, CalcJob):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    _INOUT_FORCE_CONSTANTS = 'FORCE_CONSTANTS'

    _OUTPUT_DOS = 'projected_dos.dat'
    _OUTPUT_THERMAL_PROPERTIES = 'thermal_properties.yaml'
    _OUTPUT_BAND_STRUCTURE = 'band.yaml'
    _default_parser = 'phonopy'

    _internal_retrieve_list = [_OUTPUT_DOS,
                               _OUTPUT_THERMAL_PROPERTIES,
                               _OUTPUT_BAND_STRUCTURE]
    _calculation_cmd = [['--pdos=auto'], ['--readfc', '-t']]

    @classmethod
    def define(cls, spec):
        super(PhonopyCalculation, cls).define(spec)
        spec.input('projected_dos_filename',
                   valid_type=Str, default=Str('projected_dos.dat'))
        spec.input('thermal_properties_filename',
                   valid_type=Str, default=Str('thermal_properties.yaml'))
        spec.input('band_structure_filename',
                   valid_type=Str, default=Str('band.yaml'))
        spec.input('force_constants_filename',
                   valid_type=Str, default=Str(cls._INOUT_FORCE_CONSTANTS))

        super(PhonopyCalculation, cls)._baseclass_use_methods(spec)

        spec.output('force_constants', valid_type=ArrayData,
                    required=False,
                    help='Calculated force constants')
        spec.output('dos', valid_type=PhononDosData,
                    required=False,
                    help='Calculated force constants')
        spec.output('thermal_properties', valid_type=ArrayData,
                    required=False,
                    help='Calculated force constants')
        spec.output('band_structure', valid_type=BandStructureData,
                    required=False,
                    help='Calculated force constants')
        spec.exit_code(100, 'ERROR_NO_RETRIEVED_FOLDER',
                       message=('The retrieved folder data node could not '
                                'be accessed.'))
        spec.exit_code(200, 'ERROR_MISSING_FILE',
                       message='An important file is missing.')
        spec.exit_code(300, 'ERROR_PARSING_FILE_FAILED',
                       message='Parsing a file has failed.')

    def _create_additional_files(self, folder):
        self.logger.info("create_additional_files")

        self._additional_cmdline_params = []

        if 'force_sets' in self.inputs:
            force_sets = self.inputs.force_sets
        else:
            force_sets = None
        if 'displacement_dataset' in self.inputs:
            displacement_dataset = self.inputs.displacement_dataset
        else:
            displacement_dataset = None
        if 'force_constants' in self.inputs:
            force_constants = self.inputs.force_constants
        else:
            force_constants = None

        if force_constants is not None:
            force_constants_txt = get_FORCE_CONSTANTS_txt(force_constants)
            force_constants_filename = folder.get_abs_path(
                self._INOUT_FORCE_CONSTANTS)
            with open(force_constants_filename, 'w') as infile:
                infile.write(force_constants_txt)
            self._additional_cmdline_params += ['--readfc']
        elif force_sets is not None and displacement_dataset is not None:
            force_sets_txt = get_FORCE_SETS_txt(force_sets,
                                                displacement_dataset)
            force_sets_filename = folder.get_abs_path(
                self._INPUT_FORCE_SETS)
            with open(force_sets_filename, 'w') as infile:
                infile.write(force_sets_txt)
            self._additional_cmdline_params += ['--writefc']
            self._internal_retrieve_list += ['FORCE_CONSTANTS']
        else:
            msg = ("no force_sets nor force_constants are specified for "
                   "this calculation")
            raise InputValidationError(msg)

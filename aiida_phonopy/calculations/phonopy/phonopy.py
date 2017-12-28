from aiida.common.utils import classproperty
from aiida.orm import DataFactory
from aiida_phonopy.calculations.phonopy import BasePhonopyCalculation
from aiida.orm.calculation.job import JobCalculation

from aiida.common.exceptions import InputValidationError
from aiida_phonopy.common.raw_parsers import get_FORCE_CONSTANTS_txt, get_FORCE_SETS_txt


KpointsData = DataFactory('array.kpoints')


class PhonopyCalculation(BasePhonopyCalculation, JobCalculation):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    _INOUT_FORCE_CONSTANTS = 'FORCE_CONSTANTS'

    _OUTPUT_DOS = 'partial_dos.dat'
    _OUTPUT_THERMAL_PROPERTIES = 'thermal_properties.yaml'
    _OUTPUT_BAND_STRUCTURE = 'band.yaml'

    def _init_internal_params(self):
        super(PhonopyCalculation, self)._init_internal_params()

        self._default_parser = 'phonopy'

        self._internal_retrieve_list += [self._OUTPUT_DOS,
                                         self._OUTPUT_THERMAL_PROPERTIES,
                                         self._OUTPUT_BAND_STRUCTURE]

        self._calculation_cmd = [['--pdos=0'], ['-t']]

    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        retdict = JobCalculation._use_methods
        retdict.update(BasePhonopyCalculation._baseclass_use_methods)

        #retdict['bands'] = {
        #    'valid_types': BandStructureData,
        #    'additional_parameter': None,
        #    'linkname': 'bands',
        #    'docstring': "Use the node defining the band structure to use",
        #}

        return retdict

    def _create_additional_files(self, tempfolder, inputdict):

        data_sets = inputdict.pop(self.get_linkname('data_sets'), None)
        force_constants = inputdict.pop(self.get_linkname('force_constants'), None)
        bands = inputdict.pop(self.get_linkname('bands'), None)

        if data_sets is None and force_constants is None:
            raise InputValidationError("no force_sets nor force_constants are specified for this calculation")

        if data_sets is not None:
            force_sets_txt = get_FORCE_SETS_txt(data_sets)
            force_sets_filename = tempfolder.get_abs_path(self._INPUT_FORCE_SETS)
            with open(force_sets_filename, 'w') as infile:
                infile.write(force_sets_txt)
            self._additional_cmdline_params += ['--writefc']
            self._internal_retrieve_list += ['FORCE_CONSTANTS']

        if force_constants is not None:
            force_constants_txt = get_FORCE_CONSTANTS_txt(force_constants)
            force_constants_filename = tempfolder.get_abs_path(self._INOUT_FORCE_CONSTANTS)
            with open(force_constants_filename, 'w') as infile:
                infile.write(force_constants_txt)
            self._additional_cmdline_params += ['--readfc']


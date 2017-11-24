from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.orm import DataFactory
from aiida_phonopy.calculations.phonopy import BasePhonopyCalculation


BandStructureData = DataFactory('phonopy.band_structure')
KpointsData = DataFactory('array.kpoints')

class PhonopyCalculation(BasePhonopyCalculation, JobCalculation):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    def _init_internal_params(self):
        super(PhonopyCalculation, self)._init_internal_params()

        self._default_parser = 'phonopy'

        self._internal_retrieve_list += [self._OUTPUT_DOS,
                                         self._OUTPUT_THERMAL_PROPERTIES,
                                         self._OUTPUT_BAND_STRUCTURE]

        self._calculation_cmd = ['--pdos=0', '-t']

    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        retdict = JobCalculation._use_methods
        retdict.update(BasePhonopyCalculation._baseclass_use_methods)

        retdict['bands'] = {
            'valid_types': BandStructureData,
            'additional_parameter': None,
            'linkname': 'bands',
            'docstring': "Use the node defining the band structure to use",
        }

        return retdict
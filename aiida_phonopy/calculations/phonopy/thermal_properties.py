from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.orm import DataFactory
from aiida_phonopy.calculations.phonopy import BasePhonopyCalculation

KpointsData = DataFactory('array.kpoints')


class ThermalPropertiesCalculation(BasePhonopyCalculation, JobCalculation):
    """
    A basic plugin for calculating force constants using Phonopy.
    Requirement: the node should be able to import phonopy
    """

    _OUTPUT_THERMAL_PROPERTIES = 'thermal_properties.yaml'

    def _init_internal_params(self):
        super(ThermalPropertiesCalculation, self)._init_internal_params()

        self._default_parser = 'phonopy'
        self._additional_cmdline_params += ['-t']

        self._internal_retrieve_list += [self._OUTPUT_THERMAL_PROPERTIES]

    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        retdict = JobCalculation._use_methods
        retdict.update(BasePhonopyCalculation._baseclass_use_methods)

#        retdict['kpoints'] = {
#            'valid_types': KpointsData,
#            'additional_parameter': None,
#            'linkname': 'kpoints',
#            'docstring': "Use the node defining the kpoint sampling to use",
#        }

        return retdict
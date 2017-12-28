from aiida.common.exceptions import InputValidationError
from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.orm import DataFactory
from aiida_phonopy.calculations.phonopy import BasePhonopyCalculation

from aiida_phonopy.common.raw_parsers import get_disp_fc3_txt, get_forces_txt

BandStructureData = DataFactory('phonopy.band_structure')
KpointsData = DataFactory('array.kpoints')

class Phono3pyCalculation(BasePhonopyCalculation, JobCalculation):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    _INPUT_DISP_FC3 = 'disp_fc3.yaml'
    _INPUT_FORCES_FC3 = 'FORCES_FC3'
    _OUTPUT_KAPPA = 'kappa-'

    def _init_internal_params(self):
        super(Phono3pyCalculation, self)._init_internal_params()

        self._default_parser = 'phono3py'

        self._calculation_cmd = [['--br']]

    @classproperty
    def _use_methods(cls):
        """
        Extend the parent _use_methods with further keys.
        """
        retdict = JobCalculation._use_methods
        retdict.update(BasePhonopyCalculation._baseclass_use_methods)

        return retdict

    def _create_additional_files(self, tempfolder, inputdict):

        data_sets = inputdict.pop(self.get_linkname('data_sets'), None)
        force_constants = inputdict.pop(self.get_linkname('force_constants'), None)

        if data_sets is None and force_constants is None:
            raise InputValidationError("no force_sets nor force_constants are specified for this calculation")

        parameters_data = inputdict[self.get_linkname('parameters')]
        structure = inputdict[self.get_linkname('structure')]

        if data_sets is not None:
            disp_fc3_filename = tempfolder.get_abs_path(self._INPUT_DISP_FC3)
            disp_f3c_txt = get_disp_fc3_txt(structure, parameters_data, data_sets)
            with open(disp_fc3_filename, 'w') as infile:
                infile.write(disp_f3c_txt)

            forces_txt = get_forces_txt(data_sets)
            forces_filename = tempfolder.get_abs_path(self._INPUT_FORCES_FC3)
            with open(forces_filename, 'w') as infile:
                infile.write(forces_txt)

        # Set name for output file
        kappa_filename = self._OUTPUT_KAPPA + 'm{}{}{}.hdf5'.format(*parameters_data.dict.mesh)
        self._internal_retrieve_list += [kappa_filename]
        print self._internal_retrieve_list
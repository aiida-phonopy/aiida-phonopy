from aiida.common.exceptions import InputValidationError
from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.orm import DataFactory
from aiida_phonopy.calculations.phonopy import BasePhonopyCalculation


from aiida_phonopy.common.raw_parsers import get_disp_fc3_txt, get_forces_txt, \
    write_fc2_to_hdf5_file, write_fc3_to_hdf5_file

BandStructureData = DataFactory('phonopy.band_structure')
KpointsData = DataFactory('array.kpoints')
ForceConstantsData = DataFactory('phonopy.force_constants')

class Phono3pyCalculation(BasePhonopyCalculation, JobCalculation):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    _INPUT_DISP_FC3 = 'disp_fc3.yaml'
    _INPUT_FORCES_FC3 = 'FORCES_FC3'
    _INPUT_FC3 = 'fc3.hdf5'
    _INPUT_FC2 = 'fc2.hdf5'
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

        retdict['force_constants_3'] = {
            'valid_types': ForceConstantsData,
            'additional_parameter': None,
            'linkname': 'force_constants_2nd',
            'docstring': "Use the node defining the 2nd order force constants",
        }

        return retdict

    def _create_additional_files(self, tempfolder, inputdict):

        data_sets = inputdict.pop(self.get_linkname('data_sets'), None)
        fc3 = inputdict.pop(self.get_linkname('force_constants_3'), None)
        fc2 = inputdict.pop(self.get_linkname('force_constants'), None)

        parameters_data = inputdict[self.get_linkname('parameters')]
        structure = inputdict[self.get_linkname('structure')]

        if fc3 and fc2 is not None:
            fc2_filename = tempfolder.get_abs_path(self._INPUT_FC2)
            write_fc2_to_hdf5_file(fc2, fc2_filename)

            fc3_filename = tempfolder.get_abs_path(self._INPUT_FC3)
            write_fc3_to_hdf5_file(fc3, fc3_filename)

            self._additional_cmdline_params += ['--fc2', '--fc3']

        elif data_sets is not None:
            disp_fc3_filename = tempfolder.get_abs_path(self._INPUT_DISP_FC3)
            disp_f3c_txt = get_disp_fc3_txt(structure, parameters_data, data_sets)
            with open(disp_fc3_filename, 'w') as infile:
                infile.write(disp_f3c_txt)

            forces_txt = get_forces_txt(data_sets)
            forces_filename = tempfolder.get_abs_path(self._INPUT_FORCES_FC3)
            with open(forces_filename, 'w') as infile:
                infile.write(forces_txt)
        else:
            raise InputValidationError("either force_sets or force_constants should be specified for this calculation")

        # Set name for output file
        kappa_filename = self._OUTPUT_KAPPA + 'm{}{}{}.hdf5'.format(*parameters_data.dict.mesh)
        self._internal_retrieve_list += [kappa_filename]
        print self._internal_retrieve_list
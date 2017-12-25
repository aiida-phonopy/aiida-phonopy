from aiida.common.exceptions import InputValidationError
from aiida.orm.calculation.job import JobCalculation
from aiida.common.utils import classproperty
from aiida.orm import DataFactory
from aiida_phonopy.calculations.phonopy import BasePhonopyCalculation


BandStructureData = DataFactory('phonopy.band_structure')
KpointsData = DataFactory('array.kpoints')

class Phono3pyCalculation(BasePhonopyCalculation, JobCalculation):
    """
    A basic plugin for calculating phonon properties using Phonopy.
    """

    _INPUT_DISP_FC3 = 'disp_fc3.yaml'


    def _init_internal_params(self):
        super(Phono3pyCalculation, self)._init_internal_params()

        self._default_parser = 'phono3py'

        self._internal_retrieve_list += [self._OUTPUT_DOS,
                                         self._OUTPUT_THERMAL_PROPERTIES,
                                         self._OUTPUT_BAND_STRUCTURE]

        self._calculation_cmd = ['--fc3', '--fc2', '--br']


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
            get_disp_fc3_txt(structure, parameters_data, data_sets, disp_fc3_filename)

        if force_constants is not None:
            force_constants_txt = get_FORCE_CONSTANTS_txt(force_constants)
            force_constants_filename = tempfolder.get_abs_path(self._INOUT_FORCE_CONSTANTS)
            with open(force_constants_filename, 'w') as infile:
                infile.write(force_constants_txt)
            self._additional_cmdline_params += ['--readfc']


def get_disp_fc3_txt(structure, parameters_data, data_sets, filename):
        import phono3py.file_IO
        from phono3py.file_IO import write_disp_fc3_yaml
        from phonopy.structure.cells import get_supercell
        from aiida_phonopy.workchains.phonon import phonopy_bulk_from_structure

        supercell = get_supercell(phonopy_bulk_from_structure(structure),
                                  parameters_data.dict.supercell,
                                  symprec=parameters_data.dict.symmetry_precision)

        write_disp_fc3_yaml(data_sets.get_data_sets3(),
                            supercell,
                            filename=filename)

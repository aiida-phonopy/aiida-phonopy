from aiida.common import InputValidationError
from aiida.engine import CalcJob
from aiida.common.utils import classproperty
from aiida.plugins import DataFactory
from aiida_phonopy.calcs.base import BasePhonopyCalculation
from aiida_phonopy.common.file_generators import (
    get_disp_fc3_txt, get_forces_txt, write_fc2_to_hdf5_file,
    write_fc3_to_hdf5_file, write_kappa_to_hdf5_file)

BandStructureData = DataFactory('phonopy.band_structure')
KpointsData = DataFactory('array.kpoints')
ArrayData = DataFactory('array')
ForceConstantsData = DataFactory('phonopy.force_constants')


def get_grid_data_files(grid_data):

    gp_data = {}
    for key in grid_data.get_arraynames():
        array = grid_data.get_array(key)
        num = key.split('_')[-1].replace('g', '')
        array_name = '_'.join(key.split('_')[:-1])
        try:
            gp_data[num].update({array_name: array})
        except KeyError:
            gp_data[num] = {array_name: array}

    return gp_data


class Phono3pyCalculation(BasePhonopyCalculation, CalcJob):
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

        retdict['grid_data'] = {
            'valid_types': ArrayData,
            'additional_parameter': None,
            'linkname': 'grid_data',
            'docstring': ("Use the node to include grid points data in "
                          "distributed calculation"),
        }

        return retdict

    def _create_additional_files(self, tempfolder, inputdict):

        datasets = inputdict.pop(self.get_linkname('datasets'), None)
        fc3 = inputdict.pop(self.get_linkname('force_constants_3'), None)
        fc2 = inputdict.pop(self.get_linkname('force_constants'), None)

        parameters_data = inputdict[self.get_linkname('parameters')]
        structure = inputdict[self.get_linkname('structure')]

        self._additional_cmdline_params = ['--thm']
        #self._internal_retrieve_list = []

        if fc3 and fc2 is not None:
            fc2_filename = tempfolder.get_abs_path(self._INPUT_FC2)
            write_fc2_to_hdf5_file(fc2, fc2_filename)

            fc3_filename = tempfolder.get_abs_path(self._INPUT_FC3)
            write_fc3_to_hdf5_file(fc3, fc3_filename)

            self._additional_cmdline_params += ['--fc2', '--fc3']

        elif datasets is not None:
            disp_fc3_filename = tempfolder.get_abs_path(self._INPUT_DISP_FC3)
            disp_f3c_txt = get_disp_fc3_txt(
                structure, parameters_data, datasets)
            with open(disp_fc3_filename, 'w') as infile:
                infile.write(disp_f3c_txt)

            forces_txt = get_forces_txt(datasets)
            forces_filename = tempfolder.get_abs_path(self._INPUT_FORCES_FC3)
            with open(forces_filename, 'w') as infile:
                infile.write(forces_txt)
        else:
            msg = ("either force_sets or force_constants should be specified "
                   "for this calculation")
            raise InputValidationError(msg)

        grid_data = inputdict.pop(self.get_linkname('grid_data'), None)

        if grid_data is not None:
            gp_data = get_grid_data_files(grid_data)
            for gp in gp_data:
                kappa_g_filename = self._OUTPUT_KAPPA + \
                               'm{}{}{}'.format(*parameters_data.dict.mesh) + \
                               '-g{}.hdf5'.format(gp)
                write_kappa_to_hdf5_file(
                    gp_data[gp],
                    filename=tempfolder.get_abs_path(kappa_g_filename))

            self._additional_cmdline_params += ['--read-gamma']
        else:
            if 'grid_point' in parameters_data.get_dict():
                gp_string = ','.join(
                    [str(gp) for gp in parameters_data.dict.grid_point])
                self._additional_cmdline_params += [
                    '--gp={}'.format(gp_string), '--write-gamma']
                for gp in parameters_data.dict.grid_point:
                    kappa_g_filename = (
                        self._OUTPUT_KAPPA +
                        'm{}{}{}'.format(*parameters_data.dict.mesh) +
                        '-g{}.hdf5'.format(gp))
                    self._internal_retrieve_list += [kappa_g_filename]
                return

        # Set name for output file
        kappa_filename = (self._OUTPUT_KAPPA +
                          'm{}{}{}.hdf5'.format(*parameters_data.dict.mesh))
        self._internal_retrieve_list += [kappa_filename]

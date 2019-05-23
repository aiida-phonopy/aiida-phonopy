from aiida.common import InputValidationError
from aiida.common import CalcInfo, CodeInfo
from aiida.common.utils import classproperty
from aiida.plugins import DataFactory
from aiida_phonopy.common.raw_parsers import (get_BORN_txt,
                                              get_phonopy_conf_file_txt,
                                              get_poscar_txt)

ParameterData = DataFactory('dict')
StructureData = DataFactory('structure')
ArrayData = DataFactory('array')
BandStructureData = DataFactory('phonopy.band_structure')


class BasePhonopyCalculation(object):
    """
    A basic plugin for calculating force constants using Phonopy.
    Requirement: the node should be able to import phonopy if NAC is used
    """

    _INPUT_FILE_NAME = 'phonopy.conf'
    _INPUT_CELL = 'POSCAR'
    _INPUT_FORCE_SETS = 'FORCE_SETS'
    _INPUT_NAC = 'BORN'

    # initialize with default files that should always be retrieved,
    # additional files are added in the specific plugin
    _internal_retrieve_list = []

    # Initialize list of commands to be specified for each specific
    # plugin (lists should be empty)
    _additional_cmdline_params = []
    _calculation_cmd = []

    @classmethod
    def _baseclass_use_methods(cls, spec):
        spec.input('parameters', valid_type=ParameterData,
                   help=('Use a node that specifies the phonopy '
                         'parameters for the namelists'))
        spec.input('structure', valid_type=StructureData,
                   help=('Use a node for the structure'))
        spec.input('force_sets', valid_type=ArrayData,
                   required=False, default=None,
                   help=('Use a node that specifies the force_sets '
                         'array for the namelists'))
        spec.input('displacement_dataset', valid_type=ParameterData,
                   required=False, default=None,
                   help=('Use a node that specifies the dispalcement '
                         'dataset parameters for the namelists'))
        spec.input('force_constants', valid_type=ArrayData,
                   required=False, default=None,
                   help=('Use a node that specifies the force constants '
                         'arrays for the namelists'))
        spec.input('nac_params', valid_type=ArrayData,
                   required=False, default=None,
                   help=('Use a node for the Non-analitical '
                         'corrections parameters arraies'))
        spec.input('bands', valid_type=BandStructureData,
                   required=False,
                   default=None,
                   help='Use the node defining the band structure to use')

        # return {
        #     "parameters": {
        #         'valid_types': ParameterData,
        #         'additional_parameter': None,
        #         'linkname': 'parameters',
        #         'docstring': ("Use a node that specifies the phonopy "
        #                       "parameters for the namelists"),
        #     },
        #     "force_sets": {
        #         'valid_types': ArrayData,
        #         'additional_parameter': None,
        #         'linkname': 'force_sets',
        #         'docstring': ("Use a node that specifies the force_sets "
        #                       "array for the namelists"),
        #     },
        #     "displacement_dataset": {
        #         'valid_types': ParameterData,
        #         'additional_parameter': None,
        #         'linkname': 'displacement_dataset',
        #         'docstring': ("Use a node that specifies the dispalcement "
        #                       "dataset parameters for the namelists"),
        #     },
        #     "force_constants": {
        #         'valid_types': ArrayData,
        #         'additional_parameter': None,
        #         'linkname': 'force_constants',
        #         'docstring': ("Use a node that specifies the force constants "
        #                       "arraies for the namelists"),
        #     },
        #     "structure": {
        #         'valid_types': StructureData,
        #         'additional_parameter': None,
        #         'linkname': 'structure',
        #         'docstring': "Use a node for the structure",
        #     },
        #     "nac_params": {
        #         'valid_types': ArrayData,
        #         'additional_parameter': None,
        #         'linkname': 'nac_params',
        #         'docstring': ("Use a node for the Non-analitical "
        #                       "corrections parameters arraies"),
        #     },
        #     'bands': {
        #         'valid_types': BandStructureData,
        #         'additional_parameter': None,
        #         'linkname': 'bands',
        #         'docstring': "Use the node defining the band structure to use",
        #     }
        # }

    def _create_additional_files(self, tempfolder, inputs_params):
        pass

    def prepare_for_submission(self, folder, inputdict):
        """Create the input files from the input nodes passed to this instance of the `CalcJob`.

        :param folder: an `aiida.common.folders.Folder` to temporarily write files on disk
        :return: `aiida.common.datastructures.CalcInfo` instance
        """

        parameters_data = self.inputs.parameters
        structure = self.inputs.structure
        code = self.inputs.code
        nac_params = self.inputs.nac_params
        bands = self.inputs.bands

        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################

        # ================= prepare the python input files =================

        self._create_additional_files(folder, inputdict)

        cell_txt = get_poscar_txt(structure)
        input_txt = get_phonopy_conf_file_txt(parameters_data, bands=bands)

        input_filename = folder.get_abs_path(self._INPUT_FILE_NAME)
        with open(input_filename, 'w') as infile:
            infile.write(input_txt)

        cell_filename = folder.get_abs_path(self._INPUT_CELL)
        with open(cell_filename, 'w') as infile:
            infile.write(cell_txt)

        if nac_params is not None:
            born_txt = get_BORN_txt(nac_params, parameters_data, structure)
            nac_filename = folder.get_abs_path(self._INPUT_NAC)
            with open(nac_filename, 'w') as infile:
                infile.write(born_txt)
            self._additional_cmdline_params += ['--nac']

        # ============================ calcinfo ===============================

        local_copy_list = []
        remote_copy_list = []
        # additional_retrieve_list = settings_dict.pop("ADDITIONAL_RETRIEVE_LIST",[])

        calcinfo = CalcInfo()

        calcinfo.uuid = self.uuid

        # Empty command line by default
        calcinfo.local_copy_list = local_copy_list
        calcinfo.remote_copy_list = remote_copy_list

        # Retrieve files
        calcinfo.retrieve_list = self._internal_retrieve_list

        calcinfo.codes_info = []
        for property_cmd in self._calculation_cmd:
            codeinfo = CodeInfo()
            codeinfo.cmdline_params = ([self._INPUT_FILE_NAME]
                                       + self._additional_cmdline_params
                                       + property_cmd)
            codeinfo.code_uuid = code.uuid
            codeinfo.withmpi = False
            calcinfo.codes_info += [codeinfo]

        return calcinfo

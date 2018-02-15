from aiida.orm.calculation.job import JobCalculation

from aiida.common.exceptions import InputValidationError
from aiida.common.datastructures import CalcInfo, CodeInfo
from aiida.common.utils import classproperty

from aiida.orm import DataFactory

ForceSetsData = DataFactory('phonopy.force_sets')

ForceConstantsData = DataFactory('phonopy.force_constants')
ParameterData = DataFactory('parameter')
StructureData = DataFactory('structure')
ArrayData = DataFactory('array')
BandStructureData = DataFactory('phonopy.band_structure')

from aiida_phonopy.common.raw_parsers import get_BORN_txt, get_FORCE_CONSTANTS_txt, get_FORCE_SETS_txt, \
    get_phonopy_conf_file_txt, get_poscar_txt


class BasePhonopyCalculation(object):
    """
    A basic plugin for calculating force constants using Phonopy.
    Requirement: the node should be able to import phonopy if NAC is used
    """

    _INPUT_FILE_NAME = 'phonopy.conf'
    _INPUT_CELL = 'POSCAR'
    _INPUT_FORCE_SETS = 'FORCE_SETS'
    _INPUT_NAC = 'BORN'


    # initialize with default files that should always be retrieved, additional files are added in the specific plugin
    _internal_retrieve_list = []

    # Initialize list of commands to be specified for each specific plugin (lists should be empty)
    _additional_cmdline_params = []
    _calculation_cmd = []

    @classproperty
    def _baseclass_use_methods(cls):
        """
        Added to to each subclass
        """
        return {
            "parameters": {
                'valid_types': ParameterData,
                'additional_parameter': None,
                'linkname': 'parameters',
                'docstring': ("Use a node that specifies the phonopy parameters "
                              "for the namelists"),
            },
            "data_sets": {
                'valid_types': ForceSetsData,
                'additional_parameter': None,
                'linkname': 'data_sets',
                'docstring': ("Use a node that specifies the data_sets "
                              "for the namelists"),
            },
            "force_constants": {
                'valid_types': ForceConstantsData,
                'additional_parameter': None,
                'linkname': 'force_constants',
                'docstring': ("Use a node that specifies the data_sets "
                              "for the namelists"),
            },
            "structure": {
                'valid_types': StructureData,
                'additional_parameter': None,
                'linkname': 'structure',
                'docstring': "Use a node for the structure",
            },
            "nac_data": {
                'valid_types': ArrayData,
                'additional_parameter': None,
                'linkname': 'nac_data',
                'docstring': "Use a node for the Non-analitical corrections data",
            },
            'bands': {
                'valid_types': BandStructureData,
                'additional_parameter': None,
                'linkname': 'bands',
                'docstring': "Use the node defining the band structure to use",
            }
        }


    def _create_additional_files(self, tempfolder, inputs_params):
        pass

    def _prepare_for_submission(self, tempfolder, inputdict):
        """
        This is the routine to be called when you want to create
        the input files and related stuff with a plugin.
        :param tempfolder: a aiida.common.folders.Folder subclass where
                           the plugin should put all its files.
        :param inputdict: a dictionary with the input nodes, as they would
                be returned by get_inputdata_dict (without the Code!)
        """

        try:
            parameters_data = inputdict[self.get_linkname('parameters')]
        except KeyError:
            raise InputValidationError("No parameters specified for this calculation")

        try:
            structure = inputdict[self.get_linkname('structure')]
        except KeyError:
            raise InputValidationError("no structure is specified for this calculation")

        try:
            code = inputdict.pop(self.get_linkname('code'))
        except KeyError:
            raise InputValidationError("no code is specified for this calculation")

        nac_data = inputdict.pop(self.get_linkname('nac_data'), None)
        bands = inputdict.pop(self.get_linkname('bands'), None)


        ##############################
        # END OF INITIAL INPUT CHECK #
        ##############################

        # =================== prepare the python input files =====================

        self._create_additional_files(tempfolder, inputdict)

        cell_txt = get_poscar_txt(structure)
        input_txt = get_phonopy_conf_file_txt(parameters_data, bands=bands)

        input_filename = tempfolder.get_abs_path(self._INPUT_FILE_NAME)
        with open(input_filename, 'w') as infile:
            infile.write(input_txt)

        cell_filename = tempfolder.get_abs_path(self._INPUT_CELL)
        with open(cell_filename, 'w') as infile:
            infile.write(cell_txt)

        if nac_data is not None:
            born_txt = get_BORN_txt(nac_data, structure=structure, parameters=parameters_data)
            nac_filename = tempfolder.get_abs_path(self._INPUT_NAC)
            with open(nac_filename, 'w') as infile:
                infile.write(born_txt)
            self._additional_cmdline_params += ['--nac']

        # ============================ calcinfo ================================

        local_copy_list = []
        remote_copy_list = []
        #    additional_retrieve_list = settings_dict.pop("ADDITIONAL_RETRIEVE_LIST",[])

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
            codeinfo.cmdline_params = [self._INPUT_FILE_NAME] + self._additional_cmdline_params + property_cmd
            codeinfo.code_uuid = code.uuid
            codeinfo.withmpi = False
            calcinfo.codes_info += [codeinfo]

        return calcinfo
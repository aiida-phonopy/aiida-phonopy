# -*- coding: utf-8 -*-
"""CalcJob for phonopy post-processing."""

from aiida import orm
from aiida.common import InputValidationError, datastructures
from aiida.engine import CalcJob

from aiida_phonopy.data import ForceConstantsData, PhonopyData
from aiida_phonopy.utils.mapping import _lowercase_dict, _uppercase_dict


def get_default_metadata_options():
    """Get a default metadata option for Phonopy calculation."""


class PhonopyCalculation(CalcJob):
    """Base `CalcJob` implementation for Phonopy post-processing."""

    # Mapping keys for parsers, giving respectively output file name and output node name.
    _OUTPUTS = {
        # format: .hdf5
        'fc': ('force_constants.hdf5', 'output_force_constants'),
        'band': ('band.hdf5', 'phonon_bands'),
        'qpoints': ('qpoints.hdf5', 'qpoints'),
        'mesh': ('mesh.hdf5', 'qpoints_mesh'),
        # format: .yaml
        'irreps': ('irreps.yaml', 'irreducible_representations'),
        'tprop': ('thermal_properties.yaml', 'thermal_properties'),
        'tdisp': ('thermal_displacements.yaml', 'thermal_displacements'),
        'tdispmat': ('thermal_displacement_matrices.yaml', 'thermal_displacement_matrices'),
        'mod': ('modulation.yaml', 'modulation'),
        # format: .dat
        'dos': ('total_dos.dat', 'total_phonon_dos'),
        'pdos': ('projected_dos.dat', 'projected_phonon_dos')
    }

    _INPUT_FORCE_CONSTANTS = 'force_constants.hdf5'

    _DEFAULT_INPUT_FILE = 'aiida.in'
    _DEFAULT_OUTPUT_FILE = 'aiida.out'
    _DEFAULT_PHONOPY_FILE = 'phonopy.yaml'

    _INPUT_SUBFOLDER = './'
    _OUTPUT_SUBFOLDER = './'

    # Available tags - further selection has to be made
    _AVAILABLE_TAGS = {
        # Basic tags
        'PRIMITIVE_AXES': [list],  # delete?
        'PRIMITIVE_AXIS': [list],
        'EIGENVECTORS': [bool],
        # Band structure tags
        'BAND': [str, list],
        'BAND_PATHS': [str, list],
        'BAND_POINTS': [float],
        'BAND_LABELS': [list],
        'BAND_CONNECTION': [bool],
        'BAND_INDICES': [list],
        # Mehs sampling tags
        'MESH': [list],
        'MP': [int, float],
        'MESH_NUMBERS': [int, float],  # ?
        'MP_SHIFT': [list],  # ?
        'GAMMA_CENTER': [bool],
        'WRITE_MESH': [bool],
        # Phonon density of states (DOS) tags
        'DOS': [bool],
        'DOS_RANGE': [str, list],
        'FMIN': [int, float],
        'FMAX': [int, float],
        'FPITCH': [int, float],
        'PDOS': [str],
        'PROJECTION_DIRECTION': [list],
        'XYZ_DIRECTION': [bool],
        'SIGMA': [int, float],
        'DEBYE_MODEL': [bool],
        'MOMEMT': [bool],
        'MOMENT_ORDER': [int],
        # Thermal properties and displacements related tags
        'TPROP': [bool],
        'TMIN': [int, float],
        'TMAX': [int, float],
        'TSTEP': [int, float],
        'PRETEND_REAL': [bool],
        'CUTOFF_FREQUENCY': [int, float],
        'TDISP': [bool],
        'TDISPMAT': [bool],
        'TDISPMAT_CIF': [int, float],
        # Specific q-points
        'QPOINTS': [list],
        'WRITEDM': [bool],
        # Non-analytical term correction
        'NAC_METHOD': [str],
        'Q_DIRECTION': [list],
        # Group velocity
        'GROUP_VELOCITY': [bool],
        'GV_DELTA_Q': [int, float],
        # Symmetry
        'SYMMETRY_TOLERANCE': [int, float],
        'SYMMETRY': [bool],
        'MESH_SYMMETRY': [bool],
        'FC_SYMMETRY': [bool],
        # Force constants
        'FULL_FORCE_CONSTANTS': [bool],
        'WRITE_FORCE_CONSTANTS': [bool],
        # Create animation file
        'ANIME_TYPE': [str],
        'ANIME': [list],
        # Create modulated structure
        'MODULATION': [str],
        # Characters of irreducible representations
        'IRREPS': [list],
        'SHOW_IRREPS': [bool],
        'LITTLE_COGROUP': [bool],
    }

    # Keywords that cannot be set by the user but will be set by the plugin
    _BLOCKED_TAGS = [
        'DIM',
        'ATOM_NAME',
        'MASS',  #??? from structure or also from here??
        'MAGMOM',
        'CREATE_DISPLACEMENTS',
        'DISPLACEMENT_DISTANCE',
        'DIAG',
        'PM',
        'RANDOM_DISPLACEMENTS',
        'RANDOM_SEED',
        'NAC',
        'FORCE_CONSTANTS',
        'READ_FORCE_CONSTANTS',
        'FC_FORMAT',
        'READFC_FORMAT',
        'WRITEFC_FORMAT',
        'BAND_FORMAT',
        'MESH_FORMAT',
        'QPOINTS_FORMAT',
        'HDF5',
    ]

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)

        spec.input(
            'parameters',
            valid_type=orm.Dict,
            required=True,
            help=(
                'Phonopy parameters (`setting tags`) for post processing. '
                'The following tags, along their type, are allowed:\n' +
                '\n'.join(f'{tag_name}' for tag_name in cls._AVAILABLE_TAGS)
            ),
            validator=cls._validate_parameters,
        )
        spec.input(
            'phonopy_data',
            valid_type=PhonopyData,
            required=False,
            help='The preprocess output info of a previous ForceConstantsWorkChain.'
        )
        spec.input(
            'force_constants',
            valid_type=ForceConstantsData,
            required=False,
            help='Force constants of the input structure.'
        )
        spec.input('settings', valid_type=orm.Dict, required=False, help='Settings for phonopy calculation.')
        # spec.inputs.validator = cls._validate_inputs

        spec.input('metadata.options.withmpi', valid_type=bool, default=False)
        spec.input('metadata.options.input_filename', valid_type=str, default=cls._DEFAULT_INPUT_FILE)
        spec.input('metadata.options.output_filename', valid_type=str, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input('metadata.options.parser_name', valid_type=str, default='phonopy.phonopy')
        spec.inputs['metadata']['options']['resources'].default = lambda: {'num_machines': 1}

        spec.output(
            'output_parameters', valid_type=orm.Dict, required=False, help='Sum up info of phonopy calculation.'
        )
        spec.output(cls._OUTPUTS['fc'][1], valid_type=orm.ArrayData, required=False, help='Calculated force constants.')
        spec.output(cls._OUTPUTS['dos'][1], valid_type=orm.XyData, required=False, help='Calculated total DOS.')
        spec.output(cls._OUTPUTS['pdos'][1], valid_type=orm.XyData, required=False, help='Calculated projected DOS.')
        spec.output(
            cls._OUTPUTS['band'][1], valid_type=orm.BandsData, required=False, help='Calculated phonon band structure.'
        )
        spec.output(cls._OUTPUTS['qpoints'][1], valid_type=orm.BandsData, required=False, help='Calculated qpoints.')
        spec.output(cls._OUTPUTS['mesh'][1], valid_type=orm.BandsData, required=False, help='Calculated qpoint mesh.')
        spec.output(
            cls._OUTPUTS['irreps'][1], valid_type=orm.Dict, required=False, help='Irreducible representation output.'
        )
        spec.output(
            cls._OUTPUTS['tprop'][1], valid_type=orm.XyData, required=False, help='Calculated thermal properties.'
        )
        spec.output(
            cls._OUTPUTS['tdisp'][1], valid_type=orm.Dict, required=False, help='Calculated thermal displacements.'
        )
        spec.output(
            cls._OUTPUTS['tdispmat'][1],
            valid_type=orm.Dict,
            required=False,
            help='Calculated thermal displacements matrices.'
        )
        spec.output(cls._OUTPUTS['mod'][1], valid_type=orm.Dict, required=False, help='Modulation information.')

        # Unrecoverable errors: required retrieved files could not be read, parsed or are otherwise incomplete
        spec.exit_code(
            301, 'ERROR_NO_RETRIEVED_TEMPORARY_FOLDER', message='The retrieved temporary folder could not be accessed.'
        )
        spec.exit_code(
            302,
            'ERROR_OUTPUT_STDOUT_MISSING',
            message='The retrieved folder did not contain the required stdout output file.'
        )
        spec.exit_code(
            303,
            'ERROR_OUTPUT_PHONOPY_MISSING',
            message='The retrieved folder did not contain the required phonopy file.'
        )
        spec.exit_code(
            304,
            'ERROR_OUTPUT_FILES_MISSING',
            message='The retrieved folder did not contain one or more expected output files.'
        )
        spec.exit_code(
            305, 'ERROR_BAD_INPUTS', message='No run mode has been selected.'
        )  # maybe we should check this from the input too

        spec.exit_code(310, 'ERROR_OUTPUT_STDOUT_READ', message='The stdout output file could not be read.')
        spec.exit_code(311, 'ERROR_OUTPUT_STDOUT_PARSE', message='The stdout output file could not be parsed.')
        spec.exit_code(
            312,
            'ERROR_OUTPUT_STDOUT_INCOMPLETE',
            message='The stdout output file was incomplete probably because the calculation got interrupted.'
        )

        spec.exit_code(320, 'ERROR_OUTPUT_YAML_LOAD', message='The loading of yaml file got an unexpected error.')
        spec.exit_code(321, 'ERROR_OUTPUT_NUMPY_LOAD', message='The file loading via numpy got an unexpected error.')

        spec.exit_code(350, 'ERROR_UNEXPECTED_PARSER_EXCEPTION', message='The parser raised an unexpected exception.')

        # Not implemented file parser
        spec.exit_code(
            400, 'ERROR_NOT_IMPLEMENTED_PARSER', message='The parser was not able to parse one or more files.'
        )

    @classmethod
    def _validate_parameters(cls, value, _):
        """Validate the ``parameters`` input namespace."""

        def __validate_dict(value_dict):
            enabled_dict = cls._AVAILABLE_TAGS
            unknown_tags = set(value_dict.keys()) - set(enabled_dict.keys())
            if unknown_tags:
                return (
                    f"Unknown tags in 'parameters': {unknown_tags}, "
                    f'allowed tags are {cls._AVAILABLE_TAGS.keys()}.'
                )
            invalid_values = [
                value_dict[key] for key in value_dict.keys() if not type(value_dict[key]) in enabled_dict[key]
            ]
            if invalid_values:
                return f'Parameters tags must be of the correct type; got invalid values {invalid_values}.'

        if value:
            if isinstance(value, orm.Dict):
                __validate_dict(value.get_dict())

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :py:`~aiida.common.datastructures.CalcInfo` instance.
        """
        retrieve_list = []
        retrieve_temporary_list = []
        calculation_cmds = []

        if 'settings' in self.inputs:
            settings = _lowercase_dict(self.inputs.settings.get_dict(), dict_name='settings')
        else:
            settings = {}
        # It can contain:
        # "keep_animation_files": bool
        # "symmetrize_nac": bool
        # "factor_nac": float
        # "subtract_residual_forces": bool

        # ================= prepare the phonopy input files ===================

        self.write_phonopy_info(folder)

        if 'force_constants' in self.inputs:
            self.write_force_constants(folder)

        # =============== prepare the submit and conf file ====================

        parameters = _uppercase_dict(self.inputs.parameters.get_dict(), dict_name='parameters')
        filename = self.inputs.metadata.options.input_filename
        self.write_calculation_input(folder, parameters, filename)

        # from phonopy recommendations, filename of the inputs configuration goes first
        calculation_cmds.append([filename])

        # ================= retreiving phonopy output files ===================

        # Keeping in the repository only the stdout, as deafult.
        # Animation files and `phonopy.yaml` can be retrieved preparing
        # the `settings` input namespace accordingly.
        retrieve_list.append(self.inputs.metadata.options.output_filename)
        if settings.pop('keep_animation_files', None):
            for format_ in ['jmol', 'xyz', 'xyz_jmol', 'arc', 'ascii']:
                retrieve_list.append(f'anime.{format_}')
            retrieve_list.append('APOSCAR-*')

        if settings.pop('keep_phonopy_yaml', False):
            retrieve_list.append(self._DEFAULT_PHONOPY_FILE)
        else:
            retrieve_temporary_list.append(self._DEFAULT_PHONOPY_FILE)

        # Retrieving everything and raising error from parser if something is missing.
        for value in self._OUTPUTS.values():
            retrieve_temporary_list.append(value[0])  # first value of the tuple
        if 'force_constants' in self.inputs:  # otherwise we retrieve the same thing as the input
            retrieve_temporary_list.pop(self._INPUT_FORCE_CONSTANTS)

        # ============================ calcinfo ===============================

        local_copy_list = []

        calcinfo = datastructures.CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = local_copy_list  # what we write in the folder
        calcinfo.retrieve_list = retrieve_list  # what to retrieve and keep
        # what to retrieve temporarly for just parsing purposes
        calcinfo.retrieve_temporary_list = retrieve_temporary_list

        calcinfo.codes_info = []
        for cmdline_params in calculation_cmds:
            codeinfo = datastructures.CodeInfo()  # new code info object per cmdline
            codeinfo.cmdline_params = cmdline_params
            codeinfo.stdout_name = self.inputs.metadata.options.output_filename
            codeinfo.code_uuid = self.inputs.code.uuid
            codeinfo.withmpi = self.inputs.metadata.options.withmpi
            # appending CodeInfo to CalcInfo
            calcinfo.codes_info.append(codeinfo)

        return calcinfo

    def _get_p2s_map(self):
        """Get the primitive to supercell map."""
        if 'force_constants' in self.inputs:
            return self.inputs.force_constants.get_cells_mappings()['primitive']['p2s_map']
        if 'phonopy_data' in self.inputs:
            return self.inputs.phonopy_data.get_cells_mappings()['primitive']['p2s_map']
        return None

    def write_phonopy_info(self, folder):
        """Write in `folder` the `phonopy.yaml` file."""
        from phonopy.interface.phonopy_yaml import PhonopyYaml

        kwargs = {}

        if 'settings' in self.inputs:
            the_settings = self.inputs.settings.get_dict()
            for key in ['symmetrize_nac', 'factor_nac', 'subtract_residual_forces']:
                if key in the_settings:
                    kwargs.update({key: the_settings[key]})

        if 'phonopy_data' in self.inputs:
            ph = self.inputs.phonopy_data.get_phonopy_instance(**kwargs)
        elif 'force_constants' in self.inputs:
            ph = self.inputs.force_constants.get_phonopy_instance(**kwargs)

        # Setting the phonopy yaml obtject to produce yaml lines
        # .. note: this does not write the force constants
        phpy_yaml = PhonopyYaml()
        phpy_yaml.set_phonon_info(ph)
        phpy_yaml_txt = str(phpy_yaml)

        with folder.open(self._DEFAULT_PHONOPY_FILE, 'w', encoding='utf8') as handle:
            handle.write(phpy_yaml_txt)

    def write_force_constants(self, folder):
        """Write in `folder` the force constants file."""
        from phonopy.file_IO import write_force_constants_to_hdf5

        filename = folder.get_abs_path(self._INPUT_FORCE_CONSTANTS)

        write_force_constants_to_hdf5(
            force_constants=self.inputs.force_constants.force_constants, filename=filename, p2s_map=self._get_p2s_map()
        )

    def write_calculation_input(self, folder, parameters: dict, filename: str):
        """Write in `folder` the input file containing the information regarding the calculation."""
        characters = []  # single characters to be written

        # appending eventual blocked keys to parameters
        if 'force_constants' in self.inputs:
            parameters.append({'FORCE_CONSTANTS': 'READ'})

        if 'nac_parameters' in self.inputs:
            parameters.append({'NAC': True})

        # if not "SYMMETRY_TOLERANCE" in parameters.keys():
        #     symprec = self.inputs.phonopy_data.symprec
        #     parameters.update({"SYMMETRY_TOLERANCE":symprec})

        # write potential huge outputs in `.hdf5` format, and read force constants from this format too
        parameters.update({'HDF5': True})

        # append single character
        for key, value in parameters.items():
            characters.append(key)
            characters.append('=')
            if isinstance(value, str):
                characters.append(value)
            elif isinstance(value, bool):
                characters.append(f'.{value}.')
            elif isinstance(value, (float, int)):
                characters.append(str(value))
            elif isinstance(value, list):
                string_value = ''
                for val in value:
                    string_value += str(val) + ' '
                characters.append(string_value[:-1])
            else:
                raise InputValidationError(f'value type `{type(value)}` is not supported.')
            characters.append('\n')  # one line for each key in parameters

        with folder.open(filename, 'w', encoding='utf8') as handle:
            for character in characters:
                handle.write(character)

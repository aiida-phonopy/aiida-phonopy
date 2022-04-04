# -*- coding: utf-8 -*-
"""CalcJob for phonopy post-processing."""

import os

from aiida import orm
from aiida.common import datastructures, InputValidationError

from aiida_phonopy.utils.file_generators import get_FORCE_SETS_txt, get_phonopy_yaml_txt, get_BORN_txt
from aiida_phonopy.calculations.base import BasePhonopyCalculation, _uppercase_dict


class PhonopyCalculation(BasePhonopyCalculation):
    """`CalcJob` implementation of the Phonopy post-processing calculations."""

    # Input/Output keys
    _OUTPUTS = {
        "band": "band.yaml",
        "qpoints": "qpoints.yaml",
        "irreps": "irreps.yaml",
        "dos": "total_dos.dat",  # not yaml
        "pdos": "projected_dos.dat",  # not yaml
        "tprop": "thermal_properties.yaml",
        "tdisp": "thermal_displacements.yaml",
        "tdispmat": "thermal_displacement_matrices.yaml",
        "mod": "modulation.yaml",  # we should not parse it (?)
        "mesh": "mehs.yaml",  # we should not parse it (?)
    }

    _INOUT_FORCE_CONSTANTS = "force_constants.hdf5"
    _INPUT_FORCE_SETS = "FORCE_SETS"
    _DEFAULT_INPUT_FILE = "phonopy.conf"
    _DEFAULT_CELL_FILE = "phonopy_cells.yaml"
    _DEFAULT_OUTPUT_FILE = "phonopy.yaml"

    _INPUT_SUBFOLDER = "./"
    _OUTPUT_SUBFOLDER = "./"

    # Available tags - further selection has to be made
    _AVAILABLE_TAGS = {
        # Basic tags
        "PRIMITIVE_AXES": [list],  # delete?
        "PRIMITIVE_AXIS": [list],
        "EIGENVECTORS": [bool],
        # Band structure tags
        "BAND": [str, list],
        "BAND_PATHS": [str, list],
        "BAND_POINTS": [float],
        "BAND_LABELS": [list],
        "BAND_CONNECTION": [bool],
        "BAND_INDICES": [list],
        # Mehs sampling tags
        "MESH": [list],
        "MP": [int, float],
        "MESH_NUMBERS": [int, float],  # ?
        "MP_SHIFT": [list],  # ?
        "GAMMA_CENTER": [bool],
        "WRITE_MESH": [bool],
        # Phonon density of states (DOS) tags
        "DOS": [bool],
        "DOS_RANGE": [str, list],
        "FMIN": [int, float],
        "FMAX": [int, float],
        "FPITCH": [int, float],
        "PDOS": [str],
        "PROJECTION_DIRECTION": [list],
        "XYZ_DIRECTION": [bool],
        "SIGMA": [int, float],
        "DEBYE_MODEL": [bool],
        "MOMEMT": [bool],
        "MOMENT_ORDER": [int],
        # Thermal properties and displacements related tags
        "TPROP": [str],
        "TMIN": [int, float],
        "TMAX": [int, float],
        "TSTEP": [int, float],
        "PRETEND_REAL": [bool],
        "CUTOFF_FREQUENCY": [int, float],
        "TDISP": [bool],
        "TDISPMAT": [bool],
        "TDISPMAT_CIF": [int, float],
        # Specific q-points
        "QPOINTS": [list],
        "WRITEDM": [bool],
        # Non-analytical term correction
        "NAC_METHOD": [str],
        "Q_DIRECTION": [list],
        # Group velocity
        "GROUP_VELOCITY": [bool],
        "GV_DELTA_Q": [int, float],
        # Symmetry
        "SYMMETRY_TOLERANCE": [int, float],  # !!!! ATTENTION TO INPUTS !!!!
        "SYMMETRY": [bool],
        "MESH_SYMMETRY": [bool],
        "FC_SYMMETRY": [bool],
        # Force constants
        "FULL_FORCE_CONSTANTS": [bool],
        # Create animation file
        "ANIME_TYPE": [str],
        "ANIME": [list],
        # Create modulated structure
        "MODULATION": [str],
        # Characters of irreducible representations
        "IRREPS": [list],
        "SHOW_IRREPS": [bool],
        "LITTLE_COGROUP": [bool],
    }

    # Keywords that cannot be set by the user but will be set by the plugin
    _BLOCKED_TAGS = [
        "DIM",
        "ATOM_NAME",
        "MASS",  # ??? from structure or also from here??
        "MAGMOM",
        "CREATE_DISPLACEMENTS",
        "DISPLACEMENT_DISTANCE",
        "DIAG",
        "PM",
        "RANDOM_DISPLACEMENTS",
        "RANDOM_SEED",
        "NAC",
        "FORCE_CONSTANTS",
        "READ_FORCE_CONSTANTS",
        "WRITE_FORCE_CONSTANTS",
        "FC_FORMAT",
        "READFC_FORMAT",
        "WRITEFC_FORMAT",
        "BAND_FORMAT",
        "MESH_FORMAT",
        "QPOINTS_FORMAT",
        "HDF5",
    ]

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)

        spec.input("fc_calculator", valid_type=orm.Str, required=False)
        spec.input("metadata.options.input_filename", valid_type=str, default=cls._DEFAULT_INPUT_FILE)
        spec.input("metadata.options.output_filename", valid_type=str, default=cls._DEFAULT_OUTPUT_FILE)
        spec.input("metadata.options.parser_name", valid_type=str, default="phonopy.phonopy")

        spec.output("parameters", valid_type=orm.Dict, required=False, help="Sum up info of phonopy calculation.")
        spec.output("force_constants", valid_type=orm.ArrayData, required=False, help="Calculated force constants")
        spec.output("dos", valid_type=orm.XyData, required=False, help="Calculated total DOS")
        spec.output("pdos", valid_type=orm.XyData, required=False, help="Calculated projected DOS")
        spec.output("thermal_properties", valid_type=orm.Dict, required=False, help="Calculated thermal properties")
        spec.output("band", valid_type=orm.Dict, required=False, help="Calculated phonon band structure")
        # To be yet well defined
        spec.output("irreps", valid_type=orm.Dict, required=False, help="Irreducible representation information.")
        spec.output("qpoints", valid_type=orm.Dict, required=False, help="Qpoints information.")
        spec.output("modulation", valid_type=orm.Dict, required=False, help="Modulation information.")

        spec.exit_code(
            301, "ERROR_NO_RETRIEVED_TEMPORARY_FOLDER", message="The retrieved temporary folder could not be accessed."
        )
        spec.exit_code(302, "ERROR_NO_FORCE_CONSTANTS", message="The force constants file is missing.")
        spec.exit_code(303, "ERROR_NO_PHONOPY_YAML", message="The phonopy yaml file is missing.")

    def prepare_for_submission(self, folder):
        """Prepare the calculation job for submission by transforming input nodes into input files.

        In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
        contains lists of files that need to be copied to the remote machine before job submission, as well as file
        lists that are to be retrieved after job completion.

        :param folder: a sandbox folder to temporarily write files on disk.
        :return: :py:`~aiida.common.datastructures.CalcInfo` instance.
        """

        retrieve_temporary_list = []
        # This will depend on what there is in the input `parameters`.
        # For each `dict` in `parameters` the `calculation_cmds` is updated.
        calculation_cmds = []

        # ================= prepare the phonopy input files ===================

        self._write_born(folder)
        self._write_cells_info(folder)
        self._write_force_sets(folder)

        # =============== prepare the submit and conf file ====================

        cmd_readcell = [
            "-c",
            self._DEFAULT_CELL_FILE,
            "--tolerance",
            str(self.inputs.symmetry_tolerance.value),
        ]  # in front any phonopy cmd line
        cmd_writefc = ["--writefc", "--writefc-format=hdf5"]  # first cmd if `parent_folder` is not specified
        cmd_readfc = ["--readfc", "--readfc-format=hdf5"]  # in front any phonopy post-process cmd line

        if "fc_calculator" in self.inputs:
            value = self.inputs["fc_calculator"].value.lower()
            if value == "alm":
                cmd_writefc.append("--alm")
            else:
                raise NotImplementedError(f"force calculator `{value}` is not supported.")

        if "parent_folder" not in self.inputs:
            calculation_cmds.append(cmd_readcell + cmd_writefc)

        parameters = self.inputs["parameters"]
        if isinstance(parameters, orm.Dict):
            _params = _uppercase_dict(parameters.get_dict(), dict_name="parameters")
            if _params:
                filename = self.inputs.metadata.options.input_filename
                self._write_conf_info(folder, _params, filename)
                calculation_cmds.append([filename] + cmd_readcell + cmd_readfc)
        elif isinstance(parameters, orm.List):
            _iparam = 0
            for params in parameters.get_list():
                _params = _uppercase_dict(params, dict_name=f"parameters[{_iparam}]")
                filename = self.inputs.metadata.options.input_filename + f"_{_iparam}"
                if _params:
                    self._write_conf_info(folder, _params, filename)
                    calculation_cmds.append([filename] + cmd_readcell + cmd_readfc)
                _iparam += 1

        # maybe shall we put an input specific for `nac type`?
        if "nac_params" in self.inputs:
            cmds_with_nac = []
            for params in calculation_cmds:
                cmds_with_nac.append(params + "--nac")
            calculation_cmds = cmds_with_nac

        # ================= retreiving phonopy output files ===================

        # These two are always parsed
        retrieve_temporary_list = [
            self._INOUT_FORCE_CONSTANTS,
            self.inputs.metadata.options.output_filename,
        ]

        # !!! TEMPORARY
        # !!! THEY SHOULD DEPEND ON THE INPUT PARAMETERS
        # !!! NEED TO IMPLEMENT ERRORS TOO
        for value in self._OUTPUTS.values():
            retrieve_temporary_list.append(value)

        # ============================ calcinfo ===============================

        local_copy_list = []
        remote_copy_list = []

        # Copy remote output dir
        parent_calc_folder = self.inputs.get("parent_folder", None)
        if isinstance(parent_calc_folder, orm.RemoteData):
            remote_copy_list.append(
                (
                    parent_calc_folder.computer.uuid,
                    os.path.join(parent_calc_folder.get_remote_path(), self._INPUT_SUBFOLDER),
                    self._OUTPUT_SUBFOLDER,
                )
            )
        elif isinstance(parent_calc_folder, orm.FolderData):
            for filename in parent_calc_folder.list_object_names():
                local_copy_list.append(
                    (parent_calc_folder.uuid, filename, os.path.join(self._OUTPUT_SUBFOLDER, filename))
                )

        calcinfo = datastructures.CalcInfo()
        calcinfo.uuid = self.uuid
        calcinfo.local_copy_list = local_copy_list  # what we write in the folder
        calcinfo.remote_copy_list = remote_copy_list  # what to copy from remote folder
        calcinfo.retrieve_list = []  # what to retrieve and keep - we do not keep anything
        calcinfo.retrieve_temporary_list = retrieve_temporary_list  # what to retrieve and then parse

        # We add more then one code info, i.e. more than one cmd line
        calcinfo.codes_info = []
        for cmdline_params in calculation_cmds:
            codeinfo = datastructures.CodeInfo()
            codeinfo.cmdline_params = cmdline_params
            codeinfo.code_uuid = self.inputs.code.uuid
            codeinfo.withmpi = self.inputs.metadata.options.withmpi
            calcinfo.codes_info.append(codeinfo)

        return calcinfo

    def _write_born(self, folder):
        """Write in `folder` the `BORN` file."""
        if "nac_params" in self.inputs and "primitive" in self.inputs:
            born_txt = get_BORN_txt(self.inputs.nac_params, self.inputs.primitive, self.inputs.symmetry_tolerance)

            with folder.open(self._INPUT_NAC, "w", encoding="utf8") as handle:
                handle.write(born_txt)

    def _write_cells_info(self, folder):
        """Write in `folder` the `phonopy_cells.yaml` file."""
        phpy_yaml_txt = get_phonopy_yaml_txt(
            self.inputs.structure,
            symprec=self.inputs.symmetry_tolerance.value,
            supercell_matrix=self.inputs["supercell_matrix"],
        )
        with folder.open(self._DEFAULT_CELL_FILE, "w", encoding="utf8") as handle:
            handle.write(phpy_yaml_txt)

    def _write_force_sets(self, folder):
        """Write in `folder` the `FORCE_SETS` file."""
        if "force_sets" in self.inputs:
            force_sets = self.inputs.force_sets
        else:
            force_sets = None

        if "displacement_dataset" in self.inputs:
            dataset = self.inputs["displacement_dataset"].get_dict()
        elif "dataset" in self.inputs:
            dataset = self.input["dataset"]
        elif "dataset" in self.inputs and "displacements" in self.inputs.dataset.get_arraynames():
            dataset = {"displacements": self.inputs.dataset.get_array("displacements")}
            if "forces" in self.inputs.dataset.get_arraynames():
                dataset["forces"] = self.inputs.dataset.get_array("forces")
                if force_sets is not None:
                    force_sets = None
        else:
            dataset = None

        # can work both for type-I and type-II
        force_sets_txt = get_FORCE_SETS_txt(dataset, force_sets=force_sets)
        if force_sets_txt is None:
            raise InputValidationError("Displacements or forces were not found.")

        with folder.open(self._INPUT_FORCE_SETS, "w", encoding="utf8") as handle:
            handle.write(force_sets_txt)

    @staticmethod
    def _write_conf_info(folder, parameters: dict, filename: str):
        """Write in `folder` the `phonopy.conf` file(s)."""
        characters = []

        for key, value in parameters.items():
            characters.append(key)
            characters.append("=")
            if isinstance(value, str):
                characters.append(value)
            elif isinstance(value, (float, int)):
                characters.append(str(value))
            elif isinstance(value, list):
                string_value = ""
                for val in value:
                    string_value += str(val) + " "
                characters.append(string_value[:-1])
            elif isinstance(value, bool):
                characters.append(f".{value}.")
            else:
                raise NotImplementedError(f"value type `{type(value)}` is not supported.")
            characters.append("\n")  # one line for each key in parameters

        with folder.open(filename, "w", encoding="utf8") as handle:
            for character in characters:
                handle.write(character)

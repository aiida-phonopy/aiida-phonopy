from aiida.parsers.parser import Parser
from aiida_phonopy.common.raw_parsers import (
    parse_thermal_properties, parse_FORCE_CONSTANTS, parse_partial_DOS,
    parse_band_structure)


class PhonopyParser(Parser):
    """
    Parser the DATA files from phonopy.
    """

    def __init__(self, calc):
        """
        Initialize the instance of PhonopyParser
        """
        super(PhonopyParser, self).__init__(calc)

    def parse(self, **kwargs):
        """
        Parses the datafolder, stores results.
        """

        # suppose at the start that the job is successful
        successful = True

        # select the folder object
        # Check that the retrieved folder is there
        try:
            output_folder = self.retrieved
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = output_folder.get_folder_list()

        # OUTPUT file should exist
        # if not self._calc._OUTPUT_FILE_NAME in list_of_files:
        #    successful = False
        #    self.logger.error("Output file not found")
        #    return successful, ()

        # Get files and do the parsing

        # save the outputs
        new_nodes_list = []

        if self._calc._INOUT_FORCE_CONSTANTS in list_of_files:
            with output_folder.open(self._calc._INOUT_FORCE_CONSTANTS) as f:
                self.out('force_constants', parse_FORCE_CONSTANTS(f))

        if self._calc._OUTPUT_DOS in list_of_files:
            with output_folder.open(self._calc._OUTPUT_DOS) as f:
                dos_object = parse_partial_DOS(
                    f, self.node.inputs.structure, self.node.inputs.parameters)
                self.out('dos', dos_object)

        if self._calc._OUTPUT_THERMAL_PROPERTIES in list_of_files:
            with output_folder.open(self._calc._OUTPUT_THERMAL_PROPERTIES) as f:
                self.out('thermal_properties', parse_thermal_properties(f))

        if self._calc._OUTPUT_BAND_STRUCTURE in list_of_files:
            with output_folder.open(self._calc._OUTPUT_BAND_STRUCTURE) as f:
                bs_object = parse_band_structure(f, self.node.inputs.bands)
                self.out('band_structure', bs_object)

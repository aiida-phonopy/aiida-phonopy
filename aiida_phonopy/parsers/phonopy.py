from aiida.engine import ExitCode
from aiida.common.exceptions import NotExistent
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
        self.logger.info("Parsing start.")

        # select the folder object
        # Check that the retrieved folder is there
        try:
            output_folder = self.retrieved
        except NotExistent:
            return self.exit_codes.ERROR_NO_RETRIEVED_FOLDER

        # check what is inside the folder
        list_of_files = output_folder.list_object_names()

        # OUTPUT file should exist
        # if not self._calc._OUTPUT_FILE_NAME in list_of_files:
        #    successful = False
        #    self.logger.error("Output file not found")
        #    return successful, ()

        # Get files and do the parsing

        fc_filename = self.node.inputs.force_constants_filename.value
        if fc_filename in list_of_files:
            with output_folder.open(fc_filename) as f:
                fname = f.name
            self.out('force_constants', parse_FORCE_CONSTANTS(fname))

        pdos_filename = self.node.inputs.projected_dos_filename.value
        if pdos_filename in list_of_files:
            with output_folder.open(pdos_filename) as f:
                fname = f.name
            pdos_object = parse_partial_DOS(
                fname, self.node.inputs.structure, self.node.inputs.parameters)
            self.out('dos', pdos_object)

        tp_filename = self.node.inputs.thermal_properties_filename.value
        if tp_filename in list_of_files:
            with output_folder.open(tp_filename) as f:
                fname = f.name
            self.out('thermal_properties', parse_thermal_properties(fname))

        band_filename = self.node.inputs.band_structure_filename.value
        if band_filename in list_of_files:
            with output_folder.open(band_filename) as f:
                fname = f.name
            bs_object = parse_band_structure(fname, self.node.inputs.bands)
            self.out('band_structure', bs_object)

        self.logger.info("Parsing done.")
        return ExitCode(0)

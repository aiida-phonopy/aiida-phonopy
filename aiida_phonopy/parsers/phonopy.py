from aiida.engine import ExitCode
from aiida.common.exceptions import NotExistent
from aiida.parsers.parser import Parser
from aiida_phonopy.common.raw_parsers import (
    parse_thermal_properties, parse_FORCE_CONSTANTS, parse_projected_dos,
    parse_total_dos, parse_band_structure)


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

        projected_dos_filename = self.node.inputs.projected_dos_filename.value
        if projected_dos_filename in list_of_files:
            with output_folder.open(projected_dos_filename) as f:
                fname = f.name
            self.out('pdos', parse_projected_dos(fname))

        total_dos_filename = self.node.inputs.projected_dos_filename.value
        if total_dos_filename in list_of_files:
            with output_folder.open(total_dos_filename) as f:
                fname = f.name
            self.out('dos', parse_total_dos(fname))

        tp_filename = self.node.inputs.thermal_properties_filename.value
        if tp_filename in list_of_files:
            with output_folder.open(tp_filename) as f:
                fname = f.name
            self.out('thermal_properties', parse_thermal_properties(fname))

        sym_dataset = self.node.inputs.settings['symmetry']
        label = "%s (%d)" % (sym_dataset['international'],
                             sym_dataset['number'])
        band_filename = self.node.inputs.band_structure_filename.value
        if band_filename in list_of_files:
            with output_folder.open(band_filename) as f:
                fname = f.name
            self.out('band_structure',
                     parse_band_structure(fname, label=label))

        self.logger.info("Parsing done.")
        return ExitCode(0)

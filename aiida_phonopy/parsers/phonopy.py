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
            out_folder = retrieved[self._calc._get_linkname_retrieved()]
        except KeyError:
            self.logger.error("No retrieved folder found")
            return False, ()

        # check what is inside the folder
        list_of_files = out_folder.get_folder_list()

        # OUTPUT file should exist
        # if not self._calc._OUTPUT_FILE_NAME in list_of_files:
        #    successful = False
        #    self.logger.error("Output file not found")
        #    return successful, ()

        # Get files and do the parsing

        # save the outputs
        new_nodes_list = []

        if self._calc._INOUT_FORCE_CONSTANTS in list_of_files:
            outfile = out_folder.get_abs_path(
                self._calc._INOUT_FORCE_CONSTANTS)
            object_force_constants = parse_FORCE_CONSTANTS(outfile)
            new_nodes_list.append(('force_constants', object_force_constants))

        if self._calc._OUTPUT_DOS in list_of_files:
            outfile = out_folder.get_abs_path(self._calc._OUTPUT_DOS)
            dos_object = parse_partial_DOS(outfile, self._calc.inp.structure,
                                           self._calc.inp.parameters)
            new_nodes_list.append(('dos', dos_object))

        if self._calc._OUTPUT_THERMAL_PROPERTIES in list_of_files:
            outfile = out_folder.get_abs_path(
                self._calc._OUTPUT_THERMAL_PROPERTIES)
            tp_object = parse_thermal_properties(outfile)
            new_nodes_list.append(('thermal_properties', tp_object))

        if self._calc._OUTPUT_BAND_STRUCTURE in list_of_files:
            outfile = out_folder.get_abs_path(
                self._calc._OUTPUT_BAND_STRUCTURE)
            bs_object = parse_band_structure(outfile, self._calc.inp.bands)
            new_nodes_list.append(('band_structure', bs_object))

        # look at warnings
        with open(out_folder.get_abs_path(self._calc._SCHED_ERROR_FILE)) as f:
            errors = f.readlines()
        if errors:
            for error in errors:
                self.logger.warning(error)
                successful = False

        return successful, new_nodes_list

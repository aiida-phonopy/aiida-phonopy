from aiida.parsers.parser import Parser
from aiida.parsers.exceptions import OutputParsingError
from aiida_phonopy.common.raw_parsers import parse_kappa


class Phono3pyParser(Parser):
    """
    Parser the DATA files from phonopy.
    """

    def __init__(self, calc):
        """
        Initialize the instance of PhonopyParser
        """
        super(Phono3pyParser, self).__init__(calc)

    def parse_with_retrieved(self, retrieved):
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

        #if self._calc._OUTPUT_KAPPA in list_of_files:
        for filename in filter(lambda x: x.startswith(self._calc._OUTPUT_KAPPA), list_of_files):
            outfile = out_folder.get_abs_path(filename)
            kappa_data = parse_kappa(outfile)
            new_nodes_list.append(('kappa', kappa_data))

        # look at warnings
        with open(out_folder.get_abs_path(self._calc._SCHED_ERROR_FILE)) as f:
            errors = f.readlines()
        if errors:
            for error in errors:
                self.logger.warning(error)
                successful = False

        return successful, new_nodes_list

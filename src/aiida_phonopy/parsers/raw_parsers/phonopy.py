# -*- coding: utf-8 -*-
"""Raw parsers of the phonopy output files."""

from aiida_phonopy.utils.mapping import get_logging_container


def parse_stdout(stdout):
    """Raw parser of the phonopy std output.

    :param stdout: str of the std output of phonopy
    :returns: tuple of two dictionaries, with the parsed data and log messages, respectively
    """
    stdout_lines = stdout.splitlines()
    logs = get_logging_container()

    parsed_data = {}
    is_incomplete = True

    for line in stdout_lines:
        if 'Python' in line:
            parsed_data.update({'python_version': line.split()[-1]})

        if 'Spglib' in line:
            parsed_data.update({'spglib_version': line.split()[-1]})

        if 'One of the following run modes may be specified for phonon calculations.' in line:
            logs.error.append('ERROR_BAD_INPUTS')

        if 'Summary of calculation was written in "phonopy.yaml".' in line:
            is_incomplete = False

    if is_incomplete:
        logs.error.append('ERROR_OUTPUT_STDOUT_INCOMPLETE')

    return parsed_data, logs

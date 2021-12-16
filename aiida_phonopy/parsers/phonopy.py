# -*- coding: utf-8 -*-
"""Parsers of `PhonopyPpCalculation` output files."""

import os

from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory

from .raw_parsers import parse_FORCE_CONSTANTS, parse_yaml, parse_total_dos, parse_projected_dos

PhonopyPpCalculation = CalculationFactory("phonopy.pp")


class PhonopyPpParser(Parser):
    """Parser the files produced by a phonopy post processing calculation."""

    def parse(self, **kwargs):
        """Parse retrieved files from remote folder."""
        retrieved = self.retrieved
        retrieve_temporary_list = self.node.get_attribute("retrieve_temporary_list", None)

        # If temporary files were specified, check that we have them
        if retrieve_temporary_list:
            try:
                retrieved_temporary_folder = kwargs["retrieved_temporary_folder"]
            except KeyError:
                return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER)

        # How to get the output filenames and how to open them, depends on whether they will have been retrieved in the
        # `retrieved` output node, or in the `retrieved_temporary_folder`. Instead of having a conditional with almost
        # the same loop logic in each branch, we apply a somewhat dirty trick to define an `opener` which is a callable
        # that will open a handle to the output file given a certain filename. This works since it is guaranteed that
        # these output files (excluding the standard output) will all either be in the retrieved, or in the retrieved
        # temporary folder. (Courtesy: from aiida-quantumespresso parsers)
        if retrieve_temporary_list:
            filenames = os.listdir(retrieved_temporary_folder)
            file_opener = lambda filename: open(os.path.join(retrieved_temporary_folder, filename))
        else:
            filenames = retrieved.list_object_names()
            file_opener = retrieved.open

        # We check first if `force_constants.hdf5` and `phonopy.yaml` are not missing in the retrieved files.
        # They always must be present, since they are always computed in the post-processing calcualtion.
        try:
            filename = PhonopyPpCalculation._INOUT_FORCE_CONSTANTS
            filenames.remove(filename)
        except ValueError:
            return self.exit_codes.ERROR_NO_FORCE_CONSTANTS
        # need to read in binary - this might complicate optional retrieve option
        with open(os.path.join(retrieved_temporary_folder, filename), "rb") as f:
            self.out("force_constants", parse_FORCE_CONSTANTS(f))

        try:
            filename = PhonopyPpCalculation._DEFAULT_OUTPUT_FILE
            filenames.remove(filename)
        except ValueError:
            return self.exit_codes.ERROR_NO_PHONOPY_YAML
        with file_opener(filename) as f:  # there should be "rb"
            self.out("parameters", parse_yaml(f))

        # self.node.inputs.parameters.get_dict()
        # can be used to check the if the wanted retrieved files are present in folder
        # not trivial

        projected_dos_filename = PhonopyPpCalculation._OUTPUTS["pdos"]
        if projected_dos_filename in filenames:
            with file_opener(projected_dos_filename) as f:
                self.out("pdos", parse_projected_dos(f))

        total_dos_filename = PhonopyPpCalculation._OUTPUTS["dos"]
        if total_dos_filename in filenames:
            with file_opener(total_dos_filename) as f:
                self.out("dos", parse_total_dos(f))

        tp_filename = PhonopyPpCalculation._OUTPUTS["tprop"]
        if tp_filename in filenames:
            with file_opener(tp_filename) as f:
                self.out("thermal_properties", parse_yaml(f))

        band_filename = PhonopyPpCalculation._OUTPUTS["band"]
        if band_filename in filenames:
            with file_opener(band_filename) as f:
                self.out("band", parse_yaml(f))

        # since most of them are yaml files it will be better to loop the parsing procedure
        # while exluding eventual dat files (or any of different extension).

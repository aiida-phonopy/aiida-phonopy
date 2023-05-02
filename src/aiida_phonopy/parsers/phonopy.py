# -*- coding: utf-8 -*-
"""Parsers of `PhonopyCalculation` output files."""

import os
import traceback

from aiida import orm
import h5py
import numpy as np
import yaml

from aiida_phonopy.calculations.phonopy import PhonopyCalculation
from aiida_phonopy.utils.mapping import _uppercase_dict, get_logging_container

from .base import Parser


def file_opener(folder_path, filename):
    """Function to open sistematically each expected output format."""
    filename_path = os.path.join(folder_path, filename)

    if filename.endswith('hdf5'):
        return h5py.File(filename_path, 'r')

    return open(filename_path, mode='r', encoding='utf-8')


class PhonopyParser(Parser):
    """Parser the files produced by a phonopy post processing calculation."""

    def parse(self, **kwargs):
        """Parse retrieved files from remote folder."""
        retrieved = self.retrieved
        retrieve_temporary_list = self.node.base.attributes.get('retrieve_temporary_list', None)

        # If temporary files were specified, check that we have them
        if retrieve_temporary_list:
            try:
                retrieved_temporary_folder = kwargs['retrieved_temporary_folder']
            except KeyError:
                return self.exit(self.exit_codes.ERROR_NO_RETRIEVED_TEMPORARY_FOLDER)

            # !FIXME!
            # We should implement `filenames` depending whether they were
            # put in the retrieved or in the temporary folder.
            filenames = os.listdir(retrieved_temporary_folder)

        # Parse the stdout
        parsed_stdout, logs_stdout, exit_code_stdout = self.parse_stdout()
        self.emit_logs(logs_stdout)

        if exit_code_stdout:
            return self.exit(exit_code_stdout)

        # Parse phonopy.yaml
        try:
            filename = PhonopyCalculation._DEFAULT_PHONOPY_FILE
            filenames.remove(filename)

            with file_opener(retrieved_temporary_folder, filename) as file:
                try:
                    parsed_phonopy_yaml = yaml.safe_load(file)
                except yaml.YAMLError:
                    return self.exit(self.exit_codes.ERROR_OUTPUT_YAML_LOAD)
        except ValueError:
            try:
                with retrieved.base.repository.open(filename) as file:
                    try:
                        parsed_phonopy_yaml = yaml.safe_load(file)
                    except yaml.YAMLError:
                        return self.exit(self.exit_codes.ERROR_OUTPUT_YAML_LOAD)
            except IOError:
                # This should in principle never happen!
                return self.exit(self.exit_codes.ERROR_OUTPUT_PHONOPY_MISSING)

        # Output parameters
        phonopy_yaml_keys = ['phonopy', 'physical_unit', 'space_group']
        parsed_parameters = {key: value for key, value in parsed_phonopy_yaml.items() if key in phonopy_yaml_keys}
        parsed_parameters.update(parsed_stdout)
        self.out('output_parameters', orm.Dict(parsed_parameters))

        # We use the input `parameters` to check the expected retrieved files in folder.
        expected_filenames_keys = self.get_expected_filenames_keys()

        parsers = {
            'fc': self.parse_force_constants,
            'dos': self.parse_total_dos,
            'pdos': self.parse_projected_dos,
            'band': self.parse_band_structure,
            'mesh': self.parse_qpoints,  # to do a proper one
            'qpoints': self.parse_qpoints,  # to do a proper one
            'irreps': self.parse_yaml,  # to do a proper one
            'tprop': self.parse_thermal_properties,
            'tdisp': self.parse_yaml,  # to do a proper one
            'tdispmat': self.parse_yaml,  # to do a proper one
        }

        missing_files = []

        # Loop over every expected output file, parse it and expose its outputs.
        # If one or more are missing, we put it in the `missing_files` and raise error afterwards.
        for expected_key in expected_filenames_keys:
            try:
                filename, filout = PhonopyCalculation._OUTPUTS[expected_key]
                filenames.remove(filename)
                if expected_key == 'fc':
                    fc_output = self.parse_force_constants(os.path.join(retrieved_temporary_folder, filename))
                    self.out(filout, fc_output)
                else:
                    with file_opener(retrieved_temporary_folder, filename) as file:
                        self.out(filout, parsers[expected_key](file))
            except (FileNotFoundError, ValueError):
                missing_files.append(filename)

        if missing_files:
            return self.exit(self.exit_codes.ERROR_OUTPUT_FILES_MISSING)

        # We don't parse animation files.
        for filename in filenames:
            if filename.endswith(('jmol', 'xyz', 'xyz_jmol', 'arc', 'ascii')) or filename.startswith('APOSCAR-*'):
                filenames.remove(filename)

        # !FIXME!
        # Still to work out better. Good chances are that is not necessary...
        #
        # If something is still present and it has not been parsed, we raise error.
        # if filenames:
        #     self.exit_codes.ERROR_NOT_IMPLEMENTED_PARSER

    def parse_stdout(self):
        """Parse the stdout output file.

        :param parameters: the input parameters dictionary
        :param parsed_phonopy: the collected parsed data from the yaml phonopy output
        :return: raw parsed data
        """
        from .raw_parsers.phonopy import parse_stdout as parse_stdout_

        parsed_data = {}
        exit_code_stdout = None

        filename_stdout = self.node.get_option('output_filename')
        logs = get_logging_container()  # datastructure for logs

        if filename_stdout not in self.retrieved.base.repository.list_object_names():
            exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_MISSING
            return parsed_data, logs, exit_code_stdout

        try:
            stdout = self.retrieved.base.repository.get_object_content(filename_stdout)
        except (IOError, OSError):
            exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_READ
            return parsed_data, logs, exit_code_stdout

        try:
            # Check for job completion, indicating that phonopy exited without interruption, even if there was an error.
            parsed_data, logs = parse_stdout_(stdout)
        except ValueError:
            logs.critical.append(traceback.format_exc())
            exit_code_stdout = self.exit_codes.ERROR_UNEXPECTED_PARSER_EXCEPTION

        # If the stdout was incomplete, most likely the job was interrupted before it could cleanly finish, so the
        # output files are most likely corrupt and cannot be restarted
        if 'ERROR_OUTPUT_STDOUT_INCOMPLETE' in logs['error']:
            exit_code_stdout = self.exit_codes.ERROR_OUTPUT_STDOUT_INCOMPLETE

        if 'ERROR_BAD_INPUTS' in logs['error']:
            parameters = _uppercase_dict(self.node.inputs.parameters.get_dict(), dict_name='parameters')
            if 'FORCE_CONSTANTS' in parameters:
                if str.upper(parameters['FORCE_CONSTANTS']) != 'WRITE':
                    exit_code_stdout = self.exit_codes.ERROR_BAD_INPUTS

        return parsed_data, logs, exit_code_stdout

    def get_expected_filenames_keys(self):
        """Return the retrieve file keys (that map to the filenames outputs) depending on the tags in `parameters`."""
        retrieved_set = set()
        parameters = _uppercase_dict(self.node.inputs.parameters.get_dict(), dict_name='parameters')

        maps = {
            'band': {'BAND', 'BAND_PATHS', 'BAND_POINTS', 'BAND_LABELS', 'BAND_CONNECTION', 'BAND_INDICES'},
            'mesh': {'MESH', 'MP', 'MESH_NUMBERS', 'MP_SHIFT', 'GAMMA_CENTER', 'WRITE_MESH'},
            'dos': {'DOS'},
            'pdos': {'PDOS'},
            'tprop': {'TPROP'},
            'tdisp': {'TDISP'},
            'tdispmat': {'TDISPMAT'},
            'qpoints': {'QPOINTS', 'WRITEDM'},
            'fc': {'WRITE_FORCE_CONSTANTS', 'FORCE_CONSTANTS'},
            'mod': {'MODULATION'},
            'irreps': {'IRREPS', 'SHOW_IRREPS', 'LITTLE_COGROUP'}
        }

        for tag, value in parameters.items():
            for key_map, tag_map in maps.items():
                if tag in tag_map:
                    # Usually these tags are set only if one needs them, i.e. to TRUE,
                    # thus this further control is actually pointless. Nevertheless,
                    # we want to avoid to raise error when in inputs one of these
                    # tags is set to FALSE, meaning there will be no output to parse.
                    if tag.startswith(('WRITE', 'DOS', 'TPROP', 'TDISP', 'TDISPMAT')):
                        if value:
                            retrieved_set.add(key_map)
                    else:
                        if tag != 'FORCE_CONSTANTS':
                            retrieved_set.add(key_map)
                        elif str.upper(value) == 'WRITE':
                            retrieved_set.add(key_map)

        # Double check on mesh. The writing tag must have priority.
        if 'WRITE_MESH' in parameters.keys():
            if not parameters['WRITE_MESH']:
                retrieved_set.remove('mesh')

        return retrieved_set

    def parse_force_constants(self, filepath):
        """Parse the `force_constants.hdf5` output file."""
        from phonopy.file_IO import read_force_constants_hdf5

        p2s_map = self._get_p2s_map()
        force_constants = read_force_constants_hdf5(filename=filepath, p2s_map=p2s_map)

        fc_array = orm.ArrayData()
        fc_array.set_array('force_constants', force_constants)
        fc_array.label = 'force_constants'

        return fc_array

    def load_with_numpy(self, file):
        """Load a txt file using numpy."""
        try:
            data = np.loadtxt(file)
        except ValueError:
            self.exit(self.exit_codes.ERROR_OUTPUT_NUMPY_LOAD)
        return data

    def load_with_yaml(self, file):
        """Load a yaml file using."""
        try:
            data = yaml.safe_load(file)
        except yaml.YAMLError:
            self.exit(self.exit_codes.ERROR_OUTPUT_YAML_LOAD)
        return data

    def parse_yaml(self, file):
        """Parse a `.yaml` file and return it as a Dict."""
        data = self.load_with_yaml(file=file)
        return orm.Dict(data)

    def parse_total_dos(self, file):
        """Parse `total_dos.dat` output file."""
        data = self.load_with_numpy(file=file)

        total_dos = {'frequency_points': data[:, 0], 'total_dos': data[:, 1]}
        dos = orm.XyData()
        dos.set_x(total_dos['frequency_points'], 'Frequency', 'THz')
        dos.set_y(total_dos['total_dos'], 'DOS', '1/THz')
        dos.label = 'Total DOS'

        return dos

    def parse_projected_dos(self, file):
        """Parse `projected_dos.dat` output file."""
        data = self.load_with_numpy(file=file)

        projected_dos = {'frequency_points': data[:, 0], 'projected_dos': data[:, 1:].T}
        pdos = orm.XyData()
        pdos_list = list(projected_dos['projected_dos'])
        pdos.set_x(projected_dos['frequency_points'], 'Frequency', 'THz')
        pdos.set_y(
            pdos_list,
            [
                'Projected DOS',
            ] * len(pdos_list),
            [
                '1/THz',
            ] * len(pdos_list),
        )
        pdos.label = 'Projected DOS'

        return pdos

    def parse_thermal_properties(self, file):
        """Parse the `thermal_properties.yaml`` output file."""
        data = self.load_with_yaml(file=file)

        thermal_properties = {
            'temperatures': [],
            'free_energy': [],
            'entropy': [],
            'heat_capacity': [],
        }

        for tp in data['thermal_properties']:
            thermal_properties['temperatures'].append(tp['temperature'])
            thermal_properties['entropy'].append(tp['entropy'])
            thermal_properties['free_energy'].append(tp['free_energy'])
            thermal_properties['heat_capacity'].append(tp['heat_capacity'])
        old_thermal_property = thermal_properties.copy()
        for key, item in old_thermal_property.items():
            thermal_properties[key] = np.array(item)

        tprops = orm.XyData()
        tprops.set_x(thermal_properties['temperatures'], 'Temperature', 'K')
        tprops.set_y(
            [
                thermal_properties['free_energy'],
                thermal_properties['entropy'],
                thermal_properties['heat_capacity'],
            ],
            ['Helmholtz free energy', 'Entropy', 'Cv'],
            ['kJ/mol', 'J/K/mol', 'J/K/mol'],
        )
        tprops.label = 'Thermal properties'

        return tprops

    def parse_band_structure(self, file, freqs_units='THz'):
        """Parse the `band.hdf5`` output file.

        Expected keys are:
            * **nqpoint**: array with total number of q-points in the band structure, array(1,)
            * **frequency**: array of frequencies at each q-point; array(npath, segment_nqpoint, nband)
            * **label**: array of labels, two per each path; array(npath, 2)
            * **path**: number of q-points per path; array(npath, segment_nqpoint, qpoint/qposition(3,))
            * **distance**: distance between consecutive q-points in a segment of path; array(npath, segment_nqpoint)
            * **segment_nqpoint**: number of q-points per segment of the entire path; array(npath,)
            * **eigenvector**: (optional) eigenvector of each phonon mode, i.e. each frequency of each q-point;
                array(npath, segment_nqpoint, nband, eigenvector(6,))
            * **group_velocity**: (optional) group velocity of each phonon mode (as for eigenvector);
                array(npath, segment_nqpoint, nband, group_velocity(3,))
        """
        import re

        # Loading the data
        data = {key: file[key][:] for key in file.keys()}

        # Initializing the bands object and setting the reciprocal lattice from the unitcell.
        if 'phonopy_data' in self.node.inputs:
            structure = self.node.inputs['phonopy_data'].get_unitcell()
        else:
            structure = self.node.inputs['force_constants'].get_unitcell()

        band_structure = orm.BandsData()
        band_structure.set_cell_from_structure(structuredata=structure)

        # Extracting the frequencies and qpoints.
        frequencies = []
        qpoints = []
        for freqs_per_path, qpoint_per_path in zip(data['frequency'], data['path']):
            for freqs_per_qpoint in freqs_per_path:
                frequencies.append(freqs_per_qpoint)
            for qpoint in qpoint_per_path:
                qpoints.append(qpoint)

        # Extracting the labels along with the index of the qpoint.
        labels = []
        nqpoint = [0, -1]
        for label_pair, segment in zip(
            data['label'], data['segment_nqpoint']
        ):  # they are saved in pairs, two per segment
            nqpoint[1] += segment
            for i, label in enumerate(label_pair):
                try:  # the labels are saved in dtype |S16 --> e.g. b'$\\mathrm{X}$'
                    symbol = re.search(r'\{(\w+)\}', str(label)).group(1)
                except AttributeError:  # it means we hit Gamma
                    if label == b'$\\Gamma$':
                        symbol = r'$\Gamma$'
                    else:
                        symbol = '?'
                labels.append([nqpoint[i], symbol])
            nqpoint[0] += segment

        # Here I remove the symilar symbols if they are consecutive.
        # The sorting could be smarter, but it works.
        final_labels = []
        for nq, pair in enumerate(labels):
            if not nq == 0:
                if not pair[1] == labels[nq - 1][1]:
                    final_labels.append(pair)
            else:
                final_labels.append(pair)

        # Setting the final bands structure data.
        band_structure.set_kpoints(np.array(qpoints))
        band_structure.set_bands(np.array(frequencies), units=freqs_units)
        band_structure.labels = final_labels

        return band_structure

    def parse_qpoints(self, file, freqs_units='THz'):
        """Parse the `mesh.hdf5`` and `qpoints.hdf5`` output files.

        Expected keys are:
            * **frequency**: array of frequencies at each q-point; array(nqpoint, nband)
            * **mesh**: qpoint mesh; array(3,)
            * **qpoint**: qpoints; array(nqpoint, 3)
            * **weight**: weight of qpoints; array(nqpoint,)
            * **eigenvector**: (optional) eigenvector of each phonon mode, i.e. each frequency of each q-point;
                array(npath, segment_nqpoint, nband, eigenvector(6,))
            * **group_velocity**: (optional) group velocity of each phonon mode (as for eigenvector);
                array(npath, segment_nqpoint, nband, group_velocity(3,))
        """
        # Loading the data
        data = {key: file[key][:] for key in file.keys()}

        # Initializing the bands object and setting the reciprocal lattice from the unitcell.
        if 'phonopy_data' in self.node.inputs:
            structure = self.node.inputs['phonopy_data'].get_unitcell()
        else:
            structure = self.node.inputs['force_constants'].get_unitcell()

        band_structure = orm.BandsData()
        band_structure.set_cell_from_structure(structuredata=structure)

        qpoints = np.array(data['qpoint'])
        frequencies = np.array(data['frequency'])
        if file.filename.endswith('qpoints.hdf5'):
            weights = None
        else:
            weights = np.array(data['weight'].astype(float))

        # Setting the final bands structure data.
        band_structure.set_kpoints(qpoints, weights=weights)
        band_structure.set_bands(frequencies, units=freqs_units)

        return band_structure

    def _get_p2s_map(self):
        """Get the primitive to supercell map."""
        if 'force_constants' in self.node.inputs:
            return self.node.inputs.force_constants.get_cells_mappings()['primitive']['p2s_map']
        if 'phonopy_data' in self.node.inputs:
            return self.node.inputs.phonopy_data.get_cells_mappings()['primitive']['p2s_map']
        return None

from aiida.orm import Data


class BandStructureData(Data):
    """
    Store the band structure.
    """

    def __init__(self, *args, **kwargs):
        super(BandStructureData, self).__init__(*args, **kwargs)

    def get_number_of_bands(self):

        if 'nbands' in self.get_attrs():
            return self.get_attr('nbands')
        else:
            return None

    def get_number_of_points(self):

        if 'npoints' in self.get_attrs():
            return self.get_attr('npoints')
        else:
            return None


    def set_bands(self, ranges):

        import tempfile
        import numpy

        ranges = numpy.array(ranges)

        with tempfile.NamedTemporaryFile() as f:
            numpy.save(f, ranges)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, 'bands.npy')

        self._set_attr('nbands', len(ranges))
        self._set_attr('npoints', len(ranges[0]))

    def set_labels(self, band_labels):

        import tempfile
        import numpy

        band_labels = numpy.array(band_labels)

        with tempfile.NamedTemporaryFile() as f:
            numpy.save(f, band_labels)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, 'band_labels.npy')

    def set_unitcell(self, unitcell):
        """
        :param unitcell: Numpy Array that contains the unitcell matrix (Lattice vectors in columns)
        :return:
        """

        import tempfile
        import numpy

        unitcell = numpy.array(unitcell)

        with tempfile.NamedTemporaryFile() as f:
            numpy.save(f, unitcell)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, 'unitcell.npy')

    def set_band_structure_phonopy(self, band_structure_phonopy):

        import tempfile
        import numpy

        q_points = numpy.array(band_structure_phonopy[0])
        distances = numpy.array(band_structure_phonopy[1])

        # Check consistency
        if self.get_bands() is None:
            self.set_bands(q_points)
        else:
            numpy.testing.assert_array_almost_equal(q_points, self.get_bands(), decimal=4)

        numpy.testing.assert_array_almost_equal(distances, self.get_distances(), decimal=4)

        for element in {
                        # 'q_points.npy': numpy.array(band_structure_phonopy[0]),
                        # 'distances.npy': numpy.array(band_structure_phonopy[1]),
                        'frequencies.npy':numpy.array(band_structure_phonopy[2]),
                        }.items():

            with tempfile.NamedTemporaryFile() as f:
                numpy.save(f, element[1])
                f.flush()  # Important to flush here, otherwise the next copy command
                # will just copy an empty file
                self.add_path(f.name, element[0])

    def set_band_structure_gruneisen(self, band_structure_gruneisen):

        import tempfile
        import numpy

        q_points = numpy.array([band[0] for band in band_structure_gruneisen._paths])
        distances = numpy.array([band[5] for band in band_structure_gruneisen._paths])

        # Check consistency
        if self.get_bands() is None:
            self.set_bands(q_points)
        else:
            numpy.testing.assert_array_almost_equal(q_points, self.get_bands(), decimal=4)

        numpy.testing.assert_array_almost_equal(distances, self.get_distances(), decimal=4)

        for element in {'gamma.npy': numpy.array([band[2] for band in band_structure_gruneisen._paths]),
                        'eigenvalues.npy': numpy.array([band[3] for band in band_structure_gruneisen._paths]),
                        'frequencies.npy': numpy.array([band[4] for band in band_structure_gruneisen._paths]),
                        }.items():

            with tempfile.NamedTemporaryFile() as f:
                numpy.save(f, element[1])
                f.flush()  # Important to flush here, otherwise the next copy command
                # will just copy an empty file
                self.add_path(f.name, element[0])


    def get_unitcell(self):
        """
        Return the unitcell in the node as a numpy array
        """
        import numpy

        fname = 'unitcell.npy'
        if fname not in self.get_folder_list():
            return None

        array = numpy.load(self.get_abs_path(fname))
        return array

    def get_distances(self, band=None):
        """
        Return the distances in the node as a numpy array
        """
        import numpy as np

        inverse_unitcell = np.linalg.inv(self.get_unitcell())

        array = self.get_bands()
        # nbands = array.shape[0]
        npoints = array.shape[1]

        distances = []
        for piece in array:
            if distances == []:
                band_dist = [0.0]
            else:
                band_dist = [distances[-1][-1]]
            for i in range(npoints-1):
                #band_dist.append(np.linalg.norm(np.array(band[i+1]) - np.array(band[i]))+band_dist[i])
                band_dist.append(np.linalg.norm(np.dot(piece[i+1], inverse_unitcell.T) -
                                                np.dot(piece[i], inverse_unitcell.T)) + band_dist[i])

            distances.append(band_dist)
        distances = np.array(distances)

        if band is not None:
            distances = distances[band]

        return distances


    def get_frequencies(self, band=None):
        """
        Return the frequencies in the node as a numpy array
        """
        import numpy

        fname = 'frequencies.npy'
        if fname not in self.get_folder_list():
            return None

        array = numpy.load(self.get_abs_path(fname))
        if band is not None:
            array = array[band]

        return array

    def get_gamma(self, band=None):
        """
        Return the frequencies in the node as a numpy array
        """
        import numpy

        fname = 'gamma.npy'
        if fname not in self.get_folder_list():
            return None

        array = numpy.load(self.get_abs_path(fname))

        if band is not None:
            array = array[band]

        return array

    def get_eigenvalues(self, band=None):
        """
        Return the frequencies in the node as a numpy array
        """
        import numpy

        fname = 'eigenvalues.npy'
        if fname not in self.get_folder_list():
            return None


        array = numpy.load(self.get_abs_path(fname))
        if band is not None:
            array = array[band]

        return array

    def get_bands(self, band=None):
        """
        Return the bands in the node as a numpy array
        """
        import numpy

        fname = 'bands.npy'
        if fname not in self.get_folder_list():
            return None

        array = numpy.load(self.get_abs_path(fname))
        if band is not None:
            array = array[band]


        return array

    def get_band_ranges(self, band=None):
        """
        Return the bands in the node as a numpy array
        """
        import numpy

        array = self.get_bands()

        band_ranges = numpy.array([numpy.array([i[-0], i[-1]]) for i in array])

        if band is not None:
            band_ranges = band_ranges[band]

        return band_ranges


    def get_labels(self, band=None):
        """
        Return the band labels in the node as a numpy array
        """
        import numpy

        fname = 'band_labels.npy'

        if fname not in self.get_folder_list():
            return None

        array = numpy.load(self.get_abs_path(fname))

        if band is not None:
            array = array[band]

        return array

    def get_plot_helpers(self, style='latex'):

        # Collection of helpers to plot band data

        if style == 'unicode':

            substitutions = {'GAMMA': u'\u0393'
                            }

        elif style == 'latex':
            substitutions = {'GAMMA': u'$\Gamma$'
                             }
        else:
            return Exception('label_style not supported')

        def replace_list(text_string, substitutions):
            for item in substitutions.iteritems():
                text_string = text_string.replace(item[0], item[1])
            return text_string

        distances = self.get_distances()
        labels_array = self.get_labels()

        labels = []
        indices = []
        block = [replace_list(labels_array[0][0], substitutions)]
        block_indices = [0]

        for i, freq in enumerate(distances[:-1]):
            if labels_array[i+1][0] == labels_array[i][1]:
                block.append(replace_list(labels_array[i+1][0], substitutions))
                block_indices.append(i+1)
            else:
                block.append(replace_list(labels_array[i][1], substitutions))
                labels.append(block)
                indices.append(block_indices)

                block = [replace_list(labels_array[i+1][0], substitutions)]
                block_indices = [i+1]

        block.append(replace_list(labels_array[-1][1], substitutions))
        labels.append(block)
        indices.append(block_indices)

        widths = []
        ranges = []
        positions = []
        for j, index in enumerate(indices):
            widths.append(self.get_distances(band=index[-1])[-1] - self.get_distances(band=index[0])[0])
            ranges.append([self.get_distances(band=index[0])[0], self.get_distances(band=index[-1])[-1]])
            positions.append([self.get_distances(band=i)[0] for i in index] + [self.get_distances(band=index[-1])[-1]])

        return labels, indices, widths, ranges, positions


    def get_formatted_labels_matplotlib(self):
        distances = self.get_distances()
        labels_array = self.get_labels()

        substitutions = {'GAMMA': u'\u0393'
                         }

        substitutions = {'GAMMA': u'$\Gamma$'
                         }

        def replace_list(text_string, substitutions):

            for item in substitutions.iteritems():
                text_string = text_string.replace(item[0], item[1])

            return text_string

        labels = []
        labels_positions = []
        for i, freq in enumerate(distances):
            if labels_array[i][0] == labels_array[i-1][1]:
                labels.append(replace_list(labels_array[i][0],substitutions))
            else:
                labels.append(replace_list(labels_array[i-1][1]+'/'+labels_array[i][0], substitutions))
            labels_positions.append(distances[i][0])
        labels_positions.append(distances[-1][-1])
        labels.append(replace_list(labels_array[-1][1], substitutions))
        labels[0] = replace_list(labels_array[0][0], substitutions)

        return labels, labels_positions
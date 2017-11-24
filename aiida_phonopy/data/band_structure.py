from aiida.orm.data.array import ArrayData


class BandStructureData(ArrayData):
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

    def set_bands(self, bands):

        import numpy

        bands = numpy.array(bands)

        self.set_array('bands', bands)

        self._set_attr('nbands', len(bands))
        self._set_attr('npoints', len(bands[0]))

    def set_labels(self, band_labels):

        import numpy

        band_labels = numpy.array(band_labels)

        self.set_array('band_labels', band_labels)

    def set_unitcell(self, unitcell):
        """
        :param unitcell: Numpy Array that contains the unitcell matrix (Lattice vectors in columns)
        :return:
        """

        import numpy

        unitcell = numpy.array(unitcell)

        self.set_array('unitcell', unitcell)

    def set_band_structure_phonopy(self, band_structure_phonopy):

        import numpy

        q_points = numpy.array(band_structure_phonopy[0])
        distances = numpy.array(band_structure_phonopy[1])

        # Check consistency
        if self.get_bands() is None:
            self.set_bands(q_points)
        else:
            numpy.testing.assert_array_almost_equal(q_points, self.get_bands(), decimal=4)

        numpy.testing.assert_array_almost_equal(distances, self.get_distances(), decimal=4)

        # self.set_array('q_points', numpy.array(band_structure_phonopy[0]))
        # self.set_array('distances', numpy.array(band_structure_phonopy[1]))
        self.set_array('frequencies', numpy.array(band_structure_phonopy[2]))

    def set_band_structure_gruneisen(self, band_structure_gruneisen):

        import numpy

        q_points = numpy.array([band[0] for band in band_structure_gruneisen._paths])
        distances = numpy.array([band[5] for band in band_structure_gruneisen._paths])

        # Check consistency
        if self.get_bands() is None:
            self.set_bands(q_points)
        else:
            numpy.testing.assert_array_almost_equal(q_points, self.get_bands(), decimal=4)

        numpy.testing.assert_array_almost_equal(distances, self.get_distances(), decimal=4)

        self.set_array('gamma', numpy.array([band[2] for band in band_structure_gruneisen._paths]))
        self.set_array('eigenvalues', numpy.array([band[3] for band in band_structure_gruneisen._paths]))
        self.set_array('frequencies', numpy.array([band[4] for band in band_structure_gruneisen._paths]))

    def set_frequencies(self, frequencies):
        """
        Return the frequencies in the node as a numpy array
        """

        import numpy

        frequencies = numpy.array(frequencies)

        self.set_array('frequencies', frequencies)

    def get_unitcell(self):
        """
        Return the unitcell in the node as a numpy array
        """

        return self.get_array('unitcell')

    def get_distances(self, band=None):
        """
        Return the distances between q-points calculated from bands
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

        frequencies = self.get_array('frequencies')

        if band is not None:
            frequencies = frequencies[band]

        return frequencies

    def get_gamma(self, band=None):
        """
        Return the frequencies in the node as a numpy array
        """

        gamma =  self.get_array('gamma')

        if band is not None:
            gamma = gamma[band]

        return gamma

    def get_eigenvalues(self, band=None):
        """
        Return the frequencies in the node as a numpy array
        """
        eigenvalues = self.get_array('eigenvalues')

        if band is not None:
            eigenvalues = eigenvalues[band]

        return eigenvalues

    def get_bands(self, band=None):
        """
        Return the bands in the node as a numpy array
        """
        bands = self.get_array('bands')

        if band is not None:
            bands = bands[band]

        return bands

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

        band_labels = self.get_array('band_labels')

        if band is not None:
            band_labels = band_labels[band]

        return band_labels

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
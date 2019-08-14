from aiida.orm.nodes.data.array import ArrayData


class PhononDosData(ArrayData):
    """
    Store the phonon DOS on disk as a numpy array. It requires numpy to be installed.
    """

    def __init__(self, *args, **kwargs):
        super(PhononDosData, self).__init__(*args, **kwargs)
        self._cached_arrays = {}

    def _get_equivalent_atom_list(self, with_weights=False):
        import numpy
        partial_dos = self.get_array('partial_dos')
        partial_symbols = self.get_attr("atom_labels")

        # Check atom equivalences
        delete_list = []
        weight_list = [1] * len(partial_dos)

        for i, dos_i in enumerate(partial_dos):
            for j, dos_j in enumerate(partial_dos):
                if i < j:
                    if numpy.allclose(dos_i, dos_j, rtol=1, atol=1e-8) and partial_symbols[i] == partial_symbols[j]:
                        delete_list.append(j)
                        weight_list[i] += 1

        if with_weights:
            return (numpy.delete(range(len(partial_dos)), delete_list),
                    numpy.delete(weight_list, delete_list))

        return numpy.delete(range(len(partial_dos)), delete_list)

    def get_dos(self):
        """
        Return the density of states in a 1D array (according to get_frequencies())
        """

        return self.get_array('dos')

    def get_number_of_partial_dos(self, full=False):
        """
        Return the number of atoms (equivalent to the dimension of partial density of states)

        :param full: if true, return the total number of atoms. By default return the number of
        non-equivalent atoms by symmetry. This option is in agreement with get_partial_dos()
        :return: integer
        """

        partial_dos = self.get_array('partial_dos')

        if full:
            return len(partial_dos)

        return len(partial_dos[self._get_equivalent_atom_list()])

    def get_partial_dos(self, full=False):
        """
        Return the partial density of states

        :param full: if true, return the partial DOS for all atoms in the unitcell. By default only the
        non-equivalent atoms by symmetry are returned
        :return: 2D numpy array with the partial density of states
        """

        partial_dos = self.get_array('partial_dos')

        if full:
            return partial_dos

        indices, weights = self._get_equivalent_atom_list(with_weights=True)
        partial_dos = [partial * weight for weight, partial in zip(weights, partial_dos[indices])]

        return partial_dos

    def get_frequencies(self):
        """
        Return the corresponding frequencies to get_dos() and get_partial_dos()
        """

        return self.get_array('frequencies')

    def get_atom_labels(self, full=False):
        """
        Store the atomic symbols as a numpy array.

        :param full: if True, return the atomic symbols of all atoms. By default (false) return the atomic symbols
        of non-equivalent atoms by symmetry
        """
        import numpy

        labels = self.get_attr("atom_labels")
        if full:
            return labels
        return numpy.array(labels)[self._get_equivalent_atom_list()].tolist()

    def set_atom_labels(self, labels):
        """
        Store the phonon dos as a numpy array

        :param array: The numpy array to store
        """
        self.set_attribute("atom_labels", labels)

    def set_dos(self, array):
        """
        Store the phonon dos as a numpy array
        :param array: The numpy array to store
        """

        self.set_array('dos', array)

    def set_frequencies(self, array):
        """
        Store the frequencies as a numpy array

        :param array: The numpy array to store
        """
        self.set_array('frequencies', array)

    def set_partial_dos(self, array):
        """
        Store the partial dos as a numpy array

        :param array: The numpy array to store
        """

        self.set_array('partial_dos', array)
        self.set_attribute("n_partial_dos", len(array))

    def is_stable(self, tol=1e-6):
        """
        Returns true if no imaginary modes (shown as negative frequency DOS)
        :param tol:
        :return:
        """
        import numpy as np

        dos = self.get_dos()
        freq = self.get_frequencies()

        mask_neg = np.ma.masked_less(freq, 0.0).mask
        mask_pos = np.ma.masked_greater(freq, 0.0).mask

        if mask_neg.any() == False:
            return True

        if mask_pos.any() == False:
            return False

        int_neg = -np.trapz(np.multiply(dos[mask_neg], freq[mask_neg]), x=freq[mask_neg])
        int_pos = np.trapz(np.multiply(dos[mask_pos], freq[mask_pos]), x=freq[mask_pos])

        if int_neg / int_pos > tol:
            return False
        else:
            return True


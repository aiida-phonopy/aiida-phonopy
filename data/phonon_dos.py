from aiida.orm import Data


class PhononDosData(Data):
    """
    Store the phonon DOS on disk as a numpy array. It requires numpy to be installed.
    """

    def __init__(self, *args, **kwargs):
        super(PhononDosData, self).__init__(*args, **kwargs)
        self._cached_arrays = {}

    def _get_equivalent_atom_list(self):
        import numpy
        fname = 'partial_dos.npy'
        partial_dos = numpy.load(self.get_abs_path(fname))
        partial_symbols = self.get_attr("atom_labels")

        # Check atom equivalences
        list = range(len(partial_dos))
        delete_list = []
        for i, dos_i in enumerate(partial_dos):
            for j, dos_j in enumerate(partial_dos):
                if i < j:
                    if numpy.allclose(dos_i, dos_j, rtol=1, atol=1e-8) and partial_symbols[i] == partial_symbols[j]:
                        dos_i += dos_j
                        delete_list.append(j)

        return numpy.delete(range(len(partial_dos)), delete_list)

    def get_dos(self):
        """
        Return the force constants stored in the node as a numpy array
        """
        import numpy

        fname = 'dos.npy'

        array = numpy.load(self.get_abs_path(fname))
        return array

    def get_number_of_partial_dos(self, full=False):
        """
        Return the force constants stored in the node as a numpy array
        """
        import numpy

        fname = 'partial_dos.npy'

        array = numpy.load(self.get_abs_path(fname))

        if full:
            return len(array)

        return len(array[self._get_equivalent_atom_list()])

    def get_partial_dos(self, full=False):
        """
        Return the force constants stored in the node as a numpy array
        """
        import numpy

        fname = 'partial_dos.npy'

        array = numpy.load(self.get_abs_path(fname))

        if full:
            return array
        return array[self._get_equivalent_atom_list()]

    def get_frequencies(self):
        """
        Return the frequencies stored in the node as a numpy array
        """
        import numpy

        fname = 'frequencies.npy'

        array = numpy.load(self.get_abs_path(fname))

        return array


    def get_atom_labels(self, full=False):
        """
        Store the phonon dos as a numpy array.
        :param array: The numpy array to store.
        """
        import numpy

        labels = self.get_attr("atom_labels")

        if full:
            return labels
        return numpy.array(labels)[self._get_equivalent_atom_list()].tolist()

    def set_atom_labels(self, labels):
        """
        Store the phonon dos as a numpy array.
        :param array: The numpy array to store.
        """
        self._set_attr("atom_labels", labels)



    def set_dos(self, array):
        """
        Store the phonon dos as a numpy array.
        :param array: The numpy array to store.
        """

        import tempfile
        import numpy

        fname = "dos.npy"
        with tempfile.NamedTemporaryFile() as f:
            # Store in a temporary file, and then add to the node
            numpy.save(f, array)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, fname)


    def set_frequencies(self, array):
        """
        Store the frequencies as a numpy array.
        :param array: The numpy array to store.
        """

        import tempfile
        import numpy

        fname = "frequencies.npy"
        with tempfile.NamedTemporaryFile() as f:
            # Store in a temporary file, and then add to the node
            numpy.save(f, array)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, fname)


    def set_partial_dos(self, array):
        """
        Store the partial dos as a numpy array.
        :param array: The numpy array to store.
        """

        import tempfile
        import numpy

        fname = "partial_dos.npy"
        with tempfile.NamedTemporaryFile() as f:
            # Store in a temporary file, and then add to the node
            numpy.save(f, array)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, fname)

        self._set_attr("n_partial_dos", len(array))

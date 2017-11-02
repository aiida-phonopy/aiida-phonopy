from aiida.orm import Data


class ForceConstantsData(Data):
    """
    Store the force constants on disk as a numpy array. It requires numpy to be installed.
    """

    def __init__(self, *args, **kwargs):
        super(ForceConstantsData, self).__init__(*args, **kwargs)
        self._cached_arrays = {}

    def get_shape(self):
        """
        Return the shape of an array (read from the value cached in the
        properties for efficiency reasons).
        :param name: The name of the array.
        """
        return self.get_attr("shape")

    def get_array(self):
        """
        Return the force constants stored in the node as a numpy array
        """
        import numpy

        fname = 'force_constants.npy'

        array = numpy.load(self.get_abs_path(fname))
        return array

    def set_array(self, array):
        """
        Store the force constants as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.
        :param array: The numpy array to store.
        """

        import tempfile

        import numpy

        if not (isinstance(array, numpy.ndarray)):
            raise TypeError("ArrayData can only store numpy arrays. Convert "
                            "the object to an array first")

        fname = "force_constants.npy"

        with tempfile.NamedTemporaryFile() as f:
            # Store in a temporary file, and then add to the node
            numpy.save(f, array)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, fname)

        self._set_attr("shape", list(array.shape))


    def get_epsilon(self):
        """
        Return dielectric tensor stored in the node as a numpy array
        """
        import numpy

        fname = 'epsilon.npy'

        if fname not in self.get_folder_list():
            return None

        array = numpy.load(self.get_abs_path(fname))

        return array

    def get_born_charges(self):
        """
        Return born charges stored in the node as a numpy array
        """
        import numpy

        fname = 'born_charges.npy'

        if fname not in self.get_folder_list():
            return None

        array = numpy.load(self.get_abs_path(fname))

        return array

    def epsilon_and_born_exist(self):

        """
        Check if born charges and epsion exists
        """

        return self.get_epsilon() is not None and self.get_born_charges() is not None

    def set_born_charges(self, array):
        """
        Store Born charges as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.
        :param array: The numpy array to store.
        """

        import tempfile
        import numpy

        fname = "born_charges.npy"

        array = numpy.array(array)
        with tempfile.NamedTemporaryFile() as f:
            # Store in a temporary file, and then add to the node
            numpy.save(f, array)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, fname)

    def set_epsilon(self, array):
        """
        Store the dielectric tensor as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.
        :param array: The numpy array to store.
        """

        import tempfile
        import numpy

        fname = "epsilon.npy"

        array = numpy.array(array)
        with tempfile.NamedTemporaryFile() as f:
            # Store in a temporary file, and then add to the node
            numpy.save(f, array)
            f.flush()  # Important to flush here, otherwise the next copy command
            # will just copy an empty file
            self.add_path(f.name, fname)

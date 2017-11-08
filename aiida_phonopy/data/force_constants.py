from aiida.orm.data.array import ArrayData


class ForceConstantsData(ArrayData):
    """
    Store the force constants on disk as a numpy array. It requires numpy to be installed.
    """

    def __init__(self, *args, **kwargs):
        super(ForceConstantsData, self).__init__(*args, **kwargs)
        self._cached_arrays = {}

    def get_data(self):
        """
        Return the force constants stored in the node as a numpy array
        """

        return self.get_array('force_constants')

    def set_data(self, force_constants):
        """
        Store the force constants as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.
        :param array: The numpy array to store.
        """
        import numpy

        self.set_array('force_constants', numpy.array(force_constants))

    def get_epsilon(self):
        """
        Return dielectric tensor stored in the node as a numpy array
        """

        return self.get_array('epsilon')

    def get_born_charges(self):
        """
        Return born charges stored in the node as a numpy array
        """
        return self.get_array('born_charges')

    def epsilon_and_born_exist(self):

        """
        Check if born charges and epsion exists
        """

        return 'born_charges' in self.get_arraynames() and 'epsilon' in self.get_arraynames()

    def set_born_charges(self, born_charges):
        """
        Store Born charges as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.
        :param array: The numpy array to store.
        """
        import numpy

        self.set_array('_born_charges', numpy.array(born_charges))

    def set_epsilon(self, epsilon):
        """
        Store the dielectric tensor as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.
        :param array: The numpy array to store.
        """
        import numpy

        self.set_array('_born_charges', numpy.array(epsilon))

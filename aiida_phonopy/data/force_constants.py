from aiida.orm.data.array import ArrayData
import numpy

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

        self.set_array('force_constants', numpy.array(force_constants))

    def read_from_phonopy_file(self, filename):
        """
        Read the force constants from a phonopy FORCE_CONSTANTS file
        :param filename: FORCE_CONSTANTS file name
        """

        fcfile = open(filename)
        num = int((fcfile.readline().strip().split())[0])
        force_constants = numpy.zeros((num, num, 3, 3), dtype=float)
        for i in range(num):
            for j in range(num):
                fcfile.readline()
                tensor = []
                for k in range(3):
                    tensor.append([float(x) for x in fcfile.readline().strip().split()])
                force_constants[i, j] = numpy.array(tensor)

        self.set_data(force_constants)

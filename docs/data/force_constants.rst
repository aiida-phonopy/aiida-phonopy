Force constants
===============

This object contains the second order force constants as a numpy array.
It also stores the dielectric tensor and Born effective charges to be used in the calculation
of the non analytical corrections in the phonon calculation. get_epsilon

Setters and getters are provided to store and get the data in phonopy format:

.. automodule:: aiida_phonopy.data.force_constants
.. autoclass:: ForceConstantsData()
   :members: set_data, get_data, read_from_phonopy_file



Force constants
===============

This object contains the harmonic force constansts as a numpy array using the phonopy format.
It also stores the dielectric tensor and Born effective charges to be used in the calculation
of the non analytical corrections in the phonon calculation. get_epsilon

Setters and getters are provided to store and get the data in phonopy format:

.. automodule:: data.force_constants
.. autoclass:: ForceConstantsData()
   :members: set_array, get_array, set_epsilon, get_epsilon, get_born_charges



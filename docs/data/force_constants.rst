Force constants
===============

This object contains the second order force constants as a numpy array.

Setters and getters are provided to store and get the data in phonopy format:

.. automodule:: aiida_phonopy.data.force_constants
.. autoclass:: ForceConstantsData()
   :members: set_data, get_data, read_from_phonopy_file

example of use
--------------
::

    force_constants = phonon.get_force_constants()
    force_constants_data = ForceConstantsData(data=force_constants)

    ...

    phonon.set_force_constants(force_constants_data.get_data())

Band structure
==============
This object contains the phonon band structure calculated by phonopy.
Also it can store the mode Gruneisen parameters data.

Special setters are provided to store the data directly from phonopy objects:

.. automodule:: aiida_phonopy.data.band_structure
.. autoclass:: BandStructureData()
   :members: set_bands, set_labels, set_unitcell, set_band_structure_phonopy, set_band_structure_gruneisen


Also special getters are provided to return data in phonopy format:

.. autoclass:: BandStructureData()
   :members: get_bands, get_band_ranges, get_distances, get_labels, get_unitcell, get_frequencies, get_gamma, get_eigenvalues


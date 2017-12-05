Phonon density of states
========================

This object contains the phonon band structure calculated by phonopy.
Also stores the mode Gruneisen parameters data in Gruneisen calculation.

Special setters are provided to store the data directly from phonopy objects:

.. automodule:: aiida_phonopy.data.dos
.. autoclass:: BandStructureData()
   :members: get_dos, get_partial_dos, get_frequencies, get_atom_labels, set_band_structure_gruneisen


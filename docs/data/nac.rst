Non-analytical corrections (NAC)
================================

This object contains the information needed to calculate the non-analytical corrections
in phonopy phonon calculations. This includes the crystal structure, Born effective charges
and dielectric tensor.

Setters and getters are provided to store the data from aiida objects and get the data in phonopy format:

.. automodule:: aiida_phonopy.data.nac
.. autoclass:: NacData()
   :members: set_structure, set_born_charges, set_epsilon, get_born_parameters_phonopy



example of use
--------------
Setup the nac_data object from calculation output data ::

    structure = StructureData()
    born_charges = NumpyArray()
    dielectric_tensor = NumpyArray()

    nac_data = NacData(structure=structure,
                       born_charges=born_charges,
                       epsilon=dielectric_tensor)

Use this object to setup a phonopy calculation ::

    phonon = Phonopy()
    primitive = phonon.get_primitive()
    nac_parameters = nac_data.get_born_parameters_phonopy(primitive_cell=primitive.get_cell())
    phonon.set_nac_params(nac_parameters)

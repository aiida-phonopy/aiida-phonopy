Phonopy
=======

This plugin is designed to calculate the harmonic force constants using phonopy.

.. function:: PhonopyCalculation(structure, parameters, data_sets, nac_data)

   :param structure: AiiDA StructureData object that contains the crystal unit cell information
   :param parameters: AiiDA ParametersData data object that the information about the requested machine resources
   :param data_sets: AiiDA ForceSetsData  object that contains the phonopy input parameters

The outputs of this plugin are:

* **force_constants**: ForceConstantsData object that contains the harmonic force constants.

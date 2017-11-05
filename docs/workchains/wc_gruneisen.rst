
wc_gruneisen
============

This workchain performs a mode Gruneisen parameters calculation using Phonopy. This workchain is designed to
keep the compatibility with wc_phonon workchain (and other phonon related workchains deveoped in the future).
For this reason the input parameters in the workchain have the same structure that in wc_phonon.
This workchain performs 3 phonon calculations. One with the crystal structure optimized at 0 kB external pressure
and the other two with the unit cell optimized at a little higher and lower pressure obtaining a slightly
smaller and larger unit cell respectively. This pressure deference can be set as an optional argument, by default
its value is 1e-2 kB.

.. function:: GruneisenPhonopy(structure, machine, ph_settings, es_settings [, stress_displacement=1e-2])

   :param structure: AiiDA StructureData object that contains the crystal unit cell information
   :param machine: AiiDA ParametersData data object that the information about the requested machine resources
   :param ph_settings: AiiDA ParametersData data  object that contains the phonopy input parameters
   :param es_settings: AiiDA ParameterData object that the electronic structure calculation input parameters. These parameters depends on the code used (see workchains/launcher examples)
   :param stress_displacement: AiiDA FloatData object. This determines the pressure difference between the 3 phonon calculations.


The outputs of this workchain are:

* **band_structure**: BandStructure object that contains the phonon band structure and the mode Gruneisen parameters.
* **mesh**: ArrayData object that contains the q-points sampling mesh of the mode Gruneises parameters.

In workchains/tools/plot_gruneisen.py there is an example of how to extract the information from these outputs.
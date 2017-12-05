Gruneisen
=========

This WorkChain performs a mode Gruneisen parameters calculation using Phonopy. This workchain is designed to
keep the compatibility with phonon WorkChain input structure. For this reason the input parameters in the workchain
have the same structure. Please, check phonon Workchain documentation for detailed information.
This workchain calculates the mode gruneisen parameters by performing 3 phonon calculations.
One with the crystal structure optimized at a given external stress (default: 0 kB) and the other two
with the unit cell optimized at a little higher and lower stress (defined by stress_displacement)
obtaining a slightly smaller and larger unit cell respectively.
Stress_displacement can be set as an optional argument, by default its value is 1e-2 kB.

.. function:: GruneisenPhonopy(structure, machine, ph_settings, es_settings [, stress_displacement=1e-2])

   :param structure: AiiDA StructureData object that contains the crystal unit cell information
   :param ph_settings: AiiDA ParametersData data  object that contains the phonopy input parameters
   :param es_settings: AiiDA ParameterData object that the electronic structure calculation input parameters.
   :param pressure: AiiDA FloatData object. This determines the absolute stress (in kBar) at which the reference crystal structure is optimized (default 0).
   :param stress_displacement: AiiDA FloatData object. This determines the stress difference between the 3 phonon calculations.

The outputs of this workchain are:

* **band_structure**: BandStructure object that contains the phonon band structure and the mode Gruneisen parameters.
* **mesh**: ArrayData object that contains the wave vectors sampling mesh and the mode Gruneises parameters at each wave vector.

Each one of this objects have its own methods for extracting the information. Check the individual object documentation
for more details. workchains/tools/plot_gruneisen.py contains a complete example script showing how to extract the information from these outputs.

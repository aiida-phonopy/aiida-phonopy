Phonon3
=======

This WorkChain calculated the 3rd order force constants using phono3py.
This WorkChain requires one of the plugins for VASP, QuantumESPRESSO and LAMMPS described in phonon WorkChain.
Non-analytical corrections can be calculated from the Born effective charges and dielectric tensor which
are only implemented for VASP plugin.

.. function:: PhononPhono3py(structure, ph_settings, es_settings [, optimize=True, use_nac=False, pressure= 0.0, calculate_fc=False])

   :param structure: AiiDA StructureData object that contains the crystal unit cell structure.
   :param ph_settings: AiiDA ParametersData data object that contains the phonopy input parameters.
   :param es_settings: AiiDA ParameterData object that contains the calculator input parameters. These parameters depends on the code used (see workchains/launcher examples)
   :param use_nac: (optional) AiiDA BooleanData object. Determines if non-analytical corrections will be included in the phonon calculations. By default this option is False.
   :param optimize: (optional) AiiDA BooleanData object. Determines if a crystal unit cell optimization is performed or not before the phonon calculation. By default this option is True.
   :param pressure: (optional) AiiDA FloatData object. If optimize is True, this sets the external pressure (in kB) at which the unit cell optimization is preformed. By default this option takes value 0 kB.
   :param calculate_fc: (optional) AiiDA BooleanData object. Determines if the 2on and 3rd order force constants are calculated. By default this option is False.
   :param chunks: (optional) AiiDA Int object that defines the maximum number of calculation to submit simultaneously. The next set of calculation will not be submitted until the previous set is finished.
   :param data_sets: (optional) AiiDA ForceSets object that contains the forces and displacements of a previously calculation. This data_set can be the output of either phonon3 or phonon WorkChains.


The results outputs of this WorkChain are the following :

* **data_sets**: ForceSetsData object that contains the information of supercells with displacements and forces. Check
* **force_constants_2order**: ForceConstantsData object that contains the 2ond order force constants.
* **force_constants_3order**: ForceConstantsData object that contains the 3rd order force constants.
* **final_structure**: StructureData containing the optimized structure.

The ForceSetsData object obtained as a output of phonon3 WorkChain can be also in harmonic phonon calculation without modification.

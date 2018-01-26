Thermal_conductivity
====================

This WorkChain is designed to calculate the thermal conductivity with phono3py using an interative
procedure by increasing the cutoff distance. This WorkChain makes use of phonon WorkChain to calculate
the harmonic band structure and DOS and determine the dynamical stability of the crystal structure.
If the crystal structure is stable then it calculates the thermal conductivity using phono3py using
a small cutoff distance (defined by the user) and successively increase it until reaching convergence.
The convergence criteria are given by the user as atol and rtol parameters. The convergence criteria works
in the same manner of numpy allclose function (https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.allclose.html),
applied to all values of kappa vs temperature.

.. function:: ThermalPhono3py(structure, ph_settings, es_settings [, optimize=True, use_nac=False, pressure= 0.0, calculate_fc=False, gp_chunks=1, gp_chunks=10, initial_cutoff=2.0, step=1.0, atol=0.1, rtol=0.3])

   :param structure: AiiDA StructureData object that contains the crystal unit cell structure.
   :param ph_settings: AiiDA ParametersData data object that contains the phonopy input parameters.
   :param es_settings: AiiDA ParameterData object that contains the calculator input parameters. These parameters depends on the code used (see workchains/launcher examples)
   :param use_nac: AiiDA BooleanData object. Determines if non-analytical corrections will be included in the phonon calculations. By default this option is False.
   :param optimize: AiiDA BooleanData object. Determines if a crystal unit cell optimization is performed or not before the phonon calculation. By default this option is True.
   :param pressure: AiiDA FloatData object. If optimize is True, this sets the external pressure (in kB) at which the unit cell optimization is preformed. By default this option takes value 0 kB.
   :param chunks: AiiDA Int object that defines the maximum number of calculation to submit simultaneously. The next set of calculation will not be submitted until the previous set is finished.
   :param gp_chunks:  AiiDA Int object that defines the number of computers in which the calculation will be distributed (default: 1).
   :param initial_cutoff: AiiDA Float object that defines initial cutoff distance for phono3py calculation.
   :param step: AiiDA Float object that defines increment of cutoff distance in each iteration.
   :param atol: AiiDA Float object that defines the convergence criteria of thermal conductivity, absolute value (thermal conductivity units).
   :param rtol: AiiDA Float object that defines the convergence criteria of thermal conductivity, relative value.

The results outputs of this WorkChain are the following :

* **kappa**: ArrayData object that contains the results stored in kappa-mxxx-gx.hdf5 output file.
* **final_structure**: StructureData object that contains the optimized unit cell. If no optimization is performed this contains the same StructureData object provided as a input.
* **thermal_properties**: ArrayData object that contains the harmonic thermal properties calculated using phonopy. These include the entropy, free energy and heat capacity at constant volume.
* **dos**: PhononDosData object that contains the harmonic phonon full and partial density of states.
* **band_structure**: BandStructureData object that contains the harmonic phonon band structure.

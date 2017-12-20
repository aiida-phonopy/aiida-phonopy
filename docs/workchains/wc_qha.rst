Quasi-harmonic approximation
============================

This WorkChain performs a quasi-harmonic approximation calculation using Phonopy.
This WorkChain requires one of the plugins for VASP, QuantumESPRESSO and LAMMPS described in phonon WorkChain.
Non-analytical corrections can be included from the Born effective charges and dielectric tensor which
are only implemented for VASP plugin.
Phonopy can be used either locally and remotely. To use phonopy remotely phonopy code must be setup as described
in AiiDA documentation (https://aiida-core.readthedocs.io/en/latest/get_started/index.html#code-setup-and-configuration).
using the phonopy plugin provided in this package.

.. function:: QHAPhonopy(structure, ph_settings, es_settings [, optimize=True, use_nac=False, num_expansions=10])

   :param structure: AiiDA StructureData object that contains the crystal unit cell structure.
   :param ph_settings: AiiDA ParametersData data  object that contains the phonopy input parameters.
   :param es_settings: AiiDA ParameterData object that contains the calculator input parameters. These parameters depends on the code used (see workchains/launcher examples)
   :param num_expansions (optional): AiiDA IntData object. The number of volume expansions around the optimized structure at zero pressure to perform. By default the value is 10.
   :param use_nac (optional): AiiDA BooleanData object. Determines if non-analytical corrections will be included in the phonon calculations. By default this option is False.

The results outputs of this WorkChain are the following :

* **qha_results**: ArrayData object that contains the thermal properties at constant pressure calculated using phonopy.
This ArrayData includes the following arrays: qha_temperatures, helmholtz_volume, thermal_expansion, volume_temperature, heat_capacity_P_numerical,
volume_expansion and gibbs_temperature.
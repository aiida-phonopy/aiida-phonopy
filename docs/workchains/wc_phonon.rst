
wc_phonon
=========

This performs a phonon calculation using Phonopy. At this moment, to calculate the forces
3 codes can be used: Quantum Espresso, VASP and LAMMPS. Quantum Espresso uses the plugin shipped with
AiiDA. VASP requires a modified version of the plugin developed by Mario Zic (you can find it here:
https://github.com/abelcarreras/aiida_extensions). LAMMPS requires the on development plugin found here:
https://github.com/abelcarreras/aiida-lammps . At this moment Born effective charges can only be calculated
using VASP. Phonopy by default is used locally in workfunctions, however the force constants can be calculated
remotely in cluster computers using the phonopy plugin (see examples in workchains/launchers).

.. function:: PhononPhonopy(structure, machine, ph_settings, es_settings [, optimize=True, pressure=0.0])

   :param structure: AiiDA StructureData object that contains the crystal unit cell information
   :param machine: AiiDA ParametersData data object that the information about the requested machine resources
   :param ph_settings: AiiDA ParametersData data  object that contains the phonopy input parameters
   :param es_settings: AiiDA ParameterData object that the electronic structure calculation input parameters. These parameters depends on the code used (see workchains/launcher examples)
   :param optimize: AiiDA booleanData object. Determines if a crystal unit cell optimization is performed or not before the phonon calculation
   :param pressure: AiiDA FloatData object. If optimize is True, this sets the external pressure (in kB) at which the unit cell optimization is preformed.


The outputs of this workchain are:

* **force_constants**: ForceConstantsData object that contains the harmonic force constants. If Born effective charges are calculated this object also contains the dielectric tensor and the born effective charges of each atom in the unit cell
* **thermal_properties**: ArrayData object that contains the thermal properties calculated using phonopy. These include the entropy, free energy and heat capacity at constant volume.
* **dos**: PhononDosData object that contains the phonon full and partial density of states.
* **band_structure**: BandStructureData object that contains the harmonic phonon band structure.
* **final_structure**: StructureData object that contains the optimized unit cell. If no optimization is performed this is the same StructureData object provided as a input.

In workchains/tools/plot_phonon.py there is an example of how to extract the information from these outputs.
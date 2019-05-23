Phonopy
=======

This plugin is designed to calculate the harmonic force constants, thermal properties (entropy, free energy
 and heat capacity at constant volume) and density of states using phonopy.

.. function:: PhonopyCalculation(structure, parameters, data_sets, nac_data)

   :param structure: StructureData object that contains the crystal unit cell information
   :param parameters: ParametersData data object that contains the phonopy input parameters
   :param data_sets: ForceSetsData object that contains the forces, directions and detail of all the supercells with displacements (equivalent to FORCE_SETS file in phonopy)
   :param force_constants: AiiDA ForceConstants object that contains the force constants
   :param nac_data: (optional) NacData object that contains the Born effective charges and the dielectric tensor

- input parameters should be dictionary with the following entries ::

    parameters_dict = {'supercell': [[2, 0, 0],
                                     [0, 2, 0],
                                     [0, 0, 2]],
                       'primitive': [[1.0, 0.0, 0.0],
                                     [0.0, 1.0, 0.0],
                                     [0.0, 0.0, 1.0]],
                       'distance': 0.01,
                       'mesh': [40, 40, 40],
                       'symmetry_tolerance': 1e-5}
    ParameterData = DataFactory('dict')
    parameters = ParameterData(dict=parameters_dict)

Either data_sets of force_constants should be used. If data_sets is used force constants will be calculated
and returned as a calculation output. If force_constants is used the calculation will be faster.

- data_sets can be created from a phonopy FORCE_SETS file by ::

    ForceSetsData = DataFactory('phonopy.force_sets')
    force_sets = ForceSetsData()
    force_sets.read_from_phonopy_file('FORCE_SETS')

- force_constants can be create from a phonopy FORCE_CONSTANTS file by ::

    ForceConstantsData = DataFactory('phonopy.force_constants')
    force_constants = ForceConstantsData()
    force_constants.read_from_phonopy_file('FORCE_CONSTANTS')

- nac_data is an optional parameter and can be created from single point calculation, the required information is the crystal structure(StructureData) used in the calculation, born effective charges (numpy array) for all the atoms in the crystal structure and the dielectric tensor (numpy array [dim: Natoms x 3 x 3]) ::

    nac_data = NacData(structure=crystal_structure,
                       born_charges=born_charges_numpy_array,
                       epsilon=epsilon_numy_array)

The outputs of this plugin are:

* **force_constants**: ForceConstantsData object that contains the second order force constants.
* **dos**: DosData object that contains the full and partial density of states.
* **band_structure**: BandStructureData object that contains the phonon band structure.
* **thermal_properties**: ArrayData object that contains the entropy, free energy and heat capacity at constant volume.


Take a look at the examples in examples/plugins folder for reference

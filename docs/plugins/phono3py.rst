Phono3py
========

This plugin is designed to calculate the thermal conductivity using phono3py.

.. function:: Phono3pyCalculation(structure, parameters, data_sets, nac_data)

   :param structure: AiiDA StructureData object that contains the crystal unit cell information
   :param parameters: AiiDA ParametersData data object that contains the phonopy input parameters
   :param data_sets: AiiDA ForceSetsData object that contains the forces, directions and detail of all the supercells with displacements (equivalent to FORCE_SETS file in phonopy)
   :param force_constants: AiiDA ForceConstants object that contains the force constants
   :param nac_data: AiiDA ForceSetsData object that contains the Born effective charges and the dielectric tensor

- input parameters should be dictionary with the following entries ::

    parameters_dict = {'supercell': [[2, 0, 0],
                                     [0, 2, 0],
                                     [0, 0, 2]],
                       'primitive': [[1.0, 0.0, 0.0],
                                     [0.0, 1.0, 0.0],
                                     [0.0, 0.0, 1.0]],
                       'distance': 0.01,
                       'mesh': [40, 40, 40],
                       'symmetry_precision': 1e-5}
    ParameterData = DataFactory('parameter')
    parameters = ParameterData(dict=parameters_dict)

Either data_sets of force_constants should be used. If data_sets is used force constants will be calculated
and returned as a calculation output. If force_constants is used the calculation will be faster.

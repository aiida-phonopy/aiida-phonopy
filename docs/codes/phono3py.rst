Phono3py
========

This plugin is designed to calculate the thermal conductivity using phono3py.

.. function:: Phono3pyCalculation(structure, parameters, data_sets, nac_data)

   :param structure: StructureData object that contains the crystal unit cell information
   :param parameters: Dict object that contains the phonopy input parameters
   :param data_sets: ForceSetsData object that contains the forces, directions and detail of all the supercells with displacements (equivalent to FORCE_SETS file in phonopy)
   :param force_constants: ForceConstants object that contains the 2nd order force constants
   :param force_constants_3: ForceConstants object that contains the 3rd order force constants
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
    Dict = DataFactory('dict')
    parameters = Dict(dict=parameters_dict)

Either data_sets of force_constants/force_constants_3 should be defined. If data_sets is used force constants
will be calculated during the thermal conductivity calculation.

The outputs of this plugin are:

* **kappa**: ArrayData object that contains the results stored in kappa-mxxx-gx.hdf5 output file.

Phono3pyDist
============

This WorkChain is designed to do a phono3py calculation using a distributed workload scheme as described in
   phono3py manual (https://atztogo.github.io/phono3py/workload-distribution.html). The input nodes have the
   same name as the usual phono3py calculation. An additional node that defines the number of computers in
   which distribute the calculation is defined in this WorkChain.

.. function:: PhononPhono3py(structure, ph_settings, es_settings [, optimize=True, use_nac=False, pressure= 0.0, calculate_fc=False])

   :param structure: AiiDA StructureData object that contains the crystal unit cell information
   :param parameters: AiiDA ParametersData data object that contains the phonopy input parameters
   :param data_sets: AiiDA ForceSetsData object that contains the forces, directions and detail of all the supercells with displacements (equivalent to FORCE_SETS file in phonopy)
   :param force_constants: AiiDA ForceConstants object that contains the 2nd order force constants
   :param force_constants_3: AiiDA ForceConstants object that contains the 3rd order force constants
   :param nac_data: AiiDA ForceSetsData object that contains the Born effective charges and the dielectric tensor
   :param gp_chunks: AiiDA Int object that defines the number of computers in which the calculation will be distributed

The results outputs of this WorkChain are the following :

* **kappa**: ArrayData object that contains the results stored in kappa-mxxx-gx.hdf5 output file.


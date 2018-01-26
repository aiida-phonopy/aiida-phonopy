Phono3pyDist
============

This WorkChain is designed to do a phono3py calculation using a distributed workload scheme as described in
phono3py manual (https://atztogo.github.io/phono3py/workload-distribution.html). The input nodes have the
same name as the usual phono3py calculation. An additional node that defines the number of computers in
which distribute the calculation is defined in this WorkChain.

.. function:: Phono3pyDist(structure, parameters [, data_sets=None, nac_data=None, force_constants=None, force_constants_3=None, gp_chunks=10])

   :param structure: StructureData object that contains the crystal unit cell information
   :param parameters: ParametersData data object that contains the phonopy input parameters
   :param data_sets: ForceSetsData object that contains the forces, directions and detail of all the supercells with displacements (equivalent to FORCE_SETS file in phonopy)
   :param force_constants: ForceConstantsData object that contains the 2nd order force constants
   :param force_constants_3: ForceConstantsData object that contains the 3rd order force constants
   :param nac_data: (optional) NacData object that contains the Born effective charges and the dielectric tensor
   :param gp_chunks: (optional) Int object that defines the number of computers in which the calculation will be distributed

The results outputs of this WorkChain are the following :

* **kappa**: ArrayData object that contains the results stored in kappa-mxxx-gx.hdf5 output file.


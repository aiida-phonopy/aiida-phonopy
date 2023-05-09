:py:mod:`aiida_phonopy.calculations.functions.data_utils`
=========================================================

.. py:module:: aiida_phonopy.calculations.functions.data_utils

.. autoapi-nested-parse::

   Calcfunctions Utils for aiida-phonopy DataTypes.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.calculations.functions.data_utils.CalcfunctionMixin



Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.calculations.functions.data_utils.get_unitcell
   aiida_phonopy.calculations.functions.data_utils.get_primitive
   aiida_phonopy.calculations.functions.data_utils.get_supercell
   aiida_phonopy.calculations.functions.data_utils.get_supercells_with_displacements
   aiida_phonopy.calculations.functions.data_utils.get_displacements
   aiida_phonopy.calculations.functions.data_utils.generate_preprocess_data
   aiida_phonopy.calculations.functions.data_utils.get_preprocess_with_new_displacements
   aiida_phonopy.calculations.functions.data_utils.generate_phonopy_data



.. py:function:: get_unitcell(preprocess_data: aiida_phonopy.data.preprocess.PreProcessData) -> aiida.orm.StructureData

   Get the unitcell of a PreProcessData as a StructureData.


.. py:function:: get_primitive(preprocess_data: aiida_phonopy.data.preprocess.PreProcessData) -> aiida.orm.StructureData

   Get the primitive cell of a PreProcessData as a StructureData.


.. py:function:: get_supercell(preprocess_data: aiida_phonopy.data.preprocess.PreProcessData) -> aiida.orm.StructureData

   Get the supercell (pristine) of a PreProcessData as a StructureData.


.. py:function:: get_supercells_with_displacements(preprocess_data: aiida_phonopy.data.preprocess.PreProcessData) -> dict[aiida.orm.StructureData]

   Get the supercells with displacements of a PreProcessData as a StructureData.


.. py:function:: get_displacements(preprocess_data: aiida_phonopy.data.preprocess.PreProcessData) -> aiida.orm.ArrayData

   Get the displacements of a PreProcessData as an ArrayData with array name `displacements`.


.. py:function:: generate_preprocess_data(structure: aiida.orm.StructureData, displacement_generator: orm.Dict | None = None, supercell_matrix: orm.List | None = None, primitive_matrix: orm.List | None = None, symprec: orm.Float | None = None, is_symmetry: orm.Float | None = None, distinguish_kinds: orm.Bool | None = None)

   Return a complete stored PreProcessData node.

   :param structure: structure data node representing the unitcell
   :type structure: :class:`~aiida.orm.StructureData`
   :param displacement_generator: dictionary containing the info for generating the displacements
   :type displacement_generator: :class:`~aiida.orm.Dict`
   :param supercell_matrix: supercell matrix, defaults to diag(1,1,1)
   :type supercell_matrix: :class:`~aiida.orm.List`, Optional
   :param primitive_matrix: primitive matrix, defaults to "auto"
   :type primitive_matrix: :class:`~aiida.orm.List`, Optional
   :param symprec: symmetry precision on atoms, defaults to 1e-5
   :type symprec: :class:`~aiida.orm.Float`, Optional
   :param is_symmetry: if using space group symmetry, defaults to True
   :type is_symmetry: :class:`~aiida.orm.Bool`, Optional
   :param distinguish_kinds: if distinguish names of same specie by symmetry, defaults to True
   :type distinguish_kinds: :class:`~aiida.orm.Bool`, Optional

   :return: :class:`aiida_phonopy.data.preprocess.PreProcessData` node


.. py:function:: get_preprocess_with_new_displacements(preprocess_data: aiida_phonopy.data.preprocess.PreProcessData, displacement_generator: aiida.orm.Dict) -> aiida_phonopy.data.preprocess.PreProcessData

   Get a new PreProcessData from an old one from new displacement generator settings.


.. py:function:: generate_phonopy_data(preprocess_data: aiida_phonopy.data.preprocess.PreProcessData, nac_parameters: orm.ArrayData | None = None, forces_index: orm.Int | None = None, **forces_dict) -> aiida_phonopy.data.phonopy.PhonopyData

   Create a PhonopyData node from a PreProcess(Phonopy)Data node.

   `Forces` must be passed as **kwargs**, since we are calling a calcfunction with a variable
   number of supercells forces.

   :param nac_parameters: ArrayData containing 'dielectric' and 'born_charges' as arrays
       with their correct shape
   :param forces_index: Int if a TrajectoryData is given, in order to get the correct slice of the array.
   :param forces_dict: dictionary of supercells forces as ArrayData stored as `forces`, each Data
       labelled in the dictionary in the format `forces_{suffix}`.
       The prefix is common and the suffix corresponds to the suffix number of the supercell with
       displacement label given from the `get_supercells_with_displacements` method.

       For example:
           {'forces_1':ArrayData, 'forces_2':ArrayData}
           <==>
           {'supercell_1':StructureData, 'supercell_2':StructureData}
           and forces in each ArrayData stored as 'forces',
           i.e. ArrayData.get_array('forces') must not raise error

       .. note: if residual forces would be stored, label it with 0 as suffix.


.. py:class:: CalcfunctionMixin(data_node: PreProcessData | PhonopyData)

   Set of calcfunctions to be called from the aiida-phonopy DataTypes.

   .. py:method:: get_unitcell() -> aiida.orm.StructureData

      Get the unitcell as a StructureData through a calfunction.


   .. py:method:: get_primitive_cell() -> aiida.orm.StructureData

      Get the primitive cell as a StructureData through a calfunction.


   .. py:method:: get_supercell() -> aiida.orm.StructureData

      Get the supercell (pristine) as a StructureData through a calfunction.


   .. py:method:: get_supercells_with_displacements() -> dict

      Get the supercells with displacements as a StructureData through a calfunction.


   .. py:method:: get_displacements() -> aiida.orm.ArrayData

      Get the displacements as an ArrayData through a calfunction.


   .. py:method:: get_preprocess_with_new_displacements(displacement_generator: aiida.orm.Dict) -> aiida_phonopy.data.preprocess.PreProcessData

      Create a PreProcessData node from a PreProcess/PhonopyData with a new set of displacements.

      :param displacement_generator: a `storable` dictionary


   .. py:method:: generate_phonopy_data(nac_parameters: orm.ArrayData | None = None, forces_index: orm.Int | None = None, **forces_dict) -> aiida_phonopy.data.phonopy.PhonopyData

      Create a PhonopyData node from a PreProcess(Phonopy)Data node.

      `Forces` must be passed as **kwargs**, since we are calling a calcfunction with a variable
      number of supercells forces.

      :param nac_parameters: ArrayData containing 'dielectric' and 'born_charges' as arrays
          with their correct shape
      :param forces_index: Int if a TrajectoryData is given, in order to get the correct slice of the array.
      :param forces_dict: dictionary of supercells forces as ArrayData stored as `forces`, each Data
          labelled in the dictionary in the format `{prefix}_{suffix}`.
          The prefix is common and the suffix corresponds to the suffix number of the supercell with
          displacement label given from the `get_supercells_with_displacements` method

          For example:
              {'forces_1':ArrayData, 'forces_2':ArrayData} goes along with
              {'supercell_1':StructureData, 'supercell_2':StructureData}
              and forces in each ArrayData stored as 'forces',
              i.e. ArrayData.get_array('forces') must not raise error

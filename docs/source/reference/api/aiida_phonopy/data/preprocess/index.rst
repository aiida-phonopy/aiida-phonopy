:py:mod:`aiida_phonopy.data.preprocess`
=======================================

.. py:module:: aiida_phonopy.data.preprocess

.. autoapi-nested-parse::

   Module defining the class for managing the frozen phonon structure.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.data.preprocess.PreProcessData



Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.data.preprocess._serialize



.. py:class:: PreProcessData(structure: orm.StructureData | None = None, phonopy_atoms: PhonopyAtoms | None = None, supercell_matrix: list | None = None, primitive_matrix: list | None = None, symprec: float = 1e-05, is_symmetry: bool = True, distinguish_kinds: bool = True, **kwargs)

   Bases: :py:obj:`aiida_phonopy.data.raw.RawData`

   Class for pre-processing of frozen-phonon calculations.

   This class is designed for handling the pre-process information regarding
   a frozen phonon calculation using ``Phonopy``. These regard the unitcell structure,
   the supercell and primitive matrix, as well as other symmetry information.

   .. py:property:: displacement_dataset
      :type: dict | list | None

      Get the dispacement dataset in a readible format for phonopy.

      If not set, None is returned.


   .. py:property:: displacements
      :type: list | dict

      Get the displacements to apply to the supercell.

      :returns: array with displacements; can be type-I or type-II (see :func:`phonopy.Phonopy.displacements`)


   .. py:property:: calcfunctions

      Namespace to access the calcfunction utilities.


   .. py:method:: get_displacements() -> list | dict

      Get the displacements to apply to the supercell.


   .. py:method:: set_displacements(distance: float = 0.01, is_plusminus: str = 'auto', is_diagonal: bool = True, is_trigonal: bool = False, number_of_snapshots: int | None = None, random_seed: int | None = None, temperature: float | None = None, cutoff_frequency: float | None = None)

      Set displacements for frozen phonon calculation.

      Refer to :py:func:`phonopy.Phonopy.generate_displacements` for a complete description of the inputs.

      .. note: temperature different from 0K can be given only when forces are set up.
      Thus, in the PreProcessData value different from zero will raise an error.

      :raises ValueError: if the inputs are not compatible with phonopy standards


   .. py:method:: _set_displacements(value: list | dict)

      Put in the repository the displacement dataset in json format.


   .. py:method:: set_displacements_from_dataset(dataset: dict | list)

      Set displacements for frozen phonon calculation from a dataset.

      Useful if you want to set displacements from a previously random generated
      displacement dataset, or for setting dataset for self-consistent harmonic
      approximation.

      :param dataset: dictionary or array like (numpy or list), compatible with phonopy

      :raises ValueError: if the inputs are not compatible with phonopy standards


   .. py:method:: get_phonopy_instance(symmetrize_nac: bool = None, factor_nac: float | None = None, **kwargs)

      Return a :class:`~phonopy.Phonopy` object with the current values.

      :param symmetrize_nac: whether or not to symmetrize the nac parameters using point group symmetry;
          defaults to self.is_symmetry
      :type symmetrize_nac: bool
      :param factor_nac: factor for non-analytical corrections; defaults to Hartree*Bohr
      :type factor_nac: float
      :param kwargs: for internal use to set the primitive cell


   .. py:method:: get_supercells_with_displacements() -> dict[aiida.orm.StructureData]

      Get the supercells with displacements for frozen phonon calculation.

      .. note: this is not linking in the provenance the output structures.
          Use the `self.calcfunctions.get_supercells_with_displacements` instead

      :returns: dictionary with StructureData nodes, None if
          the displacement dataset has not been set


   .. py:method:: generate_displacement_dataset(distance: float = 0.01, is_plusminus: str = 'auto', is_diagonal: bool = True, is_trigonal: bool = False, number_of_snapshots: int | None = None, random_seed: int | None = None, temperature: float | None = None, cutoff_frequency: float | None = None)

      Return the displacement dataset for frozen phonon calculation.

      Refer to :py:func:`phonopy.Phonopy.generate_displacements` for a complete description of the inputs.

      :raises ValueError: if the inputs are not compatible with phonopy standards


   .. py:method:: generate_preprocess_data(structure: aiida.orm.StructureData, displacement_generator: dict | None = None, supercell_matrix: list | None = None, primitive_matrix: list | None = None, symprec: float | None = None, is_symmetry: bool | None = None, distinguish_kinds: bool | None = None)
      :staticmethod:

      Return a complete stored PreProcessData node.

      :param structure: structure data node representing the unitcell
      :type structure: :class:`~aiida.orm.StructureData`
      :param displacement_generator: dictionary containing the info for generating the displacements,
          defaults to phonopy default (see phonopy doc)
      :type displacement_generator: :class:`~aiida.orm.Dict`
      :param supercell_matrix: supercell matrix, defaults to diag(1,1,1)
      :type supercell_matrix: :class:`~aiida.orm.List`, Optional
      :param primitive_matrix: primitive matrix, defaults to "auto"
      :type primitive_matrix: list, Optional
      :param symprec: symmetry precision on atoms, defaults to 1e-5
      :type symprec: float, Optional
      :param is_symmetry: if using space group symmetry, defaults to True
      :type is_symmetry: bool, Optional
      :param distinguish_kinds: if distinguish names of same specie by symmetry, defaults to True
      :type distinguish_kinds: bool, Optional

      :return: :class:`~aiida_phonopy.data.PreProcessData` node



.. py:function:: _serialize(data: dict | list)

   Serialize the data for displacement dataset, in case it contains numpy.ndarray.

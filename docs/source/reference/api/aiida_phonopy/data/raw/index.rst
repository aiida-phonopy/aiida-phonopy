:py:mod:`aiida_phonopy.data.raw`
================================

.. py:module:: aiida_phonopy.data.raw

.. autoapi-nested-parse::

   Module defining the base class for other Data types.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.data.raw.RawData



Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.data.raw._get_valid_matrix



.. py:function:: _get_valid_matrix(matrix: list | np.ndarray) -> numpy.ndarray

   Get and validate the `supercell_matrix` and `primitive_matrix` inputs.

   :param matrix: (3,1) or (3,3) shape array
   :type matrix: class:`list`, class:`~aiida.orm.List`, :class:`numpy.ndarray`

   :return: a (3,3) :class:`numpy.ndarray`

   :raises TypeError: if it is not a valid array type and if the array does not contain only numbers
   :raises ValueError: if the array is not of the correct shape


.. py:class:: RawData(structure: orm.StructureData | None = None, phonopy_atoms: PhonopyAtoms | None = None, supercell_matrix: list | None = None, primitive_matrix: list | None = None, symprec: float = 1e-05, is_symmetry: bool = True, distinguish_kinds: bool = True, **kwargs)

   Bases: :py:obj:`aiida.orm.nodes.data.ArrayData`

   Base class containing the information for the other phonon related data types.

   .. py:property:: phonopy_version
      :type: str

      Get the Phonopy version used.


   .. py:property:: numbers
      :type: numpy.ndarray

      Get the array corresponding to the atomic number in periodic table of the sturcture.

      :return: (nat,) shaper array


   .. py:property:: masses
      :type: numpy.ndarray

      Get the array of the atomic masses in the cell.

      :return: (nat,) shape array


   .. py:property:: positions
      :type: numpy.ndarray

      Get the array of the atomic positions in the cell.

      :return: (nat, 3) shape array


   .. py:property:: cell
      :type: numpy.ndarray

      Get the lattice matrix of the structure.

      .. important: lattice vectors as rows of the matrix

      :return: (3,3) shape array


   .. py:property:: magnetic_moments
      :type: list[int]

      Get the `magnetic_moments` array of the atoms in the cell.


   .. py:property:: symbols
      :type: list[str]

      Get the chemical `symbols` array of the atoms in the cell.


   .. py:property:: pbc
      :type: tuple[int, int, int]

      Get the periodic boundary conditions.


   .. py:property:: names
      :type: list[str]

      Get the custom `names` array of the atoms in the cell.

      .. note: if no special names are specified, this will be equal to `symbols`.


   .. py:property:: supercell_matrix
      :type: list

      Get the `supercell_matrix`.

      :return: a (3,3) shape array


   .. py:property:: primitive_matrix
      :type: list

      Get the `primitive_matrix`.

      :return: a (3,3) shape array


   .. py:property:: symprec
      :type: float

      Get the tolerance for symmetry analysis.


   .. py:property:: is_symmetry
      :type: bool

      Get `is_symmetry` value.

      It refers to whether Phonopy will use symmetries to reduce the number of
      displacements for frozen phonons.


   .. py:property:: kinds_map
      :type: dict | None

      Get the map bewtween of the `numbers` and the `symbols` and `names`.


   .. py:property:: distinguish_kinds
      :type: bool

      Get whether or not kinds with same chemical symbol will be distinguished by symmetry.


   .. py:property:: dielectric
      :type: np.ndarray | None

      Get the high-frequency dielectric tensor in Cartesian coordinates.


   .. py:property:: born_charges
      :type: numpy.ndarray

      Get the effective Born charges tensors in Cartesian coordinates.

      .. note:
          The indices refers to:
              1. Atomic index.
              2. Polarization index.
              3. Atomic displacement index.

      :returns: numpy.ndarray, shape (num primitive cell atoms, 3, 3)


   .. py:method:: _set_phonopy_version()

      Set the installed Phonopy version.


   .. py:method:: _set_unitcell_attributes(phonopy_atoms: phonopy.structure.cells.PhonopyAtoms, pbc: tuple[bool, bool, bool])

      Set the attributes for full reproducibility of the `PhonopyAtoms` class.


   .. py:method:: _set_symbols_and_names()

      Set the `symbols` and `names`.


   .. py:method:: _set_supercell_matrix(value: list | np.ndarray)

      Set the Phonopy supercell matrix.

      :param value: (3,3) or (3,1) shape array
      :type value: :class:list, :class:`~aiida.orm.List`, :class:`numpy.ndarray`

      :raises ModificationNotAllowed: if object is already stored


   .. py:method:: _set_primitive_matrix(value: list | np.ndarray)

      Set the primitive matrix.

      :param value: (3,3) or (3,1) shape array, or 'auto'
      :type value: :class:list, :class:`~aiida.orm.List`, :class:`numpy.ndarray`, :class:str
      :raises ModificationNotAllowed: if object is already stored


   .. py:method:: _set_symprec(value: float)

      Set the symmetry tolerance.

      :param value: tolerance for symmetry analysis. Check that you get
      the right symmetry of your structure before starting any calculation. Default is 1e-05.
      :type value: float

      :raises ModificationNotAllowed: if object is already stored
      :raises TypeError: if the input is not of type `float`


   .. py:method:: _set_is_symmetry(value: bool)

      Set whether to use the symmetries.

      Use with care and if you know what your are doing.

      :param value: whether to use or not the symmetries. Deafault is True.
      :type value: bool

      :raises ModificationNotAllowed: if object is already stored
      :raises TypeError: if the input is not of type `bool`


   .. py:method:: _set_kinds_map(value: dict)

      Set the kind names map between the PhonopyAtoms unitcell and a reference structure.

      This is needed since PhonopyAtoms does not support kind names.
      This attribute allows to get proper `StructureData` supercells with displacements.

      :param value: tuple with two dictionaries (numbers_to_names, numbers_to_symbols)


   .. py:method:: _set_distinguish_kinds(value: bool)

      Set whether or not to distinguish kinds.


   .. py:method:: _get_phonopy_atoms_unitcell(distinguish_kinds: bool) -> phonopy.structure.cells.PhonopyAtoms

      Get the PhonopyAtoms object using the internal attributes.


   .. py:method:: get_phonopy_instance(symmetrize_nac: bool | None = None, factor_nac: float | None = None, **kwargs) -> phonopy.Phonopy

      Return a :class:`phonopy.Phonopy` object with the current values.

      :param symmetrize_nac: whether or not to symmetrize the nac parameters using point group symmetry;
          defaults to self.is_symmetry
      :type symmetrize_nac: bool
      :param factor_nac: factor for non-analytical corrections, defaults to Hartree*Bohr
      :type factor_nac: float
      :param kwargs: for internal use to set the primitive cell


   .. py:method:: get_unitcell() -> aiida.orm.StructureData

      Get the `unitcell` as StructureData (not stored).


   .. py:method:: get_primitive_cell() -> aiida.orm.StructureData

      Get the `primitive cell` as StructureData (not stored).


   .. py:method:: get_supercell() -> aiida.orm.StructureData

      Get the pristine `supercell` as StructureData (not stored).


   .. py:method:: get_cells_mappings() -> dict[dict[list]]

      Return a dictionary containing the mappings among unit-, super- and primitive cell.

      :return: dictionary with the following key:pair structure:
          * primitive: {p2p_map: list, p2s_map: list, s2p_map: list}
          * supercell: {u2u_map: list, u2s_map: list, s2u_map: list}


   .. py:method:: set_dielectric(dielectric: list | np.ndarray)

      Set the high-frequency dielectric tensor in Cartesian coordinates.

      .. note: it is assumed that the reference system is the same of the primitive cell.

      :param dielectric: (3, 3) array like

      :raises TypeError: if the format is not compatible or of the correct type
      :raises ValueError: if the format is not compatible or of the correct type


   .. py:method:: set_born_charges(born_charges: list | np.ndarray)

      Set the Born effective charge tensors in Cartesian coordinates.

      ..note:
          The indecis refers to:
              1. Atomic index.
              2. Polarization index.
              3. Atomic displacement index.

      :param born_charges: (number of atoms in the primitive cell, 3, 3) shape array like

      :raises:
          * TypeError: if the format is not compatible or of the correct type
          * ValueError: if the format is not compatible or of the correct type


   .. py:method:: has_nac_parameters() -> bool

      Return wheter or not the Data has non-analytical constants.


   .. py:method:: _if_can_modify()

      Check if the object is stored and raise an error if so. To use in every setter.

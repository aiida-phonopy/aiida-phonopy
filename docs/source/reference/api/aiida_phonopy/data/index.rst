:py:mod:`aiida_phonopy.data`
============================

.. py:module:: aiida_phonopy.data

.. autoapi-nested-parse::

   DataTypes for handling phonopy and frozen phonons calculations.



Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   force_constants/index.rst
   phonopy/index.rst
   preprocess/index.rst
   raw/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.data.RawData
   aiida_phonopy.data.ForceConstantsData
   aiida_phonopy.data.PreProcessData
   aiida_phonopy.data.PhonopyData
   aiida_phonopy.data.RawData
   aiida_phonopy.data.PreProcessData
   aiida_phonopy.data.RawData




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



.. py:class:: ForceConstantsData(structure: StructureData | None = None, phonopy_atoms: PhonopyAtoms | None = None, supercell_matrix: list | None = None, primitive_matrix: list | None = None, symprec: float = 1e-05, is_symmetry: bool = True, distinguish_kinds: bool = True, **kwargs)

   Bases: :py:obj:`aiida_phonopy.data.raw.RawData`

   Self-contained class for force constants data and non-analytical constants.

   It stores also the structure information (unitcell, supercell, ...), for a
   complete transferable data type.

   .. py:property:: force_constants
      :type: numpy.ndarray

      Get force force constants matrix.


   .. py:method:: get_phonopy_instance(**kwargs)

      Return a :class:`~phonopy.Phonopy` object with force and nac parameters (if set).

      :param kwargs: see :func:`aiida_phonopy.data.preprocess.PreProcessData.get_phonopy_instance`
          * symmetrize_nac: whether or not to symmetrize the nac parameters
              using point group symmetry; bool, defaults to self.is_symmetry
          * factor_nac: factor for non-analytical corrections;
              float, defaults to Hartree*Bohr


   .. py:method:: set_force_constants(force_constants: list | np.ndarray)

      Set force constants matrix.

      :param force_constants: array of force constants matrix in compact or full format
      :param type: (n_patom = atoms in primitive cell, n_satom =  atoms in supercell)
          * Compact format: (n_patom, n_satom, 3, 3)
          * Full format: (n_satom, n_satom, 3, 3)

      :raises:
          * TypeError: if the format is not of the correct type
          * ValueError: if the format is not compatible
          * RuntimeError: if the displacement dataset was not initialize in input



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



.. py:class:: PhonopyData(preprocess_data: aiida_phonopy.data.preprocess.PreProcessData, **kwargs)

   Bases: :py:obj:`aiida_phonopy.data.preprocess.PreProcessData`

   Class wrapping the :class:`phonopy.Phonopy` class.

   It represents the final Data node status of a frozen phonon calculaiton.
   It stores information regarding the pre-processing, the displacements and
   forces dataset, and the (eventual) non-analytical constants.

   .. note: direct calculation of properties from this class is still not implemented.
       Use :class:`~aiida_phonopy.calculations.phonopy.PhonopyCalculation` for
       post-processing this data node, keeping the provenance and having the
       data produced in the correct data type.

   .. py:property:: residual_forces
      :type: numpy.ndarray

      Get the residual forces calculated on the pristine (i.e. no displaced) supercell structure (if set).

      ..note: if you have specified the `forces_index` this will be used as well here.


   .. py:property:: forces
      :type: numpy.ndarray

      Get forces for each supercell with displacements in the dataset as a unique array.


   .. py:property:: forces_index
      :type: int

      Return the index of the forces to use.


   .. py:method:: set_displacements()

      Set displacements cannot be accessed from PhonopyData.


   .. py:method:: set_displacements_from_dataset()

      Set displacements cannot be accessed from PhonopyData.


   .. py:method:: get_phonopy_instance(subtract_residual_forces: bool | None = None, **kwargs) -> phonopy.Phonopy

      Return a :class:`~phonopy.Phonopy` object with forces and nac parameters (if set).

      :param subtract_residual_forces: whether or not subract residual forces (if set); defaults to False
      :type subtract_residual_forces: bool

      :param kwargs: see :func:`aiida_phonopy.data.preprocess.PreProcessData.get_phonopy_instance`
          * symmetrize_nac: whether or not to symmetrize the nac parameters
              using point group symmetry; bool, defaults to self.is_symmetry
          * factor_nac: factor for non-analytical corrections;
              float, defaults to Hartree*Bohr


   .. py:method:: set_residual_forces(forces: list | np.ndarray)

      Set the residual forces of the pristine supercell.

      :param forces: (atoms in supercell, 3) array shape

      :raises:
          * TypeError: if the format is not of the correct type
          * ValueError: if the format is not compatible


   .. py:method:: set_forces(sets_of_forces: list | np.ndarray | None = None, dict_of_forces: dict | None = None, forces_index: int | None = None)

      Set forces per each supercell with displacement in the dataset.

      :param sets_of_forces: a set of atomic forces in displaced supercells. The order of
          displaced supercells has to match with that in displacement dataset.
      :param type: (supercells with displacements, atoms in supercell, 3) array shape
      :param dict_of_forces: dictionary of forces, in numpy.ndarray to store for each displacement.
          They keys for the dictionary must be passed as `forces_{num}`, where `num` corresponds to
          the associated supercell in the dataset. `num` starts from 1.
      :param forces_index: an integer storing in the database the index for forces. The `dict_of_forces`
          may be specified from `TrajectoryData` to reduce the amount of data saved in the repository.
          For example: forces_1 = [[actual array]] ==> forces_index = 0

      :raises:
          * TypeError: if the format is not of the correct type
          * ValueError: if the format is not compatible
          * RuntimeError: if the displacement dataset was not initialize in input


   .. py:method:: set_forces_index(value: int)

      Set the `forces_index` attribute.

      This is used for tacking a particular array index of each single forces set.
      This is useful to not duplicate data of e.g. trajectories.



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

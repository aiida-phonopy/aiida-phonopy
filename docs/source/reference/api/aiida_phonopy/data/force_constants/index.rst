:py:mod:`aiida_phonopy.data.force_constants`
============================================

.. py:module:: aiida_phonopy.data.force_constants

.. autoapi-nested-parse::

   Module defining the class for force constants data.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.data.force_constants.ForceConstantsData




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

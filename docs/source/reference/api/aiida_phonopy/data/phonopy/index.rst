:py:mod:`aiida_phonopy.data.phonopy`
====================================

.. py:module:: aiida_phonopy.data.phonopy

.. autoapi-nested-parse::

   Module defining the class which wraps the :class:`phonopy.Phonopy` main class.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.data.phonopy.PhonopyData




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

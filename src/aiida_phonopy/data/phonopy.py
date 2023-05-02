# -*- coding: utf-8 -*-
"""
This module defines the class which wraps the :func:`phonopy.Phonopy` main class.
"""

import numpy as np

from .preprocess import PreProcessData


class PhonopyData(PreProcessData):  # pylint: disable=too-many-ancestors
    """
    This class wraps the :class:`phonopy.Phonopy` class. It represents the final Data node status
    of a frozen phonon calculaiton. It stores information regarding the pre-processing,
    the displacements and forces dataset, and the (eventual) non-analytical constants.

    .. note: direct calculation of properties from this class is still not implemented.
        Use :class:`~aiida_phonopy.calculations.phonopy.PhonopyCalculation` for
        post-processing this data node, keeping the provenance and having the
        data produced in the correct data type.
    """

    def __init__(self, preprocess_data: PreProcessData, **kwargs):

        if isinstance(preprocess_data, PreProcessData):
            kwargs['structure'] = preprocess_data.get_unitcell()
            kwargs['supercell_matrix'] = preprocess_data.supercell_matrix
            kwargs['primitive_matrix'] = preprocess_data.primitive_matrix
            kwargs['symprec'] = preprocess_data.symprec
            kwargs['is_symmetry'] = preprocess_data.is_symmetry
            kwargs['distinguish_kinds'] = preprocess_data.distinguish_kinds
            super().__init__(**kwargs)

            # super().set_displacements_from_dataset(preprocess_data.displacement_dataset)
            dataset = preprocess_data.displacement_dataset

            if dataset is not None:
                self.base.attributes.set('displacement_dataset', dataset)
            else:
                raise ValueError('cannot instantiate object without having displacement dataset set')

    def set_displacements(self):
        raise RuntimeError('`displacements` cannot be changed for this Data')

    def set_displacements_from_dataset(self):
        raise RuntimeError('`displacements` cannot be changed for this Data')

    def get_phonopy_instance(self, subtract_residual_forces=None, **kwargs):
        """Return a `phonopy.Phonopy` object with force and nac parameters (if set).

        :param subtract_residual_forces: whether or not subract residual forces (if set)
        :type subtract_residual_forces: bool, defaults to False
        :param kwargs: see :func:`aiida_phonopy.data.preprocess.PreProcessData.get_phonopy_instance`
            * symmetrize_nac: whether or not to symmetrize the nac parameters
                using point group symmetry; bool, defaults to self.is_symmetry
            * factor_nac: factor for non-analytical corrections;
                float, defaults to Hartree*Bohr
        """
        if not isinstance(subtract_residual_forces, bool) and subtract_residual_forces is not None:
            raise TypeError('`subtract_residual_forces` not of the right type')

        if subtract_residual_forces and self.residual_forces is None:
            raise ValueError('`subtract_residual_forces` is set to True, but `residual_forces` are None')

        ph_instance = super().get_phonopy_instance(**kwargs)

        if self.displacement_dataset is not None:
            ph_instance.dataset = self.displacement_dataset

        try:
            if self.forces is not None:
                the_forces = self.forces
                if subtract_residual_forces:
                    the_forces = the_forces - self.residual_forces
                ph_instance.forces = the_forces
        except AttributeError:
            pass

        return ph_instance

    @property
    def residual_forces(self):
        """Get the residual forces calculated on the pristine (i.e. no displaced) supercell structure (if set).

        ..note: if you have specified the `forces_index` this will be used as well here.
        """
        try:
            if self.forces_index is not None:
                the_forces = self.get_array('residual_forces')[self.forces_index]
            else:
                the_forces = self.get_array('residual_forces')
        except KeyError:
            the_forces = None
        return the_forces

    def set_residual_forces(self, forces):
        """Set the residual forces of the pristine supercell.

        :param forces: (atoms in supercell, 3) array shape

        :raises:
            * TypeError: if the format is not of the correct type
            * ValueError: if the format is not compatible
        """
        self._if_can_modify()

        if not isinstance(forces, (list, np.ndarray)):
            raise TypeError('the input is not of the correct type')

        natoms = len(self.get_supercell().sites)

        the_forces = np.array(forces)

        if self.forces_index is not None:
            if the_forces[self.forces_index].shape == (natoms, 3):
                self.set_array('residual_forces', the_forces)
            else:
                raise ValueError('the array is not of the correct shape')
        else:
            if the_forces.shape == (natoms, 3):
                self.set_array('residual_forces', the_forces)
            else:
                raise ValueError('the array is not of the correct shape. Check also `forces_index`')

    @property
    def forces(self):
        """Get forces for each supercell with displacements in the dataset as a unique array."""
        try:
            the_forces = self.get_array('forces')
        except (KeyError, AttributeError):
            try:
                nsupercells = len(self.displacements)
                the_forces = np.zeros((nsupercells)).tolist()
                for i in range(nsupercells):
                    if self.forces_index is not None:
                        the_forces[i] = self.get_array(f'forces_{i+1}')[self.forces_index]
                    else:
                        the_forces[i] = self.get_array(f'forces_{i+1}')
                the_forces = np.array(the_forces)
            except (KeyError, AttributeError):
                the_forces = None
        return the_forces

    def set_forces(self, sets_of_forces=None, dict_of_forces=None, forces_index=None):
        """Set forces per each supercell with displacement in the dataset.

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
        """
        self._if_can_modify()

        if self.displacement_dataset is None:
            raise RuntimeError('the displacement dataset has not been set yet')

        nsupercells = len(self.displacements)
        natoms = len(self.get_supercell().sites)

        if sets_of_forces is not None:
            the_forces = np.array(sets_of_forces)

            if the_forces.shape == (nsupercells, natoms, 3):
                self.set_array('forces', the_forces)
            else:
                raise ValueError('the array is not of the correct shape')

        if dict_of_forces is not None:
            # First, verify the number of keys is correct
            if len(list(dict_of_forces.keys())) != nsupercells:
                raise ValueError('the dictionary does not have the correct number of forces')
            # Second, verify the keys
            for key in dict_of_forces.keys():
                if key.split('_')[0] != 'forces':
                    raise ValueError(f'{key} is not correct. Expected `forces_num` as key')
            # Third, store
            for key, value in dict_of_forces.items():
                new_key = f"forces_{int(key.split('_')[-1])}"
                self.set_array(new_key, np.array(value))

        if forces_index is not None:
            if isinstance(forces_index, int):
                self.base.attributes.set('forces_index', forces_index)
            else:
                raise ValueError('index for forces must be an integer')

    @property
    def forces_index(self):
        """Return the index of the forces to use."""
        try:
            index = self.base.attributes.get('forces_index')
        except (KeyError, AttributeError):
            index = None
        return index

    def set_forces_index(self, value):
        """Set the `forces_index` attribute."""
        if isinstance(value, int):
            self.base.attributes.set('forces_index', value)
        else:
            raise ValueError('index for forces must be an integer')

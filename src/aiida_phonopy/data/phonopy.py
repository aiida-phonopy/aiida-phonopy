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

            super().set_displacements_from_dataset(preprocess_data.displacement_dataset)

        else:
            raise ValueError('`preprocess_data` is not of the correct type')

    def set_displacements(self):
        raise RuntimeError('`displacements` cannot be changed for this Data')

    def set_displacements_from_dataset(self):
        raise RuntimeError('`displacements` cannot be changed for this Data')

    def get_phonopy_instance(self, subtract_residual_forces=None, **kwargs):
        """Return a `phonopy.Phonopy` object with force and nac parameters (if set).

        :param subtract_residual_forces: whether or not subract residual forces (if set)
        :type subtract_residual_forces: bool, defaults to False
        :param kwargs: compatible keys with the super method in `RawData`
        """
        if not isinstance(subtract_residual_forces, bool) and subtract_residual_forces is not None:
            raise TypeError('`subtract_residual_forces` not of the right type')

        if subtract_residual_forces and self.residual_forces is None:
            raise ValueError('`subtract_residual_forces` is set to True, but `residual_forces` are None')

        ph_instance = super().get_phonopy_instance(**kwargs)

        if self.displacement_dataset is not None:
            ph_instance.dataset = self.displacement_dataset

        if self.forces is not None:
            the_forces = self.forces
            if subtract_residual_forces:
                the_forces = the_forces - self.residual_forces
            ph_instance.forces = the_forces

        return ph_instance

    @property
    def residual_forces(self):
        """Get the residual forces calculated on the pristine (i.e. no displaced) supercell structure (if set)."""
        try:
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

        if the_forces.shape == (natoms, 3):
            self.set_array('residual_forces', the_forces)
        else:
            raise ValueError('the array is not of the correct shape')

    @property
    def forces(self):
        """Get forces per each supercell with displacements in the dataset as a unique array."""
        try:
            the_forces = self.get_array('forces')
        except KeyError:
            the_forces = None
        return the_forces

    def set_forces(self, sets_of_forces):
        """Set forces per each supercell with displacement in the dataset.

        :param sets_of_forces:  a set of atomic forces in displaced supercells. The order of
            displaced supercells has to match with that in displacement dataset.
        :param type: (supercells with displacements, atoms in supercell, 3) array shape

        :raises:
            * TypeError: if the format is not of the correct type
            * ValueError: if the format is not compatible
            * RuntimeError: if the displacement dataset was not initialize in input
        """
        self._if_can_modify()

        if not isinstance(sets_of_forces, (list, np.ndarray)):
            raise TypeError('the input is not of the correct type')

        if self.displacement_dataset is None:
            raise RuntimeError('the displacement dataset has not been set yet')

        nsupercells = len(self.displacements)
        natoms = len(self.get_supercell().sites)

        the_forces = np.array(sets_of_forces)

        if the_forces.shape == (nsupercells, natoms, 3):
            self.set_array('forces', the_forces)
        else:
            raise ValueError('the array is not of the correct shape')

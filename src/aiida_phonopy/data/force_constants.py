# -*- coding: utf-8 -*-
"""
This module defines the class for force constants data.
"""
import numpy as np

from .raw import RawData


class ForceConstantsData(RawData):  # pylint: disable=too-many-ancestors
    """
    This class provides the force constants data, along with the non-analytical constants
    (optional) and the structure information (unitcell, supercell, ...).
    """

    def __init__(
        self,
        structure=None,
        phonopy_atoms=None,
        supercell_matrix=None,
        primitive_matrix=None,
        symprec=1e-05,
        is_symmetry=True,
        **kwargs
    ):
        kwargs['structure'] = structure
        kwargs['phonopy_atoms'] = phonopy_atoms
        kwargs['supercell_matrix'] = supercell_matrix
        kwargs['primitive_matrix'] = primitive_matrix
        kwargs['symprec'] = symprec
        kwargs['is_symmetry'] = is_symmetry

        super().__init__(**kwargs)

    def get_phonopy_instance(self, **kwargs):
        """Return a `phonopy.Phonopy` object with force and nac parameters (if set).

        :param kwargs: see :func:`aiida_phonopy.data.preprocess.PreProcessData.get_phonopy_instance`
            * symmetrize_nac: whether or not to symmetrize the nac parameters
                using point group symmetry; bool, defaults to self.is_symmetry
            * factor_nac: factor for non-analytical corrections;
                float, defaults to Hartree*Bohr
        """
        ph_instance = super().get_phonopy_instance(**kwargs)

        if self.force_constants is not None:
            ph_instance.set_force_constants(self.force_constants)

        return ph_instance

    @property
    def force_constants(self):
        """Get force force constants matrix."""
        try:
            the_forces = self.get_array('force_constants')
        except KeyError:
            the_forces = None
        return the_forces

    def set_force_constants(self, force_constants):
        """Set force constants matrix.

        :param force_constants: array of force constants matrix in compact or full format
        :param type: (n_patom = atoms in primitive cell, n_satom =  atoms in supercell)
            * Compact format: (n_patom, n_satom, 3, 3)
            * Full format: (n_satom, n_satom, 3, 3)

        :raises:
            * TypeError: if the format is not of the correct type
            * ValueError: if the format is not compatible
            * RuntimeError: if the displacement dataset was not initialize in input
        """
        self._if_can_modify()

        if not isinstance(force_constants, (list, np.ndarray)):
            raise TypeError('the input is not of the correct type')

        n_satoms = len(self.get_supercell().sites)
        n_patoms = len(self.get_primitive_cell().sites)

        fc = np.array(force_constants)

        if fc.shape in [(n_patoms, n_satoms, 3, 3), (n_satoms, n_satoms, 3, 3)]:
            self.set_array('force_constants', fc)
        else:
            raise ValueError('the array is not of the correct shape')

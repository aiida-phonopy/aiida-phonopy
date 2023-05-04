# -*- coding: utf-8 -*-
"""Module defining the class for force constants data."""
from __future__ import annotations

from aiida.orm import StructureData
import numpy as np
from phonopy.structure.cells import PhonopyAtoms

from .raw import RawData


class ForceConstantsData(RawData):  # pylint: disable=too-many-ancestors
    """Self-contained class for force constants data and non-analytical constants.

    It stores also the structure information (unitcell, supercell, ...), for a
    complete transferable data type.
    """

    def __init__(
        self,
        structure: StructureData | None = None,
        phonopy_atoms: PhonopyAtoms | None = None,
        supercell_matrix: list | None = None,
        primitive_matrix: list | None = None,
        symprec: float = 1e-05,
        is_symmetry: bool = True,
        distinguish_kinds: bool = True,
        **kwargs
    ):
        """Instantiate the class.

        The minimal input is to define either the `structure` or the `phonopy_atoms` input.
        They cannot be specified at the same time.

        :param structure: an :class:`~aiida.orm.StructureData` node
        :param phononpy_atoms: a :class:`~phonopy.structure.cells.PhonopyAtoms` instance
        :param supercell_matrix: a (3,3) shape array describing the supercell transformation
        :param primitive_matrix: a (3,3) shape array describing the primite transformation
        :param symprec: precision tollerance for symmetry analysis
        :param is_symmetry: whether using symmetries
        :distinguish_kinds: it stores a mapping between kinds and chemical symbols;
            by default Phonopy does not support kind,
            thus useful if in the input `structure` kinds are defined
        :paramm kwargs: for internal use
        """
        kwargs['structure'] = structure
        kwargs['phonopy_atoms'] = phonopy_atoms
        kwargs['supercell_matrix'] = supercell_matrix
        kwargs['primitive_matrix'] = primitive_matrix
        kwargs['symprec'] = symprec
        kwargs['is_symmetry'] = is_symmetry
        kwargs['distinguish_kinds'] = distinguish_kinds

        super().__init__(**kwargs)

    def get_phonopy_instance(self, **kwargs):
        """Return a :class:`~phonopy.Phonopy` object with force and nac parameters (if set).

        :param kwargs: see :func:`aiida_phonopy.data.preprocess.PreProcessData.get_phonopy_instance`
            * symmetrize_nac: whether or not to symmetrize the nac parameters
                using point group symmetry; bool, defaults to self.is_symmetry
            * factor_nac: factor for non-analytical corrections;
                float, defaults to Hartree*Bohr
        """
        ph_instance = super().get_phonopy_instance(**kwargs)

        if self.force_constants is not None:
            ph_instance.force_constants = self.force_constants

        return ph_instance

    @property
    def force_constants(self) -> np.ndarray:
        """Get force force constants matrix."""
        try:
            the_forces = self.get_array('force_constants')
        except KeyError:
            the_forces = None
        return the_forces

    def set_force_constants(self, force_constants: list | np.ndarray):
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

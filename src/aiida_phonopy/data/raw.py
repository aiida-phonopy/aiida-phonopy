# -*- coding: utf-8 -*-
"""Module defining the base class for other Data types."""
from __future__ import annotations

from copy import deepcopy

from aiida import orm
from aiida.orm.nodes.data import ArrayData
import numpy as np
from phonopy import Phonopy
from phonopy.structure.cells import PhonopyAtoms

from aiida_phonopy.calculations.functions.link_structures import (
    phonopy_atoms_from_structure,
    phonopy_atoms_to_structure,
)


def _get_valid_matrix(matrix: list | np.ndarray) -> np.ndarray:
    """Get and validate the `supercell_matrix` and `primitive_matrix` inputs.

    :param matrix: (3,1) or (3,3) shape array
    :type matrix: :class:`list`, :class:`~aiida.orm.List`, :class:`numpy.ndarray`
    :return: a (3,3) numpy.ndarray
    :raises:
        * TypeError: if it is not a valid array type and if the array does not contain only numbers
        * ValueError: if the array is not of the correct shape
    """
    if not isinstance(matrix, (list, orm.List, np.ndarray)):
        raise TypeError("only `list`, `aiida.orm.List` and `numpy.ndarray` (and 'auto' as string for primitive)")

    if isinstance(matrix, np.ndarray):
        matrix = matrix.tolist()

    if not len(matrix) == 3:
        raise ValueError('need exactly (3,1) or (3,3) shape.')

    for row in matrix:
        if isinstance(row, list):
            if not len(row) in [0, 3]:
                raise ValueError('matrix need to have (3,1) or (3,3) shape.')
            for element in row:
                if not isinstance(element, (int, float)):
                    raise TypeError(
                        f'type `{type(element)}` of {element} is not an accepted '
                        'type in matrix; only `int` and `float` are valid.'
                    )

    valid_matrix = np.array(matrix)
    valid_matrix = np.diag(valid_matrix) if valid_matrix.shape == (3,) else valid_matrix

    return valid_matrix


class RawData(ArrayData):  # pylint: disable=too-many-ancestors
    """Base class containing the information for the other phonon related data types."""

    def __init__(
        self,
        structure: orm.StructureData | None = None,
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
        if structure and phonopy_atoms:
            raise ValueError('cannot pass `structure`and `phonopy_atoms` at the same time')

        if structure:
            phonopy_atoms_unitcell, mapping = phonopy_atoms_from_structure(structure)
            pbc = structure.pbc
        elif phonopy_atoms:
            phonopy_atoms_unitcell, mapping = (phonopy_atoms, None)
            pbc = [True, True, True]  # to be changed when Phonopy will support PBC
        else:
            raise ValueError('at least one between `structure` and `phonopy_atoms` must be specified')

        args = {
            'unitcell': phonopy_atoms_unitcell,
            'supercell_matrix': supercell_matrix,
            'primitive_matrix': primitive_matrix,
            'symprec': symprec,
            'is_symmetry': is_symmetry,
        }

        try:
            Phonopy(**args)
        except Exception as err:
            raise ValueError('one or more inputs are not consistent with the `phonopy.Phonopy` format') from err

        super().__init__(**kwargs)

        # Inizializing the class attributes.
        self._set_phonopy_version()
        self._set_unitcell_attributes(phonopy_atoms=phonopy_atoms_unitcell, pbc=pbc)
        self._set_kinds_map(mapping)
        self._set_distinguish_kinds(distinguish_kinds)  # crucial before setting symbols/names!
        self._set_symbols_and_names()

        if not supercell_matrix is None:
            self._set_supercell_matrix(supercell_matrix)
        else:
            self._set_supercell_matrix(np.eye(3))

        self._set_symprec(symprec)
        self._set_is_symmetry(is_symmetry)

        # Important that is the last to be specified since, if "auto" is passed,
        # a phonopy instance is generated to get the primitive structure.
        if not primitive_matrix is None:
            self._set_primitive_matrix(primitive_matrix)
        else:
            self._set_primitive_matrix('auto')

    @property
    def phonopy_version(self) -> str:
        """Get the Phonopy version used."""
        return self.base.attributes.get('phonopy_version')

    def _set_phonopy_version(self):
        """Set the installed Phonopy version."""
        from phonopy.version import __version__ as the_phonopy_version

        self._if_can_modify()
        self.base.attributes.set('phonopy_version', the_phonopy_version)

    @property
    def numbers(self) -> np.ndarray:
        """Get the array corresponding to the atomic number in periodic table of the sturcture.

        :return: (nat,) shaper array
        """
        return self.base.attributes.get('numbers')

    @property
    def masses(self) -> np.ndarray:
        """Get the array of the atomic masses in the cell.

        :return: (nat,) shape array
        """
        return self.base.attributes.get('masses')

    @property
    def positions(self) -> np.ndarray:
        """Get the array of the atomic positions in the cell.

        :return: (nat, 3) shape array
        """
        return self.base.attributes.get('positions')

    @property
    def cell(self) -> np.ndarray:
        """Get the lattice matrix of the structure.

        .. important: lattice vectors as rows of the matrix

        :return: (3,3) shape array
        """
        return self.base.attributes.get('cell')

    @property
    def magnetic_moments(self) -> list[int]:
        """Get the `magnetic_moments` array of the atoms in the cell."""
        return self.base.attributes.get('magnetic_moments')

    @property
    def symbols(self) -> list[str]:
        """Get the chemical `symbols` array of the atoms in the cell."""
        return self.base.attributes.get('symbols')

    @property
    def pbc(self) -> tuple[int, int, int]:
        """Get the periodic boundary conditions."""
        return self.base.attributes.get('pbc')

    @property
    def names(self) -> list[str]:
        """Get the custom `names` array of the atoms in the cell.

        .. note: if no special names are specified, this will be equal to `symbols`.
        """
        return self.base.attributes.get('names')

    def _set_unitcell_attributes(self, phonopy_atoms: PhonopyAtoms, pbc: tuple[bool, bool, bool]):
        """Set the attributes for full reproducibility of the `PhonopyAtoms` class."""
        self._if_can_modify()
        self.base.attributes.set_many({
            'numbers': phonopy_atoms.numbers,
            'masses': phonopy_atoms.masses,
            'positions': phonopy_atoms.positions,
            'cell': phonopy_atoms.cell,
            'magnetic_moments': phonopy_atoms.magnetic_moments,
            'pbc': pbc,
        })

    def _set_symbols_and_names(self):
        """Set the `symbols` and `names`."""
        if self.kinds_map:
            numbers_to_names, numbers_to_symbols = self.kinds_map
            symbols = []
            names = []
            for number in self.numbers:
                try:
                    names.append(numbers_to_names[number])
                    symbols.append(numbers_to_symbols[number])
                except KeyError:
                    names.append(numbers_to_names[str(number)])
                    symbols.append(numbers_to_symbols[str(number)])
            self.base.attributes.set('symbols', symbols)
            self.base.attributes.set('names', names)
        else:
            phonopy_atoms = self._get_phonopy_atoms_unitcell(distinguish_kinds=True)
            self.base.attributes.set('symbols', deepcopy(phonopy_atoms.symbols))
            self.base.attributes.set('names', deepcopy(phonopy_atoms.symbols))

    @property
    def supercell_matrix(self) -> list:
        """Get the `supercell_matrix`.

        :return: a (3,3) shape array
        """
        return self.base.attributes.get('supercell_matrix')

    def _set_supercell_matrix(self, value: list | np.ndarray):
        """Set the Phonopy supercell matrix.

        :param value: (3,3) or (3,1) shape array
        :type value: :class:`list`, :class:`~aiida.orm.List`, :class:`numpy.ndarray`

        :raises ModificationNotAllowed: if object is already stored
        """
        self._if_can_modify()
        the_supercell_matrix = _get_valid_matrix(value)

        self.base.attributes.set('supercell_matrix', deepcopy(the_supercell_matrix))

    @property
    def primitive_matrix(self) -> list:
        """Get the `primitive_matrix`.

        :return: a (3,3) shape array
        """
        return self.base.attributes.get('primitive_matrix')

    def _set_primitive_matrix(self, value: list | np.ndarray):
        """Set the primitive matrix.

        :param value: (3,3) or (3,1) shape array, or 'auto"
        :type value: :class:`list`, :class:`~aiida.orm.List`, :class:`numpy.ndarray`, :class:`str`
        :raises ModificationNotAllowed: if object is already stored
        """
        self._if_can_modify()

        if isinstance(value, str):
            the_primitive_matrix = self.get_phonopy_instance(**{'primitive': value}).primitive_matrix
        else:
            the_primitive_matrix = _get_valid_matrix(value)

        self.base.attributes.set('primitive_matrix', deepcopy(the_primitive_matrix))

    @property
    def symprec(self) -> float:
        """Get the tolerance for symmetry analysis."""
        return self.base.attributes.get('symprec')

    def _set_symprec(self, value: float):
        """Set the symmetry tolerance.

        :param value: tolerance for symmetry analysis. Check that you get
        the right symmetry of your structure before starting any calculation. Default is 1e-05.
        :type value: float
        :raises:
            * ModificationNotAllowed: if object is already stored
            * TypeError: if the input is not of type `float`
        """
        self._if_can_modify()
        if not isinstance(value, float):
            raise ValueError('only `float` type is accepted.')

        self.base.attributes.set('symprec', deepcopy(value))

    @property
    def is_symmetry(self) -> bool:
        """Get `is_symmetry` value.

        It refers to whether Phonopy will use symmetries to reduce the number of
        displacements for frozen phonons.
        """
        return self.base.attributes.get('is_symmetry')

    def _set_is_symmetry(self, value: bool):
        """Set whether to use the symmetries.

        Use with care and if you know what your are doing.

        :param value: whether to use or not the symmetries. Deafault is True.
        :type value: bool

        :raises:
            * ModificationNotAllowed: if object is already stored
            * TypeError: if the input is not of type `bool`
        """
        self._if_can_modify()
        if not isinstance(value, bool):
            raise ValueError('only `bool` type is accepted.')

        self.base.attributes.set('is_symmetry', deepcopy(value))

    @property
    def kinds_map(self) -> dict | None:
        """Get the map bewtween of the `numbers` and the `symbols` and `names`."""
        try:
            the_map = self.base.attributes.get('kinds_map')
        except AttributeError:
            the_map = None
        return the_map

    def _set_kinds_map(self, value: dict):
        """Set the kind names map between the PhonopyAtoms unitcell and a reference structure.

        This is needed since PhonopyAtoms does not support kind names.
        This attribute allows to get proper `StructureData` supercells with displacements.

        :param value: tuple with two dictionaries (numbers_to_names, numbers_to_symbols)
        """
        self._if_can_modify()
        self.base.attributes.set('kinds_map', value)

    @property
    def distinguish_kinds(self) -> bool:
        """Get whether or not kinds with same chemical symbol will be distinguished by symmetry."""
        return self.base.attributes.get('distinguish_kinds')

    def _set_distinguish_kinds(self, value: bool):
        """Set whether or not to distinguish kinds."""
        self._if_can_modify()
        if not isinstance(value, bool):
            raise ValueError('only `bool` type is accepted.')

        self.base.attributes.set('distinguish_kinds', value)

    def _get_phonopy_atoms_unitcell(self, distinguish_kinds: bool) -> PhonopyAtoms:
        """Get the PhonopyAtoms object using the internal attributes."""
        kwargs = {
            'cell': self.cell,
            'positions': self.positions,
            'masses': self.masses,
            'magmoms': self.magnetic_moments
        }
        if distinguish_kinds:
            kwargs.update({'numbers': self.numbers})
        else:
            kwargs.update({'symbols': self.symbols})

        return PhonopyAtoms(**kwargs)

    def get_phonopy_instance(
        self, symmetrize_nac: bool | None = None, factor_nac: float | None = None, **kwargs
    ) -> Phonopy:
        """Return a :py:class:`phonopy.Phonopy` object with the current values.

        :param symmetrize_nac: whether or not to symmetrize the nac parameters using point group symmetry;
            defaults to self.is_symmetry
        :type symmetrize_nac: bool
        :param factor_nac: factor for non-analytical corrections; defaults to Hartree*Bohr
        :type factor_nac: float
        :param kwargs: for internal use to set the primitive cell
        """
        from phonopy.structure.symmetry import symmetrize_borns_and_epsilon
        from phonopy.units import Bohr, Hartree

        distinguish = kwargs.pop('distinguish', self.distinguish_kinds)
        unitcell = self._get_phonopy_atoms_unitcell(distinguish)

        primitive_matrix = kwargs.pop('primitive', None)
        if primitive_matrix is None:
            primitive_matrix = self.primitive_matrix

        ph_instance = Phonopy(
            unitcell=unitcell,
            supercell_matrix=self.supercell_matrix,
            primitive_matrix=primitive_matrix,
            symprec=self.symprec,
            is_symmetry=self.is_symmetry
        )

        # Non-analytical parameters
        if not isinstance(symmetrize_nac, bool) and symmetrize_nac is not None:
            raise TypeError('`symmetrize_nac` not of the right type')

        if symmetrize_nac is not None:
            do_symmetrize = symmetrize_nac
        else:
            do_symmetrize = self.is_symmetry

        if not isinstance(factor_nac, (float, int)) and factor_nac is not None:
            raise TypeError('input not of the right type')

        if factor_nac is not None:
            the_factor = factor_nac
        else:
            the_factor = Hartree * Bohr

        if self.has_nac_parameters():
            if do_symmetrize:
                nac = symmetrize_borns_and_epsilon(
                    borns=self.born_charges, epsilon=self.dielectric, ucell=ph_instance.primitive, symprec=self.symprec
                )
                the_born = nac[0]
                the_dielectric = nac[1]
            else:
                the_born = self.born_charges
                the_dielectric = self.dielectric

            ph_instance.nac_params = {'born': the_born, 'dielectric': the_dielectric, 'factor': the_factor}

        return ph_instance

    def get_unitcell(self) -> orm.StructureData:
        """Get the `unitcell` as StructureData (not stored)."""
        kwargs = {'distinguish': True}
        phonopy_instance = self.get_phonopy_instance(**kwargs).unitcell
        return phonopy_atoms_to_structure(phonopy_instance, self.kinds_map, self.pbc)

    def get_primitive_cell(self) -> orm.StructureData:
        """Get the `primitive cell` as StructureData (not stored)."""
        kwargs = {'distinguish': True}
        phonopy_instance = self.get_phonopy_instance(**kwargs).primitive
        return phonopy_atoms_to_structure(phonopy_instance, self.kinds_map, self.pbc)

    def get_supercell(self) -> orm.StructureData:
        """Get the pristine `supercell` as StructureData (not stored)."""
        kwargs = {'distinguish': True}
        phonopy_instance = self.get_phonopy_instance(**kwargs).supercell
        return phonopy_atoms_to_structure(phonopy_instance, self.kinds_map, self.pbc)

    def get_cells_mappings(self) -> dict[dict[list]]:
        """Return a dictionary containing the mappings among unit-, super- and primitive cell.

        :return: dictionary with the following key:pair structure:
            * primitive: {p2p_map: list, p2s_map: list, s2p_map: list}
            * supercell: {u2u_map: list, u2s_map: list, s2u_map: list}
        """
        ph = self.get_phonopy_instance()

        cells_maps = {
            'primitive': {
                'p2p_map': ph.primitive.p2p_map,
                'p2s_map': ph.primitive.p2s_map,
                's2p_map': ph.primitive.s2p_map,
            },
            'supercell': {
                'u2u_map': ph.supercell.u2u_map,
                'u2s_map': ph.supercell.u2s_map,
                's2u_map': ph.supercell.s2u_map,
            },
        }

        return cells_maps

    @property
    def dielectric(self) -> np.ndarray | None:
        """Get the high-frequency dielectric tensor in Cartesian coordinates."""
        try:
            value = self.base.attributes.get('dielectric')
            value = np.array(value)
        except (KeyError, AttributeError):
            value = None
        return value

    def set_dielectric(self, dielectric: list | np.ndarray):
        """Set the high-frequency dielectric tensor in Cartesian coordinates.

        .. note: it is assumed that the reference system is the same of the primitive cell.

        :param dielectric: (3, 3) array like

        :raises:
            * TypeError: if the format is not compatible or of the correct type
            * ValueError: if the format is not compatible or of the correct type
        """
        self._if_can_modify()

        if not isinstance(dielectric, (list, np.ndarray)):
            raise TypeError('the input is not of the correct type')

        the_dielectric = np.array(dielectric)

        if the_dielectric.shape == (3, 3):
            self.base.attributes.set('dielectric', the_dielectric.tolist())
        else:
            raise ValueError('the array is not of the correct shape')

    @property
    def born_charges(self) -> np.ndarray:
        """Get the effective Born charges tensors in Cartesian coordinates.

        ..note:
            The indecis refers to:
                1. Atomic index.
                2. Polarization index.
                3. Atomic displacement index.

        :returns: numpy.ndarray, shape (num primitive cell atoms, 3, 3)
        """
        try:
            value = self.get_array('born_charges')
        except (KeyError, AttributeError):
            value = None
        return value

    def set_born_charges(self, born_charges: list | np.ndarray):
        """Set the Born effective charge tensors in Cartesian coordinates.

        ..note:
            The indecis refers to:
                1. Atomic index.
                2. Polarization index.
                3. Atomic displacement index.

        :param born_charges: (number of atoms in the primitive cell, 3, 3) shape array like

        :raises:
            * TypeError: if the format is not compatible or of the correct type
            * ValueError: if the format is not compatible or of the correct type
        """
        self._if_can_modify()

        if not isinstance(born_charges, (list, np.ndarray)):
            raise TypeError('the input is not of the correct type')

        the_born_charges = np.array(born_charges)

        nprimitive_atoms = len(self.get_primitive_cell().sites)

        if the_born_charges.shape == (nprimitive_atoms, 3, 3):
            self.set_array('born_charges', the_born_charges)
        else:
            raise ValueError('the array is not of the correct shape')

    def has_nac_parameters(self) -> bool:
        """Return wheter or not the Data has non-analytical constants."""
        return (self.dielectric is not None and self.born_charges is not None)

    def _if_can_modify(self):
        """Check if the object is stored and raise an error if so. To use in every setter."""
        from aiida.common.exceptions import ModificationNotAllowed

        if self.is_stored:
            raise ModificationNotAllowed('The PhonopyData object cannot be modified, it has already been stored')

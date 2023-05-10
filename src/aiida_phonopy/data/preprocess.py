# -*- coding: utf-8 -*-
"""Module defining the class for managing the frozen phonon structure."""
from __future__ import annotations

import copy
import json

from aiida import orm
import numpy as np
from phonopy.structure.cells import PhonopyAtoms

from aiida_phonopy.calculations.functions.link_structures import phonopy_atoms_to_structure

from .raw import RawData


class PreProcessData(RawData):  # pylint: disable=too-many-ancestors
    """Class for pre-processing of frozen-phonon calculations.

    This class is designed for handling the pre-process information regarding
    a frozen phonon calculation using ``Phonopy``. These regard the unitcell structure,
    the supercell and primitive matrix, as well as other symmetry information.
    """

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
        kwargs['structure'] = structure
        kwargs['phonopy_atoms'] = phonopy_atoms
        kwargs['supercell_matrix'] = supercell_matrix
        kwargs['primitive_matrix'] = primitive_matrix
        kwargs['symprec'] = symprec
        kwargs['is_symmetry'] = is_symmetry
        kwargs['distinguish_kinds'] = distinguish_kinds

        super().__init__(**kwargs)

        if self.__class__.__name__ == 'PreProcessData':
            self.set_displacements()

    @property
    def displacement_dataset(self) -> dict | list | None:
        """Get the dispacement dataset in a readible format for phonopy.

        If not set, None is returned.
        """
        message = '`displacement_dataset` stored in the database will be deprecated in v2.0.0'
        try:
            import warnings
            the_dataset = self.base.attributes.get('displacement_dataset')
            warnings.warn(message, DeprecationWarning)
            return the_dataset
        except AttributeError:
            filename = 'displacement_dataset.json'

            if filename in self.base.repository.list_object_names():
                with self.base.repository.open(filename, mode='rb') as handle:
                    return json.load(handle)

            return None

    @property
    def displacements(self) -> list | dict:
        """Get the displacements to apply to the supercell.

        :returns: array with displacements; can be type-I or type-II (see :func:`phonopy.Phonopy.displacements`)
        """
        # For the time being, it is important to keep this implementation for the subclasses in the package,
        # otherwise a loop is generated (i.e. in the PhonopyData class, when setting the forces).
        ph = super().get_phonopy_instance()
        ph.dataset = self.displacement_dataset
        return ph.displacements

    def get_displacements(self) -> list | dict:
        """Get the displacements to apply to the supercell."""
        return self.displacements

    def set_displacements(
        self,
        distance: float = 0.01,
        is_plusminus: str = 'auto',
        is_diagonal: bool = True,
        is_trigonal: bool = False,
        number_of_snapshots: int | None = None,
        random_seed: int | None = None,
        temperature: float | None = None,
        cutoff_frequency: float | None = None,
    ):
        """Set displacements for frozen phonon calculation.

        Refer to :py:func:`phonopy.Phonopy.generate_displacements` for a complete description of the inputs.

        .. note: temperature different from 0K can be given only when forces are set up.
        Thus, in the PreProcessData value different from zero will raise an error.

        :raises ValueError: if the inputs are not compatible with phonopy standards
        """
        self._if_can_modify()
        ph = self.get_phonopy_instance()

        try:
            ph.generate_displacements(
                distance=distance,
                is_plusminus=is_plusminus,
                is_diagonal=is_diagonal,
                is_trigonal=is_trigonal,
                number_of_snapshots=number_of_snapshots,
                random_seed=random_seed,
                temperature=temperature,
                cutoff_frequency=cutoff_frequency,
            )
        except ValueError as err:
            raise ValueError('one or more input types are not accepted') from err

        self._set_displacements(copy.deepcopy(ph.dataset))

    def _set_displacements(self, value: list | dict):
        """Put in the repository the displacement dataset in json format."""
        try:
            serialized = json.dumps(value)
        except TypeError:
            serialized = json.dumps(_serialize(value))

        self.base.repository.put_object_from_bytes(serialized.encode('utf-8'), 'displacement_dataset.json')

    def set_displacements_from_dataset(self, dataset: dict | list):
        """Set displacements for frozen phonon calculation from a dataset.

        Useful if you want to set displacements from a previously random generated
        displacement dataset, or for setting dataset for self-consistent harmonic
        approximation.

        :param dataset: dictionary or array like (numpy or list), compatible with phonopy

        :raises ValueError: if the inputs are not compatible with phonopy standards
        """
        self._if_can_modify()
        ph = self.get_phonopy_instance()

        if isinstance(dataset, dict):
            try:
                ph.dataset = dataset
            except RuntimeError as err:
                raise ValueError('non compatible format') from err
        elif isinstance(dataset, (np.ndarray, list)):
            try:
                ph.displacements = dataset
            except RuntimeError as err:
                raise ValueError('non compatible format') from err
        else:
            raise ValueError('type not accepted')

        self._set_displacements(_serialize(copy.deepcopy(ph.dataset)))

    def get_phonopy_instance(self, symmetrize_nac: bool = None, factor_nac: float | None = None, **kwargs):
        """Return a :class:`~phonopy.Phonopy` object with the current values.

        :param symmetrize_nac: whether or not to symmetrize the nac parameters using point group symmetry;
            defaults to self.is_symmetry
        :type symmetrize_nac: bool
        :param factor_nac: factor for non-analytical corrections; defaults to Hartree*Bohr
        :type factor_nac: float
        :param kwargs: for internal use to set the primitive cell
        """
        ph = super().get_phonopy_instance(symmetrize_nac, factor_nac, **kwargs)
        if self.displacement_dataset is not None:
            ph.dataset = self.displacement_dataset
        return ph

    def get_supercells_with_displacements(self) -> dict[orm.StructureData]:
        """Get the supercells with displacements for frozen phonon calculation.

        .. note: this is not linking in the provenance the output structures.
            Use the `self.calcfunctions.get_supercells_with_displacements` instead

        :returns: dictionary with StructureData nodes, None if
            the displacement dataset has not been set
        """
        # Here we distinguish the kind, but only for mapping purposes.
        # This does not affect the number of displacements, since they
        # have been set with the correct flag of not distinguishing kinds.
        ph = self.get_phonopy_instance(**{'distinguish': True})

        supercells = ph.supercells_with_displacements
        mapping = self.kinds_map
        structures_dict = {}

        digits = len(str(len(supercells)))  # this will make the labels more nice to read

        for i, scell in enumerate(supercells):  # start from 1 - better choice for gathering forces
            label = f'supercell_{str(i + 1).zfill(digits)}'
            structures_dict.update({label: phonopy_atoms_to_structure(scell, mapping, self.pbc)})

        return structures_dict

    def generate_displacement_dataset(
        self,
        distance: float = 0.01,
        is_plusminus: str = 'auto',
        is_diagonal: bool = True,
        is_trigonal: bool = False,
        number_of_snapshots: int | None = None,
        random_seed: int | None = None,
        temperature: float | None = None,
        cutoff_frequency: float | None = None,
    ):
        """Return the displacement dataset for frozen phonon calculation.

        Refer to :py:func:`phonopy.Phonopy.generate_displacements` for a complete description of the inputs.

        :raises ValueError: if the inputs are not compatible with phonopy standards
        """
        ph = self.get_phonopy_instance()

        try:
            ph.generate_displacements(
                distance=distance,
                is_plusminus=is_plusminus,
                is_diagonal=is_diagonal,
                is_trigonal=is_trigonal,
                number_of_snapshots=number_of_snapshots,
                random_seed=random_seed,
                temperature=temperature,
                cutoff_frequency=cutoff_frequency,
            )
        except ValueError as err:
            raise ValueError('one or more inputs made the phonopy instance to raise an error') from err

        return copy.deepcopy(ph.dataset)

    @property
    def calcfunctions(self):
        """Namespace to access the calcfunction utilities."""
        from aiida_phonopy.calculations.functions.data_utils import CalcfunctionMixin

        return CalcfunctionMixin(data_node=self)

    @staticmethod
    def generate_preprocess_data(
        structure: orm.StructureData,
        displacement_generator: dict | None = None,
        supercell_matrix: list | None = None,
        primitive_matrix: list | None = None,
        symprec: float | None = None,
        is_symmetry: bool | None = None,
        distinguish_kinds: bool | None = None,
    ):
        """Return a complete stored PreProcessData node.

        :param structure: structure data node representing the unitcell
        :type structure: :class:`~aiida.orm.StructureData`
        :param displacement_generator: dictionary containing the info for generating the displacements,
            defaults to phonopy default (see phonopy doc)
        :type displacement_generator: orm.Dict
        :param supercell_matrix: supercell matrix, defaults to diag(1,1,1)
        :type supercell_matrix: :class:`~aiida.orm.List`, Optional
        :param primitive_matrix: primitive matrix, defaults to "auto"
        :type primitive_matrix: :class:`~aiida.orm.List`, Optional
        :param symprec: symmetry precision on atoms, defaults to 1e-5
        :type symprec: :class:`~aiida.orm.Float`, Optional
        :param is_symmetry: if using space group symmetry, defaults to True
        :type is_symmetry: :class:`~aiida.orm.Bool`, Optional
        :param distinguish_kinds: if distinguish names of same specie by symmetry, defaults to True
        :type distinguish_kinds: :class:`~aiida.orm.Bool`, Optional

        :return: PreProcessData node
        """
        from aiida_phonopy.calculations.functions.data_utils import generate_preprocess_data

        kwargs = {}

        kwargs['structure'] = structure

        if displacement_generator is not None:
            kwargs['displacement_generator'] = displacement_generator

        if supercell_matrix is not None:
            kwargs['supercell_matrix'] = supercell_matrix

        if primitive_matrix is not None:
            kwargs['primitive_matrix'] = primitive_matrix

        if symprec is not None:
            kwargs['symprec'] = symprec

        if is_symmetry is not None:
            kwargs['is_symmetry'] = is_symmetry

        if distinguish_kinds is not None:
            kwargs['distinguish_kinds'] = distinguish_kinds

        return generate_preprocess_data(**kwargs)


def _serialize(data: dict | list):
    """Serialize the data for displacement dataset, in case it contains numpy.ndarray."""
    serialized = copy.deepcopy(data)

    try:
        for key, value in data.items():
            serialized[key] = _serialize(value)
    except AttributeError:
        if isinstance(serialized, list):
            for i, value in enumerate(data):
                serialized[i] = _serialize(value)
        else:
            if isinstance(data, np.ndarray):
                return copy.deepcopy(data.tolist())

    return serialized

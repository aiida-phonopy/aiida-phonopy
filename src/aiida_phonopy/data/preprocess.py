# -*- coding: utf-8 -*-
"""
This module defines the class for managing the frozen phonon structure.
"""

import copy

from aiida import orm
import numpy as np

from aiida_phonopy.calculations.functions.link_structures import phonopy_atoms_to_structure

from .raw import RawData


class PreProcessData(RawData):  # pylint: disable=too-many-ancestors
    """
    This class is designed for handling the pre-process information regarding
    a frozen phonon calculation using `Phonopy`. These regard the unitcell structure,
    the supercell and primitive matrix, as well as other symmetry information.
    """

    def __init__(
        self,
        structure=None,
        phonopy_atoms=None,
        supercell_matrix=None,
        primitive_matrix=None,
        symprec=1e-05,
        is_symmetry=True,
        distinguish_kinds=True,
        **kwargs
    ):
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
    def displacement_dataset(self):
        """Get the dispacement dataset in a readible format for phonopy."""
        try:
            the_displacement_dataset = self.base.attributes.get('displacement_dataset')
        except AttributeError:
            the_displacement_dataset = None
        return the_displacement_dataset

    @property
    def displacements(self):
        """Get the displacements to apply to the supercell.

        :returns: array with displacements; can be type-I or type-II (see :func:`phonopy.Phonopy.displacements`)
        """
        try:
            ph = self.get_phonopy_instance()
            ph.dataset = self.displacement_dataset

            return ph.displacements
        except (AttributeError, TypeError):
            return None

    def get_displacements(self):
        """Get the displacements to apply to the supercell."""
        return self.displacements

    def set_displacements(
        self,
        distance=0.01,
        is_plusminus='auto',
        is_diagonal=True,
        is_trigonal=False,
        number_of_snapshots=None,
        random_seed=None,
        temperature=None,
        cutoff_frequency=None,
    ):
        """Set displacements for frozen phonon calculation.
        (see `phonopy.Phonopy.generate_displacements` for a complete description of the inputs)

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

        self.base.attributes.set('displacement_dataset', copy.deepcopy(ph.dataset))

    def set_displacements_from_dataset(self, dataset):
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

        self.base.attributes.set('displacement_dataset', copy.deepcopy(ph.dataset))

    def get_supercells_with_displacements(self):
        """Get the supercells with displacements for frozen phonon calculation.

        :returns: dictionary with StructureData nodes (note: not stored), None if
            the displacement dataset has not been set
        """
        ph = self.get_phonopy_instance(**{'distinguish': True})
        try:
            ph.dataset = self.displacement_dataset

            supercells = ph.supercells_with_displacements
            mapping = self.kinds_map
            structures_dict = {}

            # this will make the labels more nice to read - not really needed though
            digits = len(str(len(supercells)))

            for i, scell in enumerate(supercells):
                # start from 1 - better choice for gathering forces
                # MAYBE JUST NUMBERS???
                label = f'supercell_{str(i + 1).zfill(digits)}'
                structures_dict.update({label: phonopy_atoms_to_structure(scell, mapping)})

            return structures_dict
        except AttributeError:
            return None

    def generate_displacement_dataset(
        self,
        distance=0.01,
        is_plusminus='auto',
        is_diagonal=True,
        is_trigonal=False,
        number_of_snapshots=None,
        random_seed=None,
        temperature=None,
        cutoff_frequency=None,
    ):
        """Return the displacement dataset for frozen phonon calculation.
        (see `phonopy.Phonopy.generate_displacements` for a complete description of the inputs)

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
        displacement_generator=None,
        supercell_matrix=None,
        primitive_matrix=None,
        symprec=None,
        is_symmetry=None,
        distinguish_kinds=None,
    ):
        """Returns a complete stored PreProcessData node.

        :param structure: structure data node representing the unitcell
        :type structure: orm.StructureData
        :param displacement_generator: dictionary containing the info for generating the displacements,
            defaults to phonopy default (see phonopy doc)
        :type displacement_generator: orm.Dict
        :param supercell_matrix: supercell matrix, defaults to diag(1,1,1)
        :type supercell_matrix: orm.List, optional
        :param primitive_matrix: primitive matrix, defaults to "auto"
        :type primitive_matrix: orm.List or orm.List, optional
        :param symprec: symmetry precision on atoms, defaults to 1e-5
        :type symprec: orm.Float, optional
        :param is_symmetry: if using space group symmetry, defaults to True
        :type is_symmetry: orm.Bool, optional
        :param distinguish_kinds: if distinguish names of same specie by symmetry, defaults to True
        :type distinguish_kinds: orm.Bool, optional
        :return: PreProcessData node
        """
        from aiida_phonopy.calculations.functions.data_utils import generate_preprocess_data_from_distance

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

        return generate_preprocess_data_from_distance(**kwargs)

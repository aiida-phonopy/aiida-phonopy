# -*- coding: utf-8 -*-
"""Functions for phonopy pre-process."""

import numpy as np

from aiida.engine import calcfunction
from aiida import orm

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms


# a similar one is desirable for phono3py --> phono3py_preprocess
@calcfunction
def phonopy_preprocess(structure, symmetry_tolerance, displacements, supercell_matrix):
    """Set up the pre-process phonopy calculation returning the needed information for frozen phonon calculation.

    :param structure: StructureData of the unit cell
    :param symmetry_tolerance: symmetry tolerance used by phonopy
    :param displacements: kwargs for the ``generate_displacements`` method of the Phonopy class
    :param supercell_matrix: list or list of lists which describes the supercell matrix.
    :return: dict containing supercells with displacements, primitive
    structure and matrix, and the displacement dataset.
    """
    cell, mapping = phonopy_atoms_from_structure(structure, if_to_map(structure))
    ph = Phonopy(
        cell,
        supercell_matrix=get_supercell_matrix(supercell_matrix),
        primitive_matrix="auto",
        symprec=symmetry_tolerance.value,
    )

    displacements = displacements.get_dict()
    if "dataset" in displacements:
        ph.dataset = displacements["dataset"]
    else:
        kwargs = {}
        for item in displacements.items():
            kwargs.update({item})
        ph.generate_displacements(**kwargs)

    structures = get_phonon_structures(ph, mapping)

    displacement_dataset = orm.Dict(dict=ph.dataset)

    primitive_matrix = orm.List(list=ph.primitive_matrix.tolist())

    return {**structures, "primitive_matrix": primitive_matrix, "displacement_dataset": displacement_dataset}


# !!!
# SHALL WE MOVE THE BELOW TO, E.G., UTILS OR COMMON ???
# !!!
def get_supercell_matrix(supercell_matrix):
    """Return a phonopy readible supercell matrix."""
    if len(np.ravel(supercell_matrix)) == 3:
        smat = np.diag(supercell_matrix)
    else:
        smat = np.array(supercell_matrix)

    return smat.tolist()


def phonopy_atoms_to_structure(cell, mapping=None):
    """Return a StructureData from a PhonopyAtoms instance.

    :param cell: a PhonopyAtoms instance
    :param mapping: a number to kinds and symbols map, defaults to None
    """
    if mapping is None:
        symbols = cell.symbols
        positions = cell.positions
        structure = orm.StructureData(cell=cell.cell)
        for symbol, position in zip(symbols, positions):
            structure.append_atom(position=position, symbols=symbol)
        return structure
    else:
        numbers_to_kinds, numbers_to_symbols = mapping
        symbols = []
        names = []
        for number in cell.get_atomic_numbers():
            names.append(numbers_to_kinds[number])
            symbols.append(numbers_to_symbols[number])
        positions = cell.positions
        structure = orm.StructureData(cell=cell.cell)
        for position, symbol, name in zip(positions, symbols, names):
            structure.append_atom(position=position, symbols=symbol, name=name)
        return structure


def phonopy_atoms_from_structure(structure, to_map=False):
    """Return a tuple containg the PhonopyAtoms and the mapping from a StructureData.

    :param structure: StructureData instance
    :param to_map: to use when kind names are used to map atoms numbers to names, defaults to False
    :return: it returns a PhonopyAtoms istance and the mapping
    """
    if not to_map:  # when kind_names=symbols are used in the structure
        cell = PhonopyAtoms(
            symbols=[site.kind_name for site in structure.sites],
            positions=[site.position for site in structure.sites],
            cell=structure.cell,
        )
        return (cell, None)
    else:  # when using custom kind_names (!=symbols) in the structure
        sites = structure.sites
        symbols = []  # chemical symbols
        positions = []
        masses = []
        names = []  # actual kind names
        numbers = []

        kind_names = structure.get_kind_names()
        kinds_to_numbers = {kind_names[i]: i for i in range(len(kind_names))}

        for site in sites:
            kind = structure.get_kind(site.kind_name)
            symbols.append(kind.symbol)
            masses.append(kind.mass)
            names.append(kind.name)
            positions.append(site.position)
            numbers.append(kinds_to_numbers[kind.name])

        numbers_to_symbols = {numbers[i]: symbols[i] for i in range(len(sites))}
        # I build a PhonopyAtoms instance using the numbers instead of the symbols.
        # This will change the symbols which will be automatically stored in the instance.
        # We keep track of this 'issue' by returning `mapping` used
        # afterwards to remap the kind_names to the atoms.
        # ===> WARNING: attention in the post processing!
        cell = PhonopyAtoms(
            numbers=numbers,
            positions=positions,
            masses=masses,
            cell=structure.cell,
        )
        numbers_to_kinds = {i: kind_names[i] for i in range(len(kind_names))}
        return (cell, [numbers_to_kinds, numbers_to_symbols])


def get_phonon_structures(ph, mapping=None):
    """Returns all the supercells and the primitive cell of a Phonopy instance
    as StructuData nodes. If there is a mapping, then it returns also the
    `phonopy_structures`, which are the one readable from phonopy.
    """
    if mapping is None:
        structures_dict = {}
        digits = len(str(len(ph.supercells_with_displacements)))
        for i, scell in enumerate(ph.supercells_with_displacements):
            structure = phonopy_atoms_to_structure(scell)
            label = "supercell_%s" % str(i + 1).zfill(digits)
            structure.label = "%s %s" % (structure.get_formula(mode="hill_compact"), label)
            structures_dict.update({label: structure})

        supercell_structure = phonopy_atoms_to_structure(ph.supercell)
        supercell_structure.label = "%s %s" % (supercell_structure.get_formula(mode="hill_compact"), "supercell")
        structures_dict.update({"supercell": supercell_structure})

        primitive_structure = phonopy_atoms_to_structure(ph.primitive)
        primitive_structure.label = "%s %s" % (primitive_structure.get_formula(mode="hill_compact"), "primitive cell")
        structures_dict.update({"primitive": primitive_structure})
        structures_dict = {"cells": structures_dict}
    else:
        structures_dict = {"phonopy_cells": {}, "cells": {}}

        digits = len(str(len(ph.supercells_with_displacements)))
        for i, scell in enumerate(ph.supercells_with_displacements):
            label = "supercell_%s" % str(i + 1).zfill(digits)
            # Mapped cells for actual calculations
            structure = phonopy_atoms_to_structure(scell, mapping)
            structure.label = "%s %s" % (structure.get_formula(mode="hill_compact"), label)
            structures_dict["cells"].update({label: structure})
            # Phonopy `original` cells
            structure_phonopy = phonopy_atoms_to_structure(scell)
            structure_phonopy.label = "%s %s" % (structure_phonopy.get_formula(mode="hill_compact"), label)
            structures_dict["phonopy_cells"].update({label: structure_phonopy})

        supercell_structure = phonopy_atoms_to_structure(ph.supercell, mapping)
        supercell_structure.label = "%s %s" % (supercell_structure.get_formula(mode="hill_compact"), "supercell")
        structures_dict["cells"].update({"supercell": supercell_structure})

        supercell_structure_phonopy = phonopy_atoms_to_structure(ph.supercell)
        supercell_structure_phonopy.label = "%s %s" % (
            supercell_structure_phonopy.get_formula(mode="hill_compact"),
            "supercell",
        )
        structures_dict["phonopy_cells"].update({"supercell": supercell_structure_phonopy})

        primitive_structure = phonopy_atoms_to_structure(ph.primitive, mapping)
        primitive_structure.label = "%s %s" % (primitive_structure.get_formula(mode="hill_compact"), "primitive cell")
        structures_dict["cells"].update({"primitive": primitive_structure})

        primitive_structure_phonopy = phonopy_atoms_to_structure(ph.primitive, mapping)
        primitive_structure_phonopy.label = "%s %s" % (
            primitive_structure_phonopy.get_formula(mode="hill_compact"),
            "primitive cell",
        )
        structures_dict["phonopy_cells"].update({"primitive": primitive_structure_phonopy})

    return structures_dict


def if_to_map(structure):
    """Return a bool according whether kind names and symbols are the same.

    :param structure: StructureData instance
    :return: True if kind names different from symbols are used, False otherwise
    """
    check_kinds = [True for kind in structure.kinds if kind.name != kind.symbol]
    if check_kinds:
        return True
    else:
        return False

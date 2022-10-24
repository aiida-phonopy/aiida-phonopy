# -*- coding: utf-8 -*-
"""Functions for linking PhonopyAtoms and StructureData."""

__all__ = ('phonopy_atoms_to_structure', 'phonopy_atoms_from_structure', 'if_to_map')


def phonopy_atoms_to_structure(cell, mapping=None):
    """Return a StructureData from a PhonopyAtoms instance.

    :param cell: a PhonopyAtoms instance
    :param mapping: a number to kinds and symbols map, defaults to None
    """
    from aiida import orm

    if mapping:
        numbers_to_kinds, numbers_to_symbols = mapping
        symbols = []
        positions = []
        names = []
        for number in cell.numbers:
            try:
                names.append(numbers_to_kinds[number])
                symbols.append(numbers_to_symbols[number])
            except KeyError:
                names.append(numbers_to_kinds[str(number)])
                symbols.append(numbers_to_symbols[str(number)])
    else:
        names = cell.symbols
        symbols = cell.symbols

    positions = cell.positions
    masses = cell.masses

    structure = orm.StructureData(cell=cell.cell)
    for position, symbol, name, mass in zip(positions, symbols, names, masses):
        structure.append_atom(position=position, symbols=symbol, name=name, mass=mass)

    return structure


def phonopy_atoms_from_structure(structure):
    """Return a tuple containg the PhonopyAtoms and the mapping from a StructureData.

    :param structure: StructureData instance
    :return: it returns a PhonopyAtoms istance and the mapping
    """
    from phonopy.structure.cells import PhonopyAtoms

    to_map = if_to_map(structure)

    if not to_map:  # when kind_names=symbols are used in the structure
        cell = PhonopyAtoms(
            symbols=[site.kind_name for site in structure.sites],
            positions=[site.position for site in structure.sites],
            cell=structure.cell,
        )
        return (cell, None)
    # when using custom kind_names (!=symbols) in the structure
    sites = structure.sites
    symbols = []  # chemical symbols
    positions = []
    masses = []
    names = []  # actual kind names
    numbers = []

    kind_names = structure.get_kind_names()
    # The numbers start from 1, i.e. from the hydrogen element.
    kinds_to_numbers = {kind_names[i]: i + 1 for i in range(len(kind_names))}

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
    cell = PhonopyAtoms(
        numbers=numbers,
        positions=positions,
        masses=masses,
        cell=structure.cell,
    )

    # Also here we start from 1.
    numbers_to_kinds = {i + 1: kind_names[i] for i in range(len(kind_names))}

    return (cell, [numbers_to_kinds, numbers_to_symbols])


def if_to_map(structure):
    """Return a bool according whether kind names and symbols are the same.

    :param structure: StructureData instance
    :return: True if kind names different from symbols are used, False otherwise
    """
    check_kinds = [True for kind in structure.kinds if kind.name != kind.symbol]
    return bool(check_kinds)

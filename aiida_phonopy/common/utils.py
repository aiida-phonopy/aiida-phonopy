import numpy as np
from aiida.work import workfunction
from aiida.orm import DataFactory


def phonopy_bulk_from_structure(structure):
    from phonopy.api_phonopy import Atoms as PhonopyAtoms
    bulk = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)
    return bulk


def get_path_using_seekpath(structure, band_resolution=30):
    import seekpath

    phonopy_structure = phonopy_bulk_from_structure(structure)
    cell = phonopy_structure.get_cell()
    scaled_positions = phonopy_structure.get_scaled_positions()
    numbers = phonopy_structure.get_atomic_numbers()

    structure = (cell, scaled_positions, numbers)
    path_data = seekpath.get_path(structure)

    labels = path_data['point_coords']

    band_ranges = []
    for set in path_data['path']:
        band_ranges.append([labels[set[0]], labels[set[1]]])

    bands = []
    for q_start, q_end in band_ranges:
        band = []
        for i in range(band_resolution+1):
            band.append(
                np.array(q_start) + (np.array(q_end) - np.array(q_start))
                / band_resolution * i)
        bands.append(band)

    band_structure = DataFactory('phonopy.band_structure')(
        bands=bands,
        labels=path_data['path'],
        unitcell=phonopy_structure.get_cell())
    return band_structure


@workfunction
def get_force_constants_from_phonopy(structure, ph_settings, force_sets):
    """Calculate the force constants locally using phonopy

    :param structure:
    :param phonopy_input: ParameterData object that contains phonopy settings
    :param force_sets: ForceSetsData object that contains the atomic forces
        and displacements info (datasets dict in phonopy)
    :return: ForceConstantsData object containing the 2nd order force
        constants calculated with phonopy
    """
    from phonopy import Phonopy

    # Generate phonopy phonon object

    phonon = Phonopy(
        phonopy_bulk_from_structure(structure),
        ph_settings.get_dict()['supercell_matrix'],
        primitive_matrix=ph_settings.get_dict()['primitive_matrix'],
        symprec=ph_settings.get_dict()['symmetry_precision'])

    phonon.generate_displacements(distance=ph_settings.dict.distance)

    # Build data_sets from forces of supercells with displacments
    phonon.set_displacement_dataset(force_sets.get_force_sets())
    phonon.produce_force_constants()

    ForceConstantsData = DataFactory('phonopy.force_constants')
    force_constants = ForceConstantsData(data=phonon.get_force_constants())

    return {'force_constants': force_constants}


def get_nac_from_data(born_charges, epsilon, structure):
    """Worfunction to extract nac ArrayData object from calc"""

    nac_data = DataFactory('phonopy.nac')(
        structure=structure,
        born_charges=born_charges.get_array('born_charges'),
        epsilon=epsilon.get_array('epsilon'))

    return nac_data


def get_primitive(structure, ph_settings):
    from phonopy import Phonopy

    phonon = Phonopy(
        phonopy_bulk_from_structure(structure),
        supercell_matrix=ph_settings.get_dict()['supercell_matrix'],
        primitive_matrix=ph_settings.get_dict()['primitive_matrix'],
        symprec=ph_settings.get_dict()['symmetry_precision'])
    primitive_phonopy = phonon.get_primitive()

    primitive_cell = primitive_phonopy.get_cell()
    symbols = primitive_phonopy.get_chemical_symbols()
    positions = primitive_phonopy.get_positions()

    primitive_structure = DataFactory('structure')(cell=primitive_cell)
    for symbol, position in zip(symbols, positions):
        primitive_structure.append_atom(position=position, symbols=symbol)

    return {'primitive_structure': primitive_structure}


def phonopy_atoms_to_structure(cell):
    symbols = cell.get_chemical_symbols()
    positions = cell.get_positions()
    structure = DataFactory('structure')(cell=cell.get_cell())
    for symbol, position in zip(symbols, positions):
        structure.append_atom(position=position, symbols=symbol)
    return structure

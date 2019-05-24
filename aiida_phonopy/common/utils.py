import numpy as np
from aiida.engine import workfunction, calcfunction
from aiida.plugins import DataFactory
from aiida.orm import load_node, Int

@calcfunction
def get_path_using_seekpath(structure, band_resolution=Int(30)):
    import seekpath

    num_division = int(band_resolution)
    phonopy_structure = phonopy_atoms_from_structure(structure)
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
        for i in range(num_division + 1):
            band.append(
                np.array(q_start) + (np.array(q_end) - np.array(q_start))
                / num_division * i)
        bands.append(band)

    band_structure = DataFactory('phonopy.band_structure')()
    band_structure.set_bands(bands)
    band_structure.set_labels(path_data['path'])
    band_structure.set_unitcell(phonopy_structure.get_cell())

    return band_structure


@workfunction
def get_forces_from_uuid(uuid):
    n = load_node(str(uuid))
    return {'output_forces': n.out.output_forces}


@workfunction
def get_born_epsilon_from_uuid(uuid):
    n = load_node(str(uuid))
    return {'output_born_charges': n.out.output_born_charges,
            'output_dielectrics': n.out.output_dielectrics}


# @workfunction
@calcfunction
def get_force_sets(**forces_dict):
    forces = []
    # for i in range(forces_dict['num_supercells']):
    for i in range(len(forces_dict)):
        label = "forces_%03d" % (i + 1)
        if label in forces_dict:
            forces.append(forces_dict[label].get_array('final'))
    force_sets = DataFactory('array')()
    force_sets.set_array('force_sets', np.array(forces))
    force_sets.label = 'force_sets'
    return force_sets


@workfunction
def get_force_sets_data(datasets, **forces_dict):
    dataset_dict = datasets.get_dict()
    forces = []
    for i in range(len(dataset_dict['first_atoms'])):
        label = "supercell_%03d" % (i + 1)
        forces.append(forces_dict[label].get_array('final'))
    ForceSetsData = DataFactory('phonopy.force_sets')
    force_sets = ForceSetsData(datasets=dataset_dict)
    force_sets.set_forces(forces)
    return force_sets


@calcfunction
def get_nac_params(born_charges, epsilon):
    """Worfunction to extract nac ArrayData object from calc"""

    nac_params = DataFactory('array')()
    nac_params.set_array('born_charges',
                         born_charges.get_array('born_charges'))
    nac_params.set_array('epsilon', epsilon.get_array('epsilon'))
    nac_params.label = 'born_charges & epsilon'

    return nac_params


@calcfunction
def get_nac_data(born_charges, epsilon, structure):
    """Worfunction to extract NacData object from calc"""

    nac_data = DataFactory('phonopy.nac')(
        structure=structure,
        born_charges=born_charges.get_array('born_charges'),
        epsilon=epsilon.get_array('epsilon'))

    return nac_data


@calcfunction
def get_force_constants(structure, phonon_settings, force_sets, dataset):
    params = {}
    phonon = get_phonopy_instance(structure, phonon_settings, params)
    phonon.dataset = dataset.get_dict()
    phonon.forces = force_sets.get_array('force_sets')
    phonon.produce_force_constants()
    force_constants = DataFactory('array')()
    force_constants.set_array('force_constants', phonon.force_constants)
    force_constants.set_array('p2s_map', phonon.primitive.p2s_map)
    force_constants.label = 'force_constants'

    return force_constants


@calcfunction
def get_phonon(structure, phonon_settings, force_constants, band_paths,
               **params):
    phonon = get_phonopy_instance(structure, phonon_settings, params)
    phonon.force_constants = force_constants.get_array('force_constants')
    phonon_settings_dict = phonon_settings.get_dict()

    # Normalization factor primitive to unit cell
    normalization_factor = (phonon.unitcell.get_number_of_atoms()
                            / phonon.primitive.get_number_of_atoms())

    # DOS
    phonon.set_mesh(phonon_settings_dict['mesh'],
                    is_eigenvectors=True,
                    is_mesh_symmetry=False)
    phonon.set_total_DOS(tetrahedron_method=True)
    phonon.set_partial_DOS(tetrahedron_method=True)
    total_dos = phonon.get_total_DOS()
    partial_dos = phonon.get_partial_DOS()
    PhononDosData = DataFactory('phonopy.phonon_dos')
    dos = PhononDosData()
    dos.set_frequencies(total_dos[0])
    dos.set_dos(total_dos[1] * normalization_factor)
    dos.set_partial_dos(np.array(partial_dos[1]) * normalization_factor)
    dos.set_atom_labels(np.array(phonon.primitive.get_chemical_symbols()))

    # THERMAL PROPERTIES (per primtive cell)
    phonon.set_thermal_properties()
    t, free_energy, entropy, cv = phonon.get_thermal_properties()

    # Stores thermal properties (per unit cell) data in DB as a workflow
    # result
    thermal_properties = DataFactory('array')()
    thermal_properties.set_array('temperature', t)
    thermal_properties.set_array('free_energy',
                                 free_energy * normalization_factor)
    thermal_properties.set_array('entropy', entropy * normalization_factor)
    thermal_properties.set_array('heat_capacity',
                                 cv * normalization_factor)
    thermal_properties.label = 'Thermal properties'

    # BAND STRUCTURE
    bands = band_paths
    phonon.set_band_structure(bands.get_bands())
    BandStructureData = DataFactory('phonopy.band_structure')
    band_structure = BandStructureData()
    band_structure.set_bands(bands.get_bands())
    band_structure.set_labels(bands.get_labels())
    band_structure.set_unitcell(bands.get_unitcell())
    band_structure.set_band_structure_phonopy(phonon.get_band_structure())

    return {'dos': dos,
            'thermal_properties': thermal_properties,
            'band_structure': band_structure}


def get_phonopy_instance(structure, phonon_settings, params):
    from phonopy import Phonopy
    phonon_settings_dict = phonon_settings.get_dict()
    if 'primitive_matrix' in phonon_settings_dict:
        phonon = Phonopy(
            phonopy_atoms_from_structure(structure),
            phonon_settings_dict['supercell_matrix'],
            primitive_matrix=phonon_settings_dict['primitive_matrix'],
            symprec=phonon_settings_dict['symmetry_tolerance'])
    else:
        phonon = Phonopy(
            phonopy_atoms_from_structure(structure),
            phonon_settings_dict['supercell_matrix'],
            symprec=phonon_settings_dict['symmetry_tolerance'])
    if 'nac_params' in params:
        from phonopy.interface import get_default_physical_units
        units = get_default_physical_units('vasp')
        factor = units['nac_factor']
        nac_params = {'born': params['nac_params'].get_array('born_charges'),
                      'dielectric': params['nac_params'].get_array('epsilon'),
                      'factor': factor}
        phonon.set_nac_params(nac_params)

    return phonon


def get_primitive(structure, ph_settings):
    from phonopy import Phonopy

    phonon = Phonopy(
        phonopy_atoms_from_structure(structure),
        supercell_matrix=ph_settings.get_dict()['supercell_matrix'],
        primitive_matrix=ph_settings.get_dict()['primitive_matrix'],
        symprec=ph_settings.get_dict()['symmetry_tolerance'])
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


def phonopy_atoms_from_structure(structure):
    from phonopy.structure.atoms import PhonopyAtoms
    cell = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)
    return cell

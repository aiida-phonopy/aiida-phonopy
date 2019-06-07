import numpy as np
from aiida.engine import calcfunction
from aiida.plugins import DataFactory


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
def get_phonon(structure, phonon_settings, force_constants, **params):
    phonon_settings_dict = phonon_settings.get_dict()
    ph = get_phonopy_instance(structure, phonon_settings_dict, params)
    ph.force_constants = force_constants.get_array('force_constants')

    # DOS
    ph.set_mesh(phonon_settings_dict['mesh'],
                is_eigenvectors=True,
                is_mesh_symmetry=False)
    ph.run_total_dos()
    ph.run_projected_dos()

    total_dos = ph.get_total_dos_dict()
    dos = DataFactory('array.xy')()
    dos.set_x(total_dos['frequency_points'], 'Frequency', 'THz')
    dos.set_y(total_dos['total_dos'], 'Total DOS', '1/THz')

    projected_dos = ph.get_projected_dos_dict()
    pdos = DataFactory('array.xy')()
    pdos_list = [pd for pd in projected_dos['projected_dos']]
    pdos.set_x(projected_dos['frequency_points'], 'Frequency', 'THz')
    pdos.set_y(pdos_list,
               ['Projected DOS', ] * len(pdos_list),
               ['1/THz', ] * len(pdos_list))

    # Thermal properties
    ph.set_thermal_properties()
    t, free_energy, entropy, cv = ph.get_thermal_properties()
    thermal_properties = DataFactory('array')()
    thermal_properties.set_array('temperature', t)
    thermal_properties.set_array('free_energy', free_energy)
    thermal_properties.set_array('entropy', entropy)
    thermal_properties.set_array('heat_capacity', cv)
    thermal_properties.label = 'Thermal properties'

    # Band structure
    qpoints_list, frequencies_list, labels_list = get_bands_data(ph)
    bs = DataFactory('array.bands')()
    bs.set_kpoints(np.array(qpoints_list))
    bs.set_bands(np.array(frequencies_list), units='THz')
    if 'primitive_matrix' not in phonon_settings_dict:
        bs.labels = labels_list

    return {'dos': dos,
            'pdos': pdos,
            'thermal_properties': thermal_properties,
            'band_structure': bs}


def get_bands_data(ph):
    ph.auto_band_structure()
    labels = [x.replace('$', '').replace('\\', '').replace('mathrm{', '').replace('}', '').upper()
              for x in ph.band_structure.labels]
    frequencies = ph.band_structure.frequencies
    qpoints = ph.band_structure.qpoints
    path_connections = ph.band_structure.path_connections

    qpoints_list = list(qpoints[0])
    frequencies_list = list(frequencies[0])
    labels_list = [(0, labels[0]), ]
    label_index = 1

    for pc, qs, fs in zip(path_connections[:-1], qpoints[1:], frequencies[1:]):
        if labels[label_index] == 'Gamma':
            labels_list.append((len(qpoints_list) - 1, labels[label_index]))
            if label_index < len(labels):
                labels_list.append((len(qpoints_list), labels[label_index]))
            label_index += 1
            qpoints_list += list(qs)
            frequencies_list += list(fs)
        elif pc:
            labels_list.append((len(qpoints_list) - 1, labels[label_index]))
            label_index += 1
            qpoints_list += list(qs[1:])
            frequencies_list += list(fs[1:])
        else:
            labels_list.append((len(qpoints_list) - 1, labels[label_index]))
            label_index += 1
            if label_index < len(labels):
                labels_list.append((len(qpoints_list), labels[label_index]))
                label_index += 1
            qpoints_list += list(qs)
            frequencies_list += list(fs)
    labels_list.append((len(qpoints_list) - 1, labels[-1]))

    return qpoints_list, frequencies_list, labels_list


def get_phonopy_instance(structure, phonon_settings_dict, params):
    from phonopy import Phonopy
    phonon = Phonopy(
        phonopy_atoms_from_structure(structure),
        phonon_settings_dict['supercell_matrix'],
        primitive_matrix='auto',
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

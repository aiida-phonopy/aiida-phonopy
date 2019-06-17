import numpy as np
from aiida.engine import calcfunction
from aiida.plugins import DataFactory
from phonopy.structure.atoms import PhonopyAtoms


@calcfunction
def get_phonon_cells(supercell_array, primitive_array):
    positions = supercell_array.get_array('positions')
    numbers = supercell_array.get_array('numbers')
    lattice = supercell_array.get_array('cell')
    cells = {}
    for i, scaled_positions in enumerate(positions):
        cell = PhonopyAtoms(numbers=numbers,
                            scaled_positions=scaled_positions,
                            cell=lattice)
        print(cell)
        structure = phonopy_atoms_to_structure(cell)
        if i == 0:
            label = "supercell"
        else:
            label = "supercell_%03d" % i
        structure.label = "%s %s" % (
            structure.get_formula(mode='hill_compact'), label)
        cells[label] = structure

    prim_scaled_positions = primitive_array.get_array('positions')
    prim_numbers = primitive_array.get_array('numbers')
    prim_lattice = primitive_array.get_array('cell')
    cell = PhonopyAtoms(numbers=prim_numbers,
                        scaled_positions=prim_scaled_positions,
                        cell=prim_lattice)
    structure = phonopy_atoms_to_structure(cell)
    structure.label = "%s %s" % (
        structure.get_formula(mode='hill_compact'), 'primitive cell')
    cells['primitive'] = structure
    return cells


@calcfunction
def get_phonon_setting_info(phonon_setting_info):
    settings = DataFactory('dict')(dict=phonon_setting_info.get_dict())
    settings.label = 'phonon_setting_info'
    return settings


@calcfunction
def get_force_sets(**forces_dict):
    forces = []
    for i in range(len(forces_dict)):
        label = "forces_%03d" % (i + 1)
        if label in forces_dict:
            forces.append(forces_dict[label].get_array('final'))

    assert len(forces) == len(forces_dict)

    force_sets = DataFactory('array')()
    force_sets.set_array('force_sets', np.array(forces))
    force_sets.label = 'force_sets'
    return force_sets


@calcfunction
def get_nac_params(born_charges, epsilon, nac_structure, **params):
    """Obtain Born effective charges and dielectric constants in primitive cell

    When Born effective charges and dielectric constants are calculated within
    phonopy workchain, those values are calculated in the primitive cell.
    However using immigrant, the cell may not be primitive cell and can be
    unit cell. In this case, conversion of data is necessary. This conversion
    needs information of the structure where those values were calcualted and
    the target primitive cell structure.

    Two kargs parameters
    primitive : StructureData
    symmetry_tolerance : Float

    """
    from phonopy.structure.symmetry import symmetrize_borns_and_epsilon

    borns = born_charges.get_array('born_charges')
    eps = epsilon.get_array('epsilon')

    nac_cell = phonopy_atoms_from_structure(nac_structure)
    kargs = {}
    if 'symmetry_tolerance' in params:
        kargs['symprec'] = params['symmetry_tolerance'].value
    if 'primitive' in params:
        pcell = phonopy_atoms_from_structure(params['primitive'])
        kargs['primitive'] = pcell
    borns_, epsilon_ = symmetrize_borns_and_epsilon(
        borns, eps, nac_cell, **kargs)

    nac_params = DataFactory('array')()
    nac_params.set_array('born_charges', borns_)
    nac_params.set_array('epsilon', epsilon_)
    nac_params.label = 'born_charges & epsilon'

    return nac_params


@calcfunction
def get_force_constants(structure, phonon_settings, force_sets):
    params = {}
    phonon = get_phonopy_instance(structure, phonon_settings, params)
    phonon.dataset = phonon_settings['displacement_dataset']
    phonon.forces = force_sets.get_array('force_sets')
    print(phonon.dataset)
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
    mesh = phonon_settings_dict['mesh']

    # Mesh
    total_dos, pdos, thermal_properties = get_mesh_property_data(ph, mesh)

    # Band structure
    bs = get_bands_data(ph)

    return {'dos': total_dos,
            'pdos': pdos,
            'thermal_properties': thermal_properties,
            'band_structure': bs}


def get_mesh_property_data(ph, mesh):
    ph.set_mesh(mesh)
    ph.run_total_dos()

    dos = get_total_dos(ph.get_total_dos_dict())

    ph.run_thermal_properties()
    tprops = get_thermal_properties(ph.get_thermal_properties_dict())

    ph.set_mesh(mesh, is_eigenvectors=True, is_mesh_symmetry=False)
    ph.run_projected_dos()
    pdos = get_projected_dos(ph.get_projected_dos_dict())

    return dos, pdos, tprops


def get_total_dos(total_dos):
    dos = DataFactory('array.xy')()
    dos.set_x(total_dos['frequency_points'], 'Frequency', 'THz')
    dos.set_y(total_dos['total_dos'], 'Total DOS', '1/THz')
    dos.label = 'Total DOS'
    return dos


def get_projected_dos(projected_dos):
    pdos = DataFactory('array.xy')()
    pdos_list = [pd for pd in projected_dos['projected_dos']]
    pdos.set_x(projected_dos['frequency_points'], 'Frequency', 'THz')
    pdos.set_y(pdos_list,
               ['Projected DOS', ] * len(pdos_list),
               ['1/THz', ] * len(pdos_list))
    pdos.label = 'Projected DOS'
    return pdos


def get_thermal_properties(thermal_properties):
    tprops = DataFactory('array.xy')()
    tprops.set_x(thermal_properties['temperatures'], 'Temperature', 'K')
    tprops.set_y([thermal_properties['free_energy'],
                  thermal_properties['entropy'],
                  thermal_properties['heat_capacity']],
                 ['Helmholtz free energy', 'Entropy', 'Cv'],
                 ['kJ/mol', 'J/K/mol', 'J/K/mol'])
    tprops.label = 'Thermal properties'
    return tprops


def get_bands_data(ph):
    ph.auto_band_structure()
    labels = [x.replace('$', '').replace('\\', '').replace('mathrm{', '').replace('}', '').upper()
              for x in ph.band_structure.labels]
    frequencies = ph.band_structure.frequencies
    qpoints = ph.band_structure.qpoints
    path_connections = ph.band_structure.path_connections
    label = "%s (%d)" % (ph.symmetry.dataset['international'],
                         ph.symmetry.dataset['number'])

    return get_bands(qpoints, frequencies, labels, path_connections,
                     label=label)


def get_bands(qpoints, frequencies, labels, path_connections, label=None):
    qpoints_list = list(qpoints[0])
    frequencies_list = list(frequencies[0])
    labels_list = [(0, labels[0]), ]
    label_index = 1

    for pc, qs, fs in zip(path_connections[:-1], qpoints[1:], frequencies[1:]):
        if labels[label_index] == 'GAMMA':
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

    bs = DataFactory('array.bands')()
    bs.set_kpoints(np.array(qpoints_list))
    bs.set_bands(np.array(frequencies_list), units='THz')
    bs.labels = labels_list
    if label is not None:
        bs.label = label

    return bs


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
    cell = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)
    return cell

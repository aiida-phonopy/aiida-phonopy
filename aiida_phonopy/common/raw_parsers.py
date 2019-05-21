import numpy as np
from aiida.plugins import DataFactory
from aiida_phonopy.common.utils import phonopy_atoms_from_structure


ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')

ForceConstantsData = DataFactory('phonopy.force_constants')
PhononDosData = DataFactory('phonopy.phonon_dos')
BandStructureData = DataFactory('phonopy.band_structure')


# Parse files to generate AIIDA OBJECTS

def parse_FORCE_CONSTANTS(filename):
    from phonopy.file_IO import parse_FORCE_CONSTANTS as parse_FC

    force_constants = parse_FC(filename=filename)
    fc_array = ArrayData()
    fc_array.set_array('force_constants', force_constants)
    fc_array.label = 'force_constants'
    return fc_array


def parse_partial_DOS(filename, structure, parameters):
    partial_dos = np.loadtxt(filename)

    from phonopy.structure.atoms import PhonopyAtoms
    from phonopy import Phonopy

    bulk = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)
    params_dict = parameters.get_dict()
    if 'primitive_matrix' in params_dict:
        phonon = Phonopy(
            bulk,
            supercell_matrix=params_dict['supercell_matrix'],
            primitive_matrix=params_dict['primitive_matrix'],
            symprec=params_dict['symmetry_tolerance'])
    else:
        phonon = Phonopy(
            bulk,
            supercell_matrix=params_dict['supercell_matrix'],
            symprec=params_dict['symmetry_tolerance'])

    dos = PhononDosData(
        frequencies=partial_dos.T[0],
        dos=np.sum(partial_dos[:, 1:], axis=1),
        partial_dos=partial_dos[:, 1:].T,
        atom_labels=phonon.get_primitive().get_chemical_symbols())

    return dos


def parse_thermal_properties(filename):
    import yaml
    temperature = []
    free_energy = []
    entropy = []
    cv = []

    with open(filename, 'r') as stream:
        thermal_properties = dict(yaml.load(stream))
        for tp in thermal_properties['thermal_properties']:
            temperature.append(tp['temperature'])
            entropy.append(tp['entropy'])
            free_energy.append(tp['free_energy'])
            cv.append(tp['heat_capacity'])

    tp_object = ArrayData()
    tp_object.set_array('temperature', np.array(temperature))
    tp_object.set_array('free_energy', np.array(free_energy))
    tp_object.set_array('entropy', np.array(entropy))
    tp_object.set_array('heat_capacity', np.array(cv))

    return tp_object


def parse_band_structure(filename, input_bands):
    import yaml

    frequencies = []
    with open(filename, 'r') as stream:
        bands = dict(yaml.load(stream))

    for k in bands['phonon']:
        frequencies.append([b['frequency'] for b in k['band']])
    frequencies = np.array(frequencies)

    nb = input_bands.get_number_of_bands()
    frequencies = frequencies.reshape((nb, -1, frequencies.shape[1]))

    band_structure = BandStructureData(frequencies=frequencies,
                                       bands=input_bands.get_bands(),
                                       labels=input_bands.get_labels(),
                                       unitcell=input_bands.get_unitcell())

    return band_structure


def parse_kappa(filename):
    import numpy as np
    import h5py
    from aiida.orm.nodes.data.array import ArrayData

    kappa = ArrayData()

    f = h5py.File(filename, 'r')
    for item in f:
        array = np.array(f[item])
        kappa.set_array(item, array)
    f.close()

    return kappa


# Generate text strings for files from AIIDA OBJECTS
def get_BORN_txt(nac_data, parameter_data, structure):
    from phonopy.file_IO import get_BORN_lines

    parameters = parameter_data.get_dict()
    born_charges = nac_data.get_array('born_charges')
    epsilon = nac_data.get_array('epsilon')
    unitcell = phonopy_atoms_from_structure(structure)
    params = {'supercell_matrix': parameters['supercell_matrix']}
    if 'symmetry_tolerance' in parameters:
        params['symprec'] = parameters['symmetry_tolerance']
    if 'primitive_matrix' in parameters:
        params['primitive_matrix'] = parameters['primitive_matrix']
    lines = get_BORN_lines(unitcell, born_charges, epsilon, **params)

    return "\n".join(lines)


def get_FORCE_SETS_txt(force_sets, displacements_dataset):
    from phonopy.file_IO import get_FORCE_SETS_lines

    forces = force_sets.get_array('force_sets')
    dataset = displacements_dataset.get_dict()
    lines = get_FORCE_SETS_lines(dataset, forces=forces)

    return "\n".join(lines)


def get_FORCE_CONSTANTS_txt(force_constants_object):
    from phonopy.file_IO import get_FORCE_CONSTANTS_lines

    force_constants = force_constants_object.get_array('force_constants')
    p2s_map = force_constants_object.get_array('p2s_map')
    lines = get_FORCE_CONSTANTS_lines(force_constants, p2s_map=p2s_map)

    return "\n".join(lines)


def get_poscar_txt(structure, use_direct=True):
    types = [site.kind_name for site in structure.sites]
    atom_type_unique = np.unique(types, return_index=True)
    sort_index = np.argsort(atom_type_unique[1])
    elements = np.array(atom_type_unique[0])[sort_index]
    elements_count = np.diff(np.append(
        np.array(atom_type_unique[1])[sort_index], [len(types)]))

    poscar_txt = '# VASP POSCAR generated using aiida workflow '
    poscar_txt += '\n1.0\n'
    cell = structure.cell
    for row in cell:
        poscar_txt += '{0: 22.16f} {1: 22.16f} {2: 22.16f}\n'.format(*row)
    poscar_txt += ' '.join([str(e) for e in elements]) + '\n'
    poscar_txt += ' '.join([str(e) for e in elements_count]) + '\n'
    if use_direct:
        poscar_txt += 'Direct\n'
    else:
        poscar_txt += 'Cartesian\n'
    for site in structure.sites:
        if use_direct:
            coordinates = np.dot(site.position, np.linalg.inv(cell))
        else:
            coordinates = site.position
        poscar_txt += '{0: 22.16f} {1: 22.16f} {2: 22.16f}\n'.format(*coordinates)

    return poscar_txt


def get_phonopy_conf_file_txt(parameters_object, bands=None):
    parameters = parameters_object.get_dict()
    supercell_matrix = np.array(parameters['supercell_matrix']).ravel()
    vals = " {}" * len(supercell_matrix)
    input_file = ('DIM =' + vals + '\n').format(*supercell_matrix)
    if 'primitive_matrix' in parameters:
        input_file += 'PRIMITIVE_AXIS = {} {} {}  {} {} {}  {} {} {}\n'.format(
            *np.array(parameters['primitive_matrix']).ravel())
    input_file += 'MESH = {} {} {}\n'.format(*parameters['mesh'])

    if bands is not None:
        input_file += 'BAND = '
        for i, band in enumerate(bands.get_band_ranges()):
            input_file += ' '.join(np.array(band.flatten(), dtype=str))
            if i < bands.get_number_of_bands() - 1:
                input_file += ','

        input_file += '\n'
        input_file += 'BAND_POINTS = {}\n'.format(bands.get_number_of_points())

    return input_file


def get_disp_fc3_txt(structure, parameters_data, force_sets):
    from phono3py.file_IO import write_cell_yaml
    from phonopy.structure.cells import get_supercell
    from aiida_phonopy.common.utils import phonopy_atoms_from_structure

    import StringIO

    dataset = force_sets.get_datasets3()
    supercell = get_supercell(phonopy_atoms_from_structure(structure),
                              parameters_data.dict.supercell,
                              symprec=parameters_data.dict.symmetry_tolerance)

    w = StringIO.StringIO()
    w.write("natom: %d\n" % dataset['natom'])

    num_first = len(dataset['first_atoms'])
    w.write("num_first_displacements: %d\n" % num_first)
    if 'cutoff_distance' in dataset:
        w.write("cutoff_distance: %f\n" % dataset['cutoff_distance'])

    num_second = 0
    num_disp_files = 0
    for d1 in dataset['first_atoms']:
        num_disp_files += 1
        num_second += len(d1['second_atoms'])
        for d2 in d1['second_atoms']:
            if 'included' in d2:
                if d2['included']:
                    num_disp_files += 1
            else:
                num_disp_files += 1

    w.write("num_second_displacements: %d\n" % num_second)
    w.write("num_displacements_created: %d\n" % num_disp_files)

    w.write("first_atoms:\n")
    count1 = 1
    count2 = num_first + 1
    for disp1 in dataset['first_atoms']:
        disp_cart1 = disp1['displacement']
        w.write("- number: %5d\n" % (disp1['number'] + 1))
        w.write("  displacement:\n")
        w.write("    [%20.16f,%20.16f,%20.16f ] # %05d\n" %
                (disp_cart1[0], disp_cart1[1], disp_cart1[2], count1))
        w.write("  second_atoms:\n")
        count1 += 1

        included = None
        atom2 = -1
        for disp2 in disp1['second_atoms']:
            if atom2 != disp2['number']:
                atom2 = disp2['number']
                if 'included' in disp2:
                    included = disp2['included']
                pair_distance = disp2['pair_distance']
                w.write("  - number: %5d\n" % (atom2 + 1))
                w.write("    distance: %f\n" % pair_distance)
                if included is not None:
                    if included:
                        w.write("    included: %s\n" % "true")
                    else:
                        w.write("    included: %s\n" % "false")
                w.write("    displacements:\n")

            disp_cart2 = disp2['displacement']
            w.write("    - [%20.16f,%20.16f,%20.16f ] # %05d\n" %
                    (disp_cart2[0], disp_cart2[1], disp_cart2[2], count2))
            count2 += 1

    write_cell_yaml(w, supercell)
    w.seek(0)
    lines = w.read()
    w.close()
    return lines


def get_forces_txt(force_sets):
    import StringIO

    w = StringIO.StringIO()

    from phono3py.file_IO import write_FORCES_FC3

    write_FORCES_FC3(force_sets.get_datasets3(),
                     force_sets.get_forces3(),
                     fp=w)
    w.seek(0)
    lines = w.read()
    w.close()
    return lines


def write_fc2_to_hdf5_file(force_constants, filename):
    from phono3py.file_IO import write_fc2_to_hdf5
    write_fc2_to_hdf5(force_constants.get_data(), filename)


def write_fc3_to_hdf5_file(force_constants, filename):
    from phono3py.file_IO import write_fc3_to_hdf5
    write_fc3_to_hdf5(force_constants.get_data(), filename)


def write_kappa_to_hdf5_file(gp, filename='kappa'):
    import h5py
    with h5py.File(filename, 'w') as w:
        for key in gp:
            w.create_dataset(key, data=gp[key])

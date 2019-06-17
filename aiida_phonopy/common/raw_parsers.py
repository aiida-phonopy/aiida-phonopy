import numpy as np
from aiida.plugins import DataFactory
from aiida_phonopy.common.utils import (phonopy_atoms_from_structure,
                                        get_total_dos, get_projected_dos,
                                        get_thermal_properties, get_bands)


ArrayData = DataFactory('array')
XyData = DataFactory('array.xy')
BandsData = DataFactory('array.bands')


# Parse files to generate AIIDA OBJECTS

def parse_FORCE_CONSTANTS(filename):
    from phonopy.file_IO import parse_FORCE_CONSTANTS as parse_FC

    force_constants = parse_FC(filename=filename)
    fc_array = ArrayData()
    fc_array.set_array('force_constants', force_constants)
    fc_array.label = 'force_constants'
    return fc_array


def parse_total_dos(filename):
    data = np.loadtxt(filename)
    total_dos = {'frequency_points': data[:, 0],
                 'total_dos': data[:, 1]}
    dos = get_total_dos(total_dos)
    return dos


def parse_projected_dos(filename):
    data = np.loadtxt(filename)
    projected_dos = {'frequency_points': data[:, 0],
                     'projected_dos': data[:, 1:].T}
    pdos = get_projected_dos(projected_dos)
    return pdos


def parse_thermal_properties(filename):
    import yaml

    thermal_properties = {'temperatures': [],
                          'free_energy': [],
                          'entropy': [],
                          'heat_capacity': []}

    with open(filename, 'r') as stream:
        try:
            data = yaml.load(stream, Loader=yaml.FullLoader)
        except AttributeError:
            data = yaml.load(stream)
        for tp in data['thermal_properties']:
            thermal_properties['temperatures'].append(tp['temperature'])
            thermal_properties['entropy'].append(tp['entropy'])
            thermal_properties['free_energy'].append(tp['free_energy'])
            thermal_properties['heat_capacity'].append(tp['heat_capacity'])
        for key in thermal_properties:
            thermal_properties[key] = np.array(thermal_properties[key])

    tprops = get_thermal_properties(thermal_properties)

    return tprops


def parse_band_structure(filename, label=None):
    import yaml

    with open(filename, 'r') as stream:
        try:
            bands = yaml.load(stream, Loader=yaml.FullLoader)
        except AttributeError:
            bands = yaml.load(stream)

    frequencies_flat = []
    qpoints_flat = []
    for k in bands['phonon']:
        frequencies_flat.append([b['frequency'] for b in k['band']])
        qpoints_flat.append(k['q-position'])

    frequencies = []
    qpoints = []
    for n in bands['segment_nqpoint']:
        qpoints.append([qpoints_flat.pop(0) for i in range(n)])
        frequencies.append([frequencies_flat.pop(0) for i in range(n)])

    label_pairs = []
    for pair in bands['labels']:
        label_pairs.append(
            [x.replace('$', '').replace('\\', '').replace('mathrm{', '').replace('}', '').upper()
             for x in pair])

    labels = [label_pairs[0][0], label_pairs[0][1]]
    path_connections = []
    for i, pairs in enumerate(label_pairs[1:]):
        if pairs[0] == label_pairs[i][1]:
            labels.append(pairs[1])
            path_connections.append(True)
        else:
            labels += pairs
            path_connections.append(False)
    path_connections.append(False)

    bs = get_bands(qpoints, frequencies, labels, path_connections, label=label)

    return bs


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
def get_BORN_txt(nac_data, structure, symmetry_tolerance):
    """Returns a string of BORN file.

    nac_data : ArrayData
        Born effective charges and dielectric constants
    structure : StructureData
        This is assumed to be the primitive cell in workchain.
    symmetry_tolerance : float
        Symmetry tolerance.
    """

    from phonopy.file_IO import get_BORN_lines

    born_charges = nac_data.get_array('born_charges')
    epsilon = nac_data.get_array('epsilon')
    pcell = phonopy_atoms_from_structure(structure)
    lines = get_BORN_lines(pcell, born_charges, epsilon,
                           symprec=symmetry_tolerance)

    return "\n".join(lines)


def get_FORCE_SETS_txt(force_sets, displacements_dataset):
    from phonopy.file_IO import get_FORCE_SETS_lines

    forces = force_sets.get_array('force_sets')
    lines = get_FORCE_SETS_lines(displacements_dataset, forces=forces)

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
    lines = [('DIM =' + vals).format(*supercell_matrix)]
    lines.append('PRIMITIVE_AXIS = AUTO')
    mesh = parameters['mesh']
    try:
        length = float(mesh)
        lines.append('MESH = {}'.format(length))
    except TypeError:
        lines.append('MESH = {} {} {}'.format(*mesh))
    lines.append('WRITE_MESH = .FALSE.')

    return '\n'.join(lines)


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

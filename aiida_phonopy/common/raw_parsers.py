from aiida.orm import DataFactory

import numpy as np

ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')

ForceConstantsData = DataFactory('phonopy.force_constants')
PhononDosData = DataFactory('phonopy.phonon_dos')
BandStructureData = DataFactory('phonopy.band_structure')


# PARSE FILES TO DATA

def parse_FORCE_CONSTANTS(filename):
    fcfile = open(filename)
    num = int((fcfile.readline().strip().split())[0])
    force_constants = np.zeros((num, num, 3, 3), dtype=float)
    for i in range(num):
        for j in range(num):
            fcfile.readline()
            tensor = []
            for k in range(3):
                tensor.append([float(x) for x in fcfile.readline().strip().split()])
            force_constants[i, j] = np.array(tensor)
    return ForceConstantsData(data=force_constants)


def parse_partial_DOS(filename, structure, parameters):
    partial_dos = np.loadtxt(filename)

    from phonopy.structure.atoms import Atoms as PhonopyAtoms
    from phonopy import Phonopy

    bulk = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)

    phonon = Phonopy(bulk,
                     supercell_matrix=parameters.dict.supercell,
                     primitive_matrix=parameters.dict.primitive,
                     symprec=parameters.dict.symmetry_precision)

    dos = PhononDosData(frequencies=partial_dos.T[0],
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

    nb = input_bands.get_number_of_bands()

    frequencies = np.array(frequencies)
    print (frequencies.shape)
    frequencies = frequencies.reshape((nb, -1, frequencies.shape[1]))

    band_structure = BandStructureData(frequencies=frequencies,
                                       bands=input_bands.get_bands(),
                                       labels=input_bands.get_labels(),
                                       unitcell=input_bands.get_unitcell())

    return band_structure


# WRITE DATA TO TEXT


def get_BORN_txt(nac_data, symprec=1.e-5, parameters=None, structure=None):
    from phonopy.structure.cells import get_primitive, get_supercell
    from phonopy.structure.symmetry import Symmetry
    from phonopy.structure.atoms import Atoms as PhonopyAtoms

    print ('inside born parameters')

    born_charges = nac_data.get_array('born_charges')
    epsilon = nac_data.get_array('epsilon')
    structure_born = nac_data.get_structure()

    ucell = PhonopyAtoms(symbols=[site.kind_name for site in structure_born.sites],
                         positions=[site.position for site in structure_born.sites],
                         cell=structure_born.cell)

    if structure is not None:
        pmat = parameters.dict.primitive
        inv_pmat = np.linalg.inv(pmat)
        inv_smat = np.linalg.inv(structure.cell)
        rmat = np.dot(inv_smat, structure_born.cell)
        pmat = np.dot(inv_pmat, rmat)
        smat = parameters.dict.supercell

    else:
        # NAC parameters are assumed to be calculated using the primitive cell defined in phonopy parameters
        pmat = np.identity(3)
        smat = np.identity(3)

    num_atom = len(born_charges)
    assert num_atom == ucell.get_number_of_atoms(), \
        "num_atom %d != len(borns) %d" % (ucell.get_number_of_atoms(),
                                          len(born_charges))

    inv_smat = np.linalg.inv(smat)
    scell = get_supercell(ucell, smat, symprec=symprec)
    pcell = get_primitive(scell, np.dot(inv_smat, pmat), symprec=symprec)

    p2s = np.array(pcell.get_primitive_to_supercell_map(), dtype='intc')
    p_sym = Symmetry(pcell, is_symmetry=True, symprec=symprec)
    s_indep_atoms = p2s[p_sym.get_independent_atoms()]
    u2u = scell.get_unitcell_to_unitcell_map()
    u_indep_atoms = [u2u[x] for x in s_indep_atoms]
    reduced_borns = born_charges[u_indep_atoms].copy()

    # factor = get_default_physical_units('vasp')['nac_factor']  # born charges in VASP units
    factor = 14.0
    born_txt = ('{}\n'.format(factor))
    for num in epsilon.flatten():
        born_txt += ('{0:4.8f} '.format(num))
    born_txt += ('\n')

    for atom in reduced_borns:
        for num in atom.flatten():
            born_txt += ('{0:4.8f} '.format(num))
        born_txt += ('\n')
    # born_txt += ('{}\n'.format(factor))

    return born_txt


def get_FORCE_SETS_txt(data_sets_object):
    data_sets = data_sets_object.get_force_sets()

    #    data_list = []
    #    for name in names:
    #        data_list.append({'direction': data_sets_object.get_array(name)[0],
    #                          'number': data_sets_object.get_array(name)[1],
    #                          'displacement': data_sets_object.get_array(name)[2],
    #                          'forces': data_sets_object.get_array(name)[3]})
    #    data_sets = {'natom': num_atom, 'first_atoms': data_list}

    displacements = data_sets['first_atoms']
    forces = [x['forces'] for x in data_sets['first_atoms']]

    # Write FORCE_SETS
    force_sets_txt = "%-5d\n" % data_sets['natom']
    force_sets_txt += "%-5d\n" % len(displacements)
    for count, disp in enumerate(displacements):
        force_sets_txt += "\n%-5d\n" % (disp['number'] + 1)
        force_sets_txt += "%20.16f %20.16f %20.16f\n" % (tuple(disp['displacement']))

        for f in forces[count]:
            force_sets_txt += "%15.10f %15.10f %15.10f\n" % (tuple(f))
    return force_sets_txt


def get_FORCE_CONSTANTS_txt(force_constants_object):

    force_constants = force_constants_object.get_data()

    # Write FORCE CONSTANTS
    force_constants_txt = '{0}\n'.format(len(force_constants))
    for i, fc in enumerate(force_constants):
        for j, atomic_fc in enumerate(fc):
            force_constants_txt += '{0} {1}\n'.format(i, j)
            for line in atomic_fc:
                force_constants_txt += '{0:20.16f} {1:20.16f} {2:20.16f}\n'.format(*line)

    return force_constants_txt


def get_poscar_txt(structure):
    types = [site.kind_name for site in structure.sites]
    atom_type_unique = np.unique(types, return_index=True)
    sort_index = np.argsort(atom_type_unique[1])
    elements = np.array(atom_type_unique[0])[sort_index]
    elements_count = np.diff(np.append(np.array(atom_type_unique[1])[sort_index], [len(types)]))

    poscar_txt = '# VASP POSCAR generated using aiida workflow '
    poscar_txt += '\n1.0\n'
    cell = structure.cell
    for row in cell:
        poscar_txt += '{0: 22.16f} {1: 22.16f} {2: 22.16f}\n'.format(*row)
    poscar_txt += ' '.join([str(e) for e in elements]) + '\n'
    poscar_txt += ' '.join([str(e) for e in elements_count]) + '\n'
    poscar_txt += 'Cartesian\n'
    for site in structure.sites:
        poscar_txt += '{0: 22.16f} {1: 22.16f} {2: 22.16f}\n'.format(*site.position)

    return poscar_txt


def get_phonopy_conf_file_txt(parameters_object, bands=None):
    parameters = parameters_object.get_dict()
    input_file = 'DIM = {} {} {}\n'.format(*np.diag(parameters['supercell']))
    input_file += 'PRIMITIVE_AXIS = {} {} {}  {} {} {}  {} {} {}\n'.format(
        *np.array(parameters['primitive']).reshape((1, 9))[0])
    input_file += 'MESH = {} {} {}\n'.format(*parameters['mesh'])

    print bands
    if bands is not None:
        input_file += 'BAND = '
        for i, band in enumerate(bands.get_band_ranges()):
            input_file += ' '.join(np.array(band.flatten(), dtype=str))
            if i < bands.get_number_of_bands() - 1:
                input_file += ','

        input_file += '\n'
        input_file += 'BAND_POINTS = {}'.format(bands.get_number_of_points())

    return input_file


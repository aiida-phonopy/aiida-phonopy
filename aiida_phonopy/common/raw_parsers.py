import numpy as np
import yaml
from aiida.plugins import DataFactory
from aiida_phonopy.common.utils import (
    get_total_dos, get_projected_dos, get_thermal_properties, get_bands)
from phonopy.file_IO import parse_FORCE_CONSTANTS as parse_FC


ArrayData = DataFactory('array')
XyData = DataFactory('array.xy')
BandsData = DataFactory('array.bands')


# Parse files to generate AIIDA OBJECTS

def parse_FORCE_CONSTANTS(filename):
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
    import h5py

    kappa = ArrayData()

    f = h5py.File(filename, 'r')
    for item in f:
        array = np.array(f[item])
        kappa.set_array(item, array)
    f.close()

    return kappa

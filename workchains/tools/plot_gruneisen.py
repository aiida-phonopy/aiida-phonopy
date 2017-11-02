from aiida import load_dbenv
load_dbenv()

from aiida.orm import load_node, load_workflow
from aiida.orm import Code, DataFactory

import matplotlib.pyplot as plt
from matplotlib import gridspec

StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
KpointsData = DataFactory('array.kpoints')

import sys
import numpy as np

if len(sys.argv) < 2:
    print ('use: plot_gruneisen pk_number')
    exit()

# Set WorkChain GruneisenPhonopy PK number
################################
wc = load_node(int(sys.argv[1]))
################################

gamma_cutoff = 0.01

# Phonon Band structure
bs = wc.out.band_structure

labels, indices, widths, ranges, positions = bs.get_plot_helpers()
gs = gridspec.GridSpec(1, len(widths), width_ratios=widths, wspace=0.05)

plt.figure(1, figsize=(8, 5))

plt.rcParams.update({'mathtext.default': 'regular'})

for j, index in enumerate(indices):
    ax1 = plt.subplot(gs[j])

    plt.gca().set_color_cycle(None)
    for i in index:
        ax1.plot(bs.get_distances(band=i),
                 bs.get_frequencies(band=i),
                 #color='r'
                 )
    if j != 0:
        ax1.axes.get_yaxis().set_visible(False)
        ax1.set_ylim(ylim)

    ylim = ax1.get_ylim()

    plt.axhline(y=0.0, color='black', linestyle='--', linewidth=0.1)
    plt.ylabel('Frequency (THz)')
    plt.xlim(ranges[j])
    plt.xticks(positions[j], labels[j], rotation='horizontal')

plt.suptitle('Phonon band structure')
plt.figtext(0.5, 0.02, 'Wave vector', ha='center')

# Mode Gruneisen parameter band structure
plt.figure(2, figsize=(8, 5))

plt.rcParams.update({'mathtext.default': 'regular'})

for j, index in enumerate(indices):
    ax1 = plt.subplot(gs[j])

    plt.gca().set_color_cycle(None)
    for i in index:
        q_points = bs.get_bands(band=i)
        mask = np.where(np.linalg.norm(q_points, axis=1) > gamma_cutoff)

        ax1.plot(bs.get_distances(band=i)[mask],
                 bs.get_gamma(band=i)[mask],
                 # color='r'
                 )

    if j != 0:
        ax1.axes.get_yaxis().set_visible(False)
        ax1.set_ylim(ylim)

    ylim = ax1.get_ylim()

    plt.axhline(y=0.0, color='black', linestyle='--', linewidth=0.1)
    plt.ylabel('$\gamma$')
    plt.xlim(ranges[j])
    plt.xticks(positions[j], labels[j], rotation='horizontal')

plt.suptitle('Mode Gruneisen parameter')
plt.figtext(0.5, 0.02, 'Wave vector', ha='center')


# Mode Gruneisen parameter mesh (freq vs gamma)
mesh = wc.out.mesh

plt.figure(3)
q_points = mesh.get_array('q_points')
mask = np.where(np.linalg.norm(q_points, axis=1) > gamma_cutoff)

for gamma, freq in zip( mesh.get_array('gruneisen').T,
                        mesh.get_array('frequencies').T):
    plt.plot(freq[mask], gamma[mask],
             marker='o',
             linestyle='None',
             markeredgecolor='black',
             # color='red'
             )
plt.xlabel('Frequency [THz]')
plt.ylabel('$\gamma$')
plt.title('Mode Gruneisen parameter (mesh)')

plt.show()

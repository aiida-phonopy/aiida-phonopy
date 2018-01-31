###################################################
# This script plots the thermal conductivity      #
# calculated using thermal_conductivity WorkChain #
# To run this script use the pk number:           #
# $ python plot_thermal_conductivity.py pknumber  #
###################################################

from aiida import load_dbenv
load_dbenv()

from aiida.orm import load_node, load_workflow
from aiida.orm import Code, DataFactory

from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

import sys

if len(sys.argv) < 2:
    print ('use: python plot_phonon.py {pk_number}')
    exit()

# Set WorkChain PhononPhonopy PK number
################################
wc = load_node(int(sys.argv[1]))
################################


def kdeplot(fig, array_data, xmax=None, ymax=None, ymin=0, zmax=None,
            ylabel=True, title=None,
            cmap='rainbow', cutoff=None, density_ratio=0.1, fmax=None,
            gv_norm=False, nbins=100, nu=False, temperature=300,
            scatter=True, smear=True, color='w'):
    """

      description   :
        cutoff --- if None, 90 % of gamma is plotted
                   warning: it doesn't consider about weight
                   the mesh points whose lifetime > cutoff is ignored.
        fmax --- you don't have to set this parameter
        scatter --- if True, ax.scatter is evoked
        color --- the color of plot
        smear --- if True, color map is evoked
    """

    ax = fig.add_axes((0.1, 0.15, 0.8, 0.7))
    epsilon = 1.0e-8

    def collect_data(gamma, weights, frequencies,
                     t_index, cutoff, max_freq):
        freqs = []
        mode_prop = []

        for w, freq, g in zip(weights, frequencies, gamma):

            tau = 1.0 / np.where(g > 10e-10, g, -1) / (2 * 2 * np.pi)

            if not cutoff:
                idx = int(len(gamma.flatten()) * 0.1)
                cutoff = 1.0 / np.sort(gamma.flatten())[idx]\
                             / (2 * 2 * np.pi)

            tau = np.where(tau < cutoff, tau, -1)
            condition = tau > 0
            _tau = np.extract(condition, tau)
            _freq = np.extract(condition, freq)

            if max_freq is None:
                freqs += list(_freq) * w
                mode_prop += list(_tau) * w
            else:
                freqs += list(np.extract(freq < max_freq, freq))
                mode_prop += list(np.extract(freq < max_freq, tau))
        x = np.array(mode_prop)
        y = np.array(freqs)

        return x, y

    def run_KDE(x, y, nbins, x_max=None, y_max=None, density_ratio=0.1):
        # Running Gaussian-KDE by scipy
        x_min = 0
        if x_max is None:
            _x_max = np.rint(x.max())
        else:
            _x_max = x_max
        y_min = 0
        if y_max is None:
            _y_max = np.rint(y.max() * 1.1)
        else:
            _y_max = y_max
        values = np.vstack([x.ravel(), y.ravel()])
        kernel = stats.gaussian_kde(values)

        xi, yi = np.mgrid[x_min:_x_max:nbins * 1j, y_min:_y_max:nbins * 1j]
        positions = np.vstack([xi.ravel(), yi.ravel()])
        zi = np.reshape(kernel(positions).T, xi.shape)

        if x_max is None:
            zi_max = np.max(zi)
            indices = []
            for i, r_zi in enumerate((zi.T)[::-1]):
                if indices:
                    indices.append(nbins - i - 1)
                elif np.max(r_zi) > zi_max * density_ratio:
                    indices = [nbins - i - 1]
            short_nbinds = len(indices)
            xnbins = nbins ** 2 // short_nbinds
            xi, yi = np.mgrid[x_min:_x_max:xnbins * 1j,
                              y_min:_y_max:nbins * 1j]
            positions = np.vstack([xi.ravel(), yi.ravel()])
            zi = np.reshape(kernel(positions).T, xi.shape)
        else:
            short_nbinds = nbins

        return xi, yi, zi, short_nbinds

    def plot(ax, xi, yi, zi, x, y, short_nbinds, nbins,
             x_max=None, z_max=None, cmap=None, color=color):
        #
        # Plotting
        #
        xmax = np.max(x)
        ymax = np.max(y)
        x_cut = []
        y_cut = []
        threshold = xmax / nbins * short_nbinds / nbins * (nbins - 1)
        for _x, _y in zip(x, y):
            if epsilon < _x and _x < threshold and epsilon < _y and _y < ymax - epsilon:
                x_cut.append(_x)
                y_cut.append(_y)

        if smear:
            plt.pcolormesh(xi[:, :nbins], yi[:, :nbins], zi[:, :nbins],
                           vmax=z_max, cmap=cmap)
            plt.colorbar()

        if scatter:
            ax.scatter(x_cut, y_cut, s=5, c=color, marker='.', linewidth=0)

        ax.set_ylim(ymin=ymin, ymax=yi.max())
        if x_max is None:
            fig_x_lst = list(x_cut)
            ax.set_xlim(xmin=0, xmax=(max(fig_x_lst) + epsilon))
        else:
            ax.set_xlim(xmin=0, xmax=(x_max + epsilon))

        ax.set_xlabel('Lifetime (ps)')

        if ylabel:
            ax.set_ylabel('Phonon frequency (THz)')

        else:
            ax.set_yticklabels([])
            ax.set_ylabel('')

        ax.plot([xmax, 0], [0, 0], color="black", linestyle=":")


    temperatures = array_data.get_array('temperature')[:]
    weights = array_data.get_array('weight')[:]
    frequencies = array_data.get_array('frequency')[:]
    group_velocity = array_data.get_array('group_velocity')[:,:,:]

    if len(temperatures) > 29:
        t_index = 30
    else:
        t_index = 0
    for i, t in enumerate(temperatures):
        if np.abs(t - temperature) < epsilon:
            t_index = i
            break
    plt.title('Temperature: {} K'.format(temperatures[t_index]))

    symbols = ['']
    if gv_norm:
        gv_norm =\
            np.sqrt((group_velocity ** 2).sum(axis=2))
        gammas = [gv_norm]
    else:
        gammas = [array_data.get_array('gamma')[t_index]]
        if nu:
            if 'gamma_N' in array_data.get_arraynames():
                gammas.append(array_data.get_array('gamma_N')[t_index])
                symbols.append('N')
            if 'gamma_U' in array_data.get_arraynames():
                gammas.append(array_data.get_array('gamma_U')[t_index])
                symbols.append('U')

    for gamma, s in zip(gammas, symbols):
        x, y = collect_data(gamma, weights, frequencies,
                            t_index, cutoff, fmax)
        xi, yi, zi, short_nbinds = run_KDE(x, y, nbins,
                                           x_max=xmax,
                                           y_max=ymax,
                                           density_ratio=density_ratio)
        plot(ax, xi, yi, zi, x, y, short_nbinds, nbins,
             x_max=xmax, z_max=zmax, cmap=cmap)

    ax.plot([xmax, 0], [0, 0], color="black", linestyle=":")


kappa = wc.out.kappa
print kappa.get_arraynames()
#fig = plt.figure(0)
#kdeplot(fig, kappa, temperature=300, ylabel=True, cutoff=None,
#        nbins=50, scatter=True, smear=True, color='black')

fig = plt.figure(0)
frequency = wc.out.kappa.get_array('frequency')
gamma = wc.out.kappa.get_array('gamma')
weight =  wc.out.kappa.get_array('weight')
temperature = wc.out.kappa.get_array('temperature')

print frequency.shape
print gamma.shape
print temperature.shape


def get_tindex(temperatures, temperature=300, epsilon=1.0e-8):
    if len(temperatures) > 29:
        t_index = 30
    else:
        t_index = 0
    for i, t in enumerate(temperatures):
        if np.abs(t - temperature) < epsilon:
            t_index = i
            break
    return t_index

t_index = get_tindex(temperature, temperature=300)
plt.scatter(gamma[t_index], frequency, s=weight)


plt.show()


from phonopy.phonon.tetrahedron_mesh import TetrahedronMesh

epsilon = 1.0e-8

class KappaDOS(object):
    def __init__(self,
                 mode_kappa,
                 cell,
                 frequencies,
                 mesh,
                 grid_address,
                 grid_mapping_table,
                 ir_grid_points,
                 grid_order=None,
                 num_sampling_points=100):
        self._mode_kappa = mode_kappa
        self._tetrahedron_mesh = TetrahedronMesh(
            cell,
            frequencies,
            mesh,
            grid_address,
            grid_mapping_table,
            ir_grid_points)

        min_freq = min(frequencies.ravel())
        max_freq = max(frequencies.ravel()) + epsilon
        self._frequency_points = np.linspace(min_freq,
                                             max_freq,
                                             num_sampling_points)
        self._kdos = np.zeros(
            (len(mode_kappa), len(self._frequency_points), 2, 6),
            dtype='double')
        self._run_tetrahedron_method()

    def get_kdos(self):
        return self._frequency_points, self._kdos

    def _run_tetrahedron_method(self):
        num_freqs = len(self._frequency_points)
        thm = self._tetrahedron_mesh
        for j, value in enumerate(('J', 'I')):
            thm.set(value=value, frequency_points=self._frequency_points)
            for i, iw in enumerate(thm):
                print iw.shape
                print mode_kappa.shape
                # kdos[temp, freq_points, IJ, tensor_elem]
                # iw[freq_points, band]
                # mode_kappa[temp, ir_gp, band, tensor_elem]
                self._kdos[:, :, j] += np.transpose(
                    np.dot(iw, self._mode_kappa[:, i]), axes=(1, 0, 2))


mode_kappa = wc.out.kappa.get_array('mode_kappa')[:]
mesh = wc.out.kappa.get_array('mesh')
qpoint = wc.out.kappa.get_array('qpoint')

print mode_kappa.shape

structure = wc.out.final_structure

from phono3py.phonon3.triplets import (get_ir_grid_points,
                                       get_grid_points_by_rotations,
                                       get_grid_address,
                                       from_coarse_to_dense_grid_points)

# May change
grid_address = get_grid_address(mesh)
ir_grid_points = np.arange(np.prod(mesh), dtype='intc')
grid_mapping_table = np.arange(np.prod(mesh), dtype='intc')


from aiida_phonopy.workchains.phonon import phonopy_bulk_from_structure
kappados = KappaDOS(mode_kappa=mode_kappa,
                    cell=phonopy_bulk_from_structure(structure),
                    frequencies=frequency,
                    mesh=mesh,
                    grid_address=grid_address,
                    ir_grid_points=ir_grid_points,
                    grid_mapping_table=grid_mapping_table,
                    num_sampling_points=100)

print kappados.get_kdos()

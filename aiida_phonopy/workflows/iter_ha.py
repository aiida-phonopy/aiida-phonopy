import numpy as np
from phonopy import Phonopy
from phonopy.structure.dataset import get_displacements_and_forces
from aiida.engine import WorkChain
from aiida.plugins import WorkflowFactory, DataFactory
from aiida.orm import Float, Int, QueryBuilder, Group, load_node
from aiida.engine import while_, if_, calcfunction
from aiida_phonopy.common.utils import phonopy_atoms_from_structure

Dict = DataFactory('dict')
PhonopyWorkChain = WorkflowFactory('phonopy.phonopy')


def _collect_dataset(data, linear_decay=True):
    """Collect supercell displacements, forces, and energies

    With linear_decay=True, numbers of snapshots to be taken
    are biased. Older snapshots are taken lesser. The fraction
    of the number of snapshots in each previous phonon calculation
    (many snapshots in one phonon calculation) to be taken is defined
    linearly by

        ratios = (np.arange(max_items, dtype=float) + 1) / max_items

    where max_items is the number of previous phonon calculations to be
    included at maximum.

    """

    nitems = max([int(key.split('_')[-1])
                  for key in data.keys() if 'forces' in key])
    max_items = data['max_previous_steps'].value
    displacements = []
    forces = []
    energies = []

    if linear_decay:
        ratios = (np.arange(max_items, dtype=float) + 1) / max_items
    else:
        ratios = np.ones(max_items, dtype=int)
    ratios = ratios[-nitems:]

    for i in range(nitems):
        force_sets = data['forces_%d' % (i + 1)].get_array('force_sets')
        phonon_setting_info = data['ph_info_%d' % (i + 1)]
        dataset = phonon_setting_info['displacement_dataset']
        disps, _ = get_displacements_and_forces(dataset)

        n_include = int(ratios[i] * len(force_sets)) + 1
        if n_include < len(force_sets):
            n_include += 1
        if len(disps) < n_include:
            n_include = len(disps)

        forces.append(force_sets[:n_include])
        if 'energies' in data['forces_%d' % (i + 1)].get_arraynames():
            energy_sets = data['forces_%d' % (i + 1)].get_array('energies')
            energies.append(energy_sets[:n_include])
        displacements.append(disps[:n_include])

    return displacements, forces, energies


def _remove_high_energy_snapshots(d, f, e, ratio):
    """Remove snapshots that have high energies with a given ratio"""

    num_include = int(np.ceil(ratio * len(e)))
    if num_include > len(e):
        num_include = len(e)
    idx = np.argsort(e)[:num_include]
    d = d[idx]
    f = f[idx]
    return d, f, e, idx


def _modify_force_constants(ph):
    # Trick (push up low frequency modes)
    # Create RandomDisplacements class instance
    # temperature=300 can be any number just to invoke finite temperature
    # random_displacements.
    ph.generate_displacements(number_of_snapshots=1, temperature=300)
    rd = ph.random_displacements  # Get RandomDisplacements class instance
    freqs = np.abs(rd.frequencies)  # Make frequencies have absolute values

    # As default phonon frequencies < 0.01 are ignored.
    # For phonon modes where 0.01 < |freqs| < 0.5  --> |freqs| + 1
    condition = np.logical_and(freqs > 0.01, freqs < 0.5)
    # if condition.any():
    #     self.report("Some phonon frequencies are very low.")

    freqs = np.where(condition, freqs + 1, freqs)
    # Set modified frequencies to RandomDisplacements class instance
    rd.frequencies = freqs
    rd.run_d2f()  # Create force constants using modified frequencies
    # Set modified force constants to Phonopy class instance
    ph.force_constants = rd.force_constants


@calcfunction
def get_random_displacements(structure,
                             number_of_snapshots,
                             temperature,
                             **data):
    """


    """

    displacements, forces, energies = _collect_dataset(data)

    # Concatenate the data
    d = np.concatenate(displacements, axis=0)
    f = np.concatenate(forces, axis=0)
    if energies:
        e = np.concatenate(energies)
    else:
        e = None

    # Remove snapshots that have high energies when include_ratio is given.
    idx = None
    if e is not None and len(e) == len(f) and 'include_ratio' in data:
        ratio = data['include_ratio'].value
        if 0 < ratio and ratio < 1:
            d, f, e, idx = _remove_high_energy_snapshots(d, f, e, ratio)

    # Calculate force constants by fitting using ALM
    phonon_setting_info = data['ph_info_1']
    smat = phonon_setting_info['supercell_matrix']
    ph = Phonopy(phonopy_atoms_from_structure(structure),
                 supercell_matrix=smat,
                 primitive_matrix='auto')
    ph.dataset = {'displacements': d, 'forces': f}
    ph.produce_force_constants(fc_calculator='alm')

    # Treatment of imaginary modes
    _modify_force_constants(ph)

    # Generate random displacements at a given temperature
    if 'random_seed' in data:
        _random_seed = data['random_seed'].value
    else:
        _random_seed = None
    ph.generate_displacements(
        number_of_snapshots=number_of_snapshots.value,
        random_seed=_random_seed,
        temperature=temperature.value)

    ret_dict = {'displacement_dataset': Dict(dict=ph.dataset)}
    e_dict = {'supercell_energies': energies}
    if idx is not None:
        e_dict['included_supercell_indices'] = idx.tolist()
    ret_dict['supercell_energies'] = Dict(dict=e_dict)

    return ret_dict


class IterHarmonicApprox(WorkChain):
    """ Workchain for harmonic force constants by iterative approach

    By default, the calculation starts with normal phonon calculation,
    i.e., in this context, which corresponds to roughly 0K force constants.
    Then the iteration loop starts. The first run is the iteration-1.
    The iteration stops after finishing that of max_iteration. Each phonon
    calculation is named 'step'.

    Steps
    -----
    0. Initial phonon calculation at 0K
    1. First phonon calculation at specified temperature. Random
       displacements are created from step-0.
    2. Second phonon calculation at specified temperature. Random
       displacements are created from step-1.
    3. Third phonon calculation at specified temperature. Random
       displacements are created from steps 1 and 2 if
       number_of_snapshots >= 2. Otherwise only the result from
       step-2 is used.
    4. Four phonon calculation at specified temperature. Random
       displacements are created from number_of_snapshots previous
       existing steps excluding step-0.
    *. Continue until iteration number = max_iteration number.

    Manual termination of iteration loop
    ------------------------------------
    It is possible to terminate at the initial point of each iteration.
    This option is not very recommended to use because not reproducible
    mechanically, but can be useful for experimental calculation.

    This is achieved by just creating AiiDA Group whose label is its
    uuid string that ejected by AiiDA, i.e., self.uuid.

    inputs
    ------
    Most of inputs are imported from PhonopyWorkChain. Specific inputs
    of this workchain are as follows:

    max_iteration : Int
        Maximum number of iterations.
    number_of_snapshots : Int
        Number of generated supercell snapshots with random displacements
        at a temperature.
    number_of_steps_for_fitting : Int
        Displacements and respective forces of supercells in the previous
        number_of_steps_for_fitting are used to simultaneously fit to
        force constants.
    random_seed : Int
        Random seed used to sample in canonical ensemble harmonic oscillator
        space. The value must be 32bit unsigned int. Unless specified,
        random seed will not be fixed.
    temperature : Float
        Temperature (K).
    initial_nodes : Dict, optional
        This gives the initial nodes that contain sets of forces, which are
        provided by PKs or UUIDs.
    include_ratio : Float
        How much supercell forces are included from lowest supercell energies.

    """

    @classmethod
    def define(cls, spec):
        super(IterHarmonicApprox, cls).define(spec)
        spec.expose_inputs(PhonopyWorkChain,
                           exclude=['immigrant_calculation_folders',
                                    'calculation_nodes', 'dry_run'])
        spec.input('max_iteration', valid_type=Int, required=False,
                   default=lambda: Int(10))
        spec.input('number_of_snapshots', valid_type=Int, required=False,
                   default=lambda: Int(100))
        spec.input('number_of_steps_for_fitting', valid_type=Int,
                   required=False, default=lambda: Int(4))
        spec.input('random_seed', valid_type=Int, required=False)
        spec.input('temperature', valid_type=Float, required=False,
                   default=lambda: Float(300.0))
        spec.input('initial_nodes', valid_type=Dict, required=False)
        spec.input('include_ratio', valid_type=Float, required=False)
        spec.outline(
            cls.initialize,
            if_(cls.import_initial_nodes)(
                cls.set_initial_nodes,
                cls.run_phonon,
            ).else_(
                cls.run_initial_phonon,
            ),
            while_(cls.is_loop_finished)(
                cls.run_phonon,
            ),
        )

    def import_initial_nodes(self):
        return 'initial_nodes' in self.inputs

    def initialize(self):
        self.report("initialize (%s)" % self.uuid)
        self.ctx.iteration = 0
        self.ctx.prev_nodes = []

    def is_loop_finished(self):
        qb = QueryBuilder()
        qb.append(Group, filters={'label': {'==': self.uuid}})
        if qb.count() == 1:
            self.report("Iteration loop is manually terminated at step %d."
                        % self.ctx.iteration)
            return False

        self.ctx.iteration += 1
        return self.ctx.iteration <= self.inputs.max_iteration.value

    def set_initial_nodes(self):
        self.report("set_initial_phonon")
        node_ids = self.inputs.initial_nodes['nodes']
        self.ctx.prev_nodes = [load_node(node_id) for node_id in node_ids]
        self.ctx.iteration = 1

    def run_initial_phonon(self):
        self.report("run_initial_phonon")
        inputs = self._get_phonopy_inputs(None, None)
        inputs['metadata'].label = "Initial phonon calculation"
        inputs['metadata'].description = "Initial phonon calculation"
        future = self.submit(PhonopyWorkChain, **inputs)
        self.ctx.initial_node = future
        self.report('{} pk = {}'.format(inputs['metadata'].label, future.pk))
        self.to_context(**{'phonon_0': future})

    def run_phonon(self):
        self.report("run_phonon_%d" % self.ctx.iteration)

        n_ave = self.inputs.number_of_steps_for_fitting.value
        if len(self.ctx.prev_nodes) == 0:
            nodes = [self.ctx.initial_node, ]
        elif len(self.ctx.prev_nodes) < n_ave:
            nodes = self.ctx.prev_nodes
        else:
            nodes = self.ctx.prev_nodes[-n_ave:]
        dataset = self._set_ave_fc(nodes)
        inputs = self._get_phonopy_inputs(dataset, False)
        label = "Phonon calculation %d" % self.ctx.iteration
        inputs['metadata'].label = label
        inputs['metadata'].description = label
        future = self.submit(PhonopyWorkChain, **inputs)
        self.ctx.prev_nodes.append(future)
        self.report('{} pk = {}'.format(label, future.pk))
        label = "phonon_%d" % self.ctx.iteration
        self.to_context(**{label: future})

    def _set_ave_fc(self, nodes):
        data = {}
        for i, node in enumerate(nodes):
            data['forces_%d' % (i + 1)] = node.outputs.force_sets
            data['ph_info_%d' % (i + 1)] = node.outputs.phonon_setting_info

        if 'random_seed' in self.inputs:
            data['random_seed'] = self.inputs.random_seed
            self.report("Random seed is %d." % self.inputs.random_seed.value)
        if 'include_ratio' in self.inputs:
            data['include_ratio'] = self.inputs.include_ratio
            self.report("Include ratio is %f."
                        % self.inputs.include_ratio.value)
        data['max_previous_steps'] = self.inputs.number_of_steps_for_fitting

        displacements = get_random_displacements(
            nodes[-1].inputs.structure,
            self.inputs.number_of_snapshots,
            self.inputs.temperature,
            **data)

        return displacements['displacement_dataset']

    def _get_phonopy_inputs(self, dataset, is_nac):
        inputs_in = self.exposed_inputs(PhonopyWorkChain)
        inputs = inputs_in.copy()
        phonon_settings = inputs_in['phonon_settings'].get_dict()
        if is_nac is not None:
            phonon_settings['is_nac'] = is_nac
        inputs['phonon_settings'] = Dict(dict=phonon_settings)
        if dataset is not None:
            inputs['displacement_dataset'] = dataset
        return inputs

import numpy as np
from phonopy import Phonopy
from phonopy.structure.dataset import get_displacements_and_forces
from aiida.engine import WorkChain
from aiida.plugins import WorkflowFactory, DataFactory
from aiida.orm import Bool, Float, Int, QueryBuilder, Group, load_node, Code
from aiida.engine import while_, if_, calcfunction
from aiida_phonopy.common.utils import phonopy_atoms_from_structure

Dict = DataFactory('dict')
ArrayData = DataFactory('array')
PhonopyWorkChain = WorkflowFactory('phonopy.phonopy')

"""

Parameters in get_random_displacements and collect_dataset
----------------------------------------------------------
To control how to include previous snapshots, several parameters
can be used.

1) number_of_steps_for_fitting : Int
    Maximum number of previous phonon calculation steps included. When
    number of the previous phonon calculations is smaller than this
    number, all the previous phonon calculations are included. But
    this number is used for (2) in any case.
2) linear_decay : True
    This will be an option, but currently always True.
    This controls weights of previous phonon calculations. One previous
    phonon calculation is included 100%. (number_of_steps_for_fitting + 1)
    previous phonon calculation is not include. Between them, those are
    included with linear scalings. The snapshots in each phonon calculation
    are included from the first element of the list to the specified
    number, i.e., [:num_included_snapshots].
3) include_ratio : Float
    After collecting snapshots in (2), all the snapshots are sorted by
    total energies. Then only lowest energy snapshots with 'include_ratio'
    are included to calculate force constants by fitting using ALM.
4) random_seed : Int
    Using force constants created in (3), phonons are calculated at
    commensurate points, and these phonons are used to generate atomic
    displacements by sampling harmonic oscillator distribution function.
    For this 'random_seed' is used when provided.

"""


def get_random_displacements(structure,
                             number_of_snapshots,
                             temperature,
                             phonon_setting_info,
                             dataset_for_fc,
                             random_seed=None):
    """Generate supercells with random displacemens

    The random displacements are generated from phonons and harmonic
    oscillator distribution function of canonical ensemble. The input
    phonons are calculated from force constants calculated from
    forces and displacemens of the supercell snapshots in previous
    phonon calculation steps.

    Returns
    -------
    dataset : Dict
        Displacement datasets to run force calculations.

    """

    # Calculate force constants by fitting using ALM
    smat = phonon_setting_info['supercell_matrix']
    ph = Phonopy(phonopy_atoms_from_structure(structure),
                 supercell_matrix=smat,
                 primitive_matrix='auto')
    d = dataset_for_fc.get_array('displacements')
    f = dataset_for_fc.get_array('forces')
    ph.dataset = {'displacements': d, 'forces': f}
    ph.produce_force_constants(fc_calculator='alm')

    if random_seed is not None:
        _random_seed = random_seed.value
    else:
        _random_seed = None
    dataset = _generate_random_displacements(ph,
                                             number_of_snapshots.value,
                                             temperature.value,
                                             random_seed=_random_seed)
    return Dict(dict=dataset)


@calcfunction
def collect_dataset(number_of_steps_for_fitting,
                    include_ratio,
                    linear_decay,
                    **data):
    """Collect supercell displacements, forces, and energies

    Returns
    -------
    dataset : ArrayData
        Displacements and forces used for calculating force constants
        to generate random displacements.
    supercell_energies : Dict
        'supercell_energies' : List of list
            Total energies of snapshots before step (3) above.
        'included' : List of list
            Finally included snapshots after step (4). The snapshots are
            indexed by concatenating list of list in step (3).

    """

    nitems = max([int(key.split('_')[-1])
                  for key in data.keys() if 'forces' in key])

    forces_in_db = []
    ph_info_in_db = []
    for i in range(nitems):
        force_sets = data['forces_%d' % (i + 1)]
        phonon_setting_info = data['phonon_setting_info_%d' % (i + 1)]
        forces_in_db.append(force_sets)
        ph_info_in_db.append(phonon_setting_info)

    d, f, energies = _extract_dataset_from_db(forces_in_db, ph_info_in_db)

    displacements, forces, included = _create_dataset(
        d, f, energies,
        number_of_steps_for_fitting.value,
        include_ratio.value,
        linear_decay.value)
    dataset = ArrayData()
    dataset.set_array('forces', forces)
    dataset.set_array('displacements', displacements)
    supercell_energies = Dict(dict={'energies': energies,
                                    'included': included})

    return {'dataset': dataset, 'supercell_energies': supercell_energies}


def _extract_dataset_from_db(forces_in_db, ph_info_in_db):
    nitems = len(forces_in_db)
    displacements = []
    forces = []
    energies = []

    for i in range(nitems):
        force_sets = forces_in_db[i].get_array('force_sets')
        dataset = ph_info_in_db[i]['displacement_dataset']
        disps, _ = get_displacements_and_forces(dataset)

        forces.append(force_sets)
        if 'energies' in forces_in_db[i].get_arraynames():
            energy_sets = forces_in_db[i].get_array('energies')
            energies.append(energy_sets)
        displacements.append(disps)

    return displacements, forces, energies


def _create_dataset(displacements, forces, energies,
                    max_items, ratio, linear_decay):
    included = _choose_snapshots_by_linear_decay(
        displacements, forces, max_items, linear_decay=linear_decay)

    # Remove snapshots that have high energies when include_ratio is given.
    if energies is not None and ratio is not None:
        if 0 < ratio and ratio < 1:
            included = _remove_high_energy_snapshots(energies, included, ratio)

    _displacements, _forces, _energies = _include_snapshots(
        displacements, forces, energies, included)

    # Concatenate the data
    d = np.concatenate(_displacements, axis=0)
    f = np.concatenate(_forces, axis=0)

    return d, f, included


def _choose_snapshots_by_linear_decay(displacements, forces, max_items,
                                      linear_decay=True):
    """Choose snapshots by linear_decay

    With linear_decay=True, numbers of snapshots to be taken
    are biased. Older snapshots are taken lesser. The fraction
    of the number of snapshots in each previous phonon calculation
    (many snapshots in one phonon calculation) to be taken is defined
    linearly by

        ratios = (np.arange(max_items, dtype=float) + 1) / max_items

    where max_items is the number of previous phonon calculations to be
    included at maximum.

    """

    assert len(forces) == len(displacements)

    nitems = len(forces)

    if linear_decay:
        ratios = (np.arange(max_items, dtype=float) + 1) / max_items
    else:
        ratios = np.ones(max_items, dtype=int)
    ratios = ratios[-nitems:]
    included = []

    for i in range(nitems):
        n = len(forces[i])

        assert n == len(displacements[i])

        n_in = int(ratios[i] * n + 0.5)
        if n < n_in:
            n_in = n
        included.append([True, ] * n_in + [False, ] * (n - n_in))

    return included


def _include_snapshots(displacements, forces, energies, included):
    _forces = [forces[i][included_batch]
               for i, included_batch in enumerate(included)]
    _displacements = [np.array(displacements[i])[included_batch]
                      for i, included_batch in enumerate(included)]
    if energies:
        _energies = [energies[i][included_batch]
                     for i, included_batch in enumerate(included)]
    else:
        _energies = []
    return _displacements, _forces, _energies


def _remove_high_energy_snapshots(energies, included, ratio):
    """Reject high energy snapshots

    Parameters
    ----------
    energies : list of ndarray
        List of supercell total energies in each batch
    included : list of list of bool
        List of list of True/False as included snapshots in each batch
        Rejected elements are turned to False.
    ratio : float
        How much ratio of lowest energy snapshots is included after
        sorting by energy.

    Returns
    -------
    ret_included :
        Rejected snapshots are turned to False from 'included'.

    """

    concat_included = np.concatenate(included)
    concat_energies = np.concatenate(energies)
    included_energies = concat_energies[concat_included]
    included_indices = np.arange(len(concat_included))[concat_included]

    num_include = int(ratio * len(included_energies) + 0.5)
    if len(included_energies) < num_include:
        num_include = len(included_energies)
    _indices = np.argsort(included_energies)[:num_include]
    included_indices_after_energy = included_indices[_indices]

    bool_list = [False, ] * len(concat_included)
    for i in included_indices_after_energy:
        bool_list[i] = True
    ret_included = []
    count = 0
    for included_batch in included:
        ret_included.append(bool_list[count:(count + len(included_batch))])
        count += len(included_batch)
    return ret_included


def _modify_force_constants(ph):
    """Apply treatment to imaginary modes

    This method modifies force constants to make phonon frequencies
    be real from imaginary. This treatment is expected to be finally
    forgotten after many iterations. Therefore it is unnecessary
    to be physical and can be physically dirty. If it works, it is OK,
    though good treatment may contribute to quick convergence.

    1) All frequencies at commensurate points are converted to their
       absolute values. freqs -> |freqs|.
    2) Phonon frequencies in the interval 0.01 < |freqs| < 0.5
       are shifted by |freqs| + 1.

    """

    # temperature=300 is used just to invoke this feature.
    ph.generate_displacements(number_of_snapshots=1, temperature=300)
    rd = ph.random_displacements
    freqs = np.abs(rd.frequencies)
    condition = np.logical_and(freqs > 0.01, freqs < 0.5)
    rd.frequencies = np.where(condition, freqs + 1, freqs)
    rd.run_d2f()
    ph.force_constants = rd.force_constants


def _generate_random_displacements(ph,
                                   number_of_snapshots,
                                   temperature,
                                   random_seed=None):
    # Treatment of imaginary modes
    _modify_force_constants(ph)

    ph.generate_displacements(
        number_of_snapshots=number_of_snapshots,
        random_seed=random_seed,
        temperature=temperature)

    return ph.dataset


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
    temperature : Float
        Temperature (K).
    include_ratio : Float
        How much supercell forces are included from lowest supercell energies.
        Default is 1.0.
    lienar_decay : Bool
        This controls weights of previous phonon calculations. One previous
        phonon calculation is included 100%. (number_of_steps_for_fitting + 1)
        previous phonon calculation is not include. Between them, those are
        included with linear scalings. The snapshots in each phonon calculation
        are included from the first element of the list to the specified
        number, i.e., [:num_included_snapshots]. Default is False.
    random_seed : Int, optional
        Random seed used to sample in canonical ensemble harmonic oscillator
        space. The value must be 32bit unsigned int. Unless specified,
        random seed will not be fixed.
    initial_nodes : Dict, optional
        This gives the initial nodes that contain sets of forces, which are
        provided by PKs or UUIDs.

    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(PhonopyWorkChain,
                           exclude=['immigrant_calculation_folders',
                                    'calculation_nodes', 'dry_run'])
        spec.input('max_iteration', valid_type=Int, default=lambda: Int(10))
        spec.input('number_of_snapshots', valid_type=Int,
                   default=lambda: Int(100))
        spec.input('number_of_steps_for_fitting', valid_type=Int,
                   default=lambda: Int(4))
        spec.input('temperature', valid_type=Float,
                   default=lambda: Float(300.0))
        spec.input('include_ratio', valid_type=Float, default=lambda: Float(1))
        spec.input('linear_decay', valid_type=Bool,
                   default=lambda: Bool(False))
        spec.input('random_seed', valid_type=Int, required=False)
        spec.input('initial_nodes', valid_type=Dict, required=False)
        spec.outline(
            cls.initialize,
            if_(cls.import_initial_nodes)(
                cls.set_initial_nodes,
            ).else_(
                cls.run_initial_phonon,
            ),
            while_(cls.is_loop_finished)(
                cls.collect_displacements_and_forces,
                if_(cls.remote_phonopy)(
                    cls.run_force_constants_calculation_remote,
                    cls.generate_displacements,
                ).else_(
                    cls.generate_displacements_local,
                ),
                cls.run_phonon,
            ),
            cls.finalize,
        )

    def remote_phonopy(self):
        return self.inputs.remote_phonopy

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

    def run_initial_phonon(self):
        self.report("run_initial_phonon")
        inputs = self._get_phonopy_inputs(None, None)
        inputs['metadata'].label = "Initial phonon calculation"
        inputs['metadata'].description = "Initial phonon calculation"
        future = self.submit(PhonopyWorkChain, **inputs)
        self.ctx.initial_node = future
        self.report('{} pk = {}'.format(inputs['metadata'].label, future.pk))
        self.to_context(**{'phonon_0': future})

    def collect_displacements_and_forces(self):
        self.report("run_generate_displacements_%d" % self.ctx.iteration)

        num_batches = self.inputs.number_of_steps_for_fitting.value

        # Initial nodes are not specified, 0K phonon is in
        # self.ctx.initial_node. This is only once used to generate random
        # displacements. In the following steps, phonons calculated at
        # specified temperature are used to generate random displacements.
        if len(self.ctx.prev_nodes) == 0:
            nodes = [self.ctx.initial_node, ]
        elif len(self.ctx.prev_nodes) < num_batches:
            nodes = self.ctx.prev_nodes
        else:
            nodes = self.ctx.prev_nodes[-num_batches:]
        data_for_fc = {}
        for i, node in enumerate(nodes):
            data_for_fc['forces_%d' % (i + 1)] = node.outputs.force_sets
            ph_info = node.outputs.phonon_setting_info
            data_for_fc['phonon_setting_info_%d' % (i + 1)] = ph_info
        ret_Dict = collect_dataset(self.inputs.number_of_steps_for_fitting,
                                   self.inputs.include_ratio,
                                   self.inputs.linear_decay,
                                   **data_for_fc)
        self.ctx.dataset_for_fc = ret_Dict['dataset']

    def generate_displacements_local(self):
        if 'random_seed' in self.inputs:
            data_rd = {'random_seed': self.inputs.random_seed}
        else:
            data_rd = {}
        dataset = get_random_displacements(
            self.inputs.structure,
            self.inputs.number_of_snapshots,
            self.inputs.temperature,
            self.inputs.phonon_settings,
            self.ctx.dataset_for_fc,
            **data_rd)
        self.ctx.dataset = dataset

    def run_force_constants_calculation_remote(self):
        """Run force constants calculation by PhonopyCalculation"""

        self.report('remote force constants calculation %d' %
                    self.ctx.iteration)

        code_string = self.inputs.code_string.value
        builder = Code.get_from_string(code_string).get_builder()
        builder.structure = self.inputs.structure
        builder.settings = self.inputs.phonon_settings
        builder.metadata.options.update(self.inputs.options)
        builder.metadata.label = ("Force constants calculation %d" %
                                  self.ctx.iteration)
        builder.dataset = self.ctx.dataset_for_fc
        builder.fc_only = Bool(True)
        future = self.submit(builder)

        self.report('Force constants remote calculation: {}'.format(future.pk))
        label = 'force_constants_%d' % self.ctx.iteration
        self.to_context(**{label: future})

    def generate_displacements(self):
        label = 'force_constants_%d' % self.ctx.iteration
        fc_array = self.ctx[label].outputs.force_constants
        fc = fc_array.get_array('force_constants')
        phonon_setting_info = self.inputs.phonon_settings
        smat = phonon_setting_info['supercell_matrix']
        ph = Phonopy(phonopy_atoms_from_structure(self.inputs.structure),
                     supercell_matrix=smat,
                     primitive_matrix='auto')
        ph.force_constants = fc

        if 'random_seed' in self.inputs:
            random_seed = self.inputs.random_seed.value
        else:
            random_seed = None
        dataset = _generate_random_displacements(
            ph,
            self.inputs.number_of_snapshots.value,
            self.inputs.temperature.value,
            random_seed=random_seed)
        self.ctx.dataset = Dict(dict=dataset)

    def run_phonon(self):
        self.report("run_phonon_%d" % self.ctx.iteration)

        inputs = self._get_phonopy_inputs(self.ctx.dataset, False)
        label = "Phonon calculation %d" % self.ctx.iteration
        inputs['metadata'].label = label
        inputs['metadata'].description = label
        future = self.submit(PhonopyWorkChain, **inputs)
        self.ctx.prev_nodes.append(future)
        self.report('{} pk = {}'.format(label, future.pk))
        label = "phonon_%d" % self.ctx.iteration
        self.to_context(**{label: future})

    def finalize(self):
        self.report("IterHarmonicApprox finished at %d" %
                    (self.ctx.iteration - 1))

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

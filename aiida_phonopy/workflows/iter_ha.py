import numpy as np
import phonopy
from phonopy.harmonic.displacement import get_displacements_and_forces
from aiida.engine import WorkChain
from aiida.plugins import WorkflowFactory, DataFactory
from aiida.orm import Float, Int, QueryBuilder, Group
from aiida.engine import while_, calcfunction
from aiida_phonopy.common.utils import phonopy_atoms_from_structure

Dict = DataFactory('dict')
PhonopyWorkChain = WorkflowFactory('phonopy.phonopy')


@calcfunction
def get_random_displacements(structure,
                             number_of_snapshots,
                             temperature,
                             random_seed=None,
                             **data):
    displacements = []
    forces = []

    for i in range(len(data) // 2):
        forces.append(data['forces_%d' % (i + 1)].get_array('force_sets'))
        phonon_setting_info = data['ph_info_%d' % (i + 1)]
        dataset = phonon_setting_info['displacement_dataset']
        disps, _ = get_displacements_and_forces(dataset)
        displacements.append(disps)
    d = np.concatenate(displacements, axis=0)
    f = np.concatenate(forces, axis=0)

    phonon_setting_info = data['ph_info_1']
    smat = phonon_setting_info['supercell_matrix']
    ph = phonopy.load(unitcell=phonopy_atoms_from_structure(structure),
                      supercell_matrix=smat,
                      primitive_matrix='auto')
    ph.dataset = {'displacements': d, 'forces': f}
    ph.produce_force_constants(fc_calculator='alm')

    _modify_force_constants(ph)

    if random_seed is None:
        _random_seed = None
    else:
        _random_seed = random_seed.value

    ph.generate_displacements(
        number_of_snapshots=number_of_snapshots.value,
        random_seed=_random_seed,
        temperature=temperature.value)

    return Dict(dict=ph.dataset)


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

    """

    @classmethod
    def define(cls, spec):
        super(IterHarmonicApprox, cls).define(spec)
        spec.expose_inputs(PhonopyWorkChain,
                           exclude=['immigrant_calculation_folders',
                                    'calculation_nodes', 'dry_run'])
        spec.input('max_iteration',
                   valid_type=Int, required=False, default=Int(10))
        spec.input('number_of_snapshots',
                   valid_type=Int, required=False, default=Int(100))
        spec.input('number_of_steps_for_fitting',
                   valid_type=Int, required=False, default=Int(4))
        spec.input('random_seed', valid_type=Int, required=False)
        spec.input('temperature',
                   valid_type=Float, required=False, default=Float(300.0))
        spec.outline(
            cls.initialize,
            cls.run_initial_phonon,
            while_(cls.is_loop_finished)(
                cls.run_phonon,
            ),
        )

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
            displacements = get_random_displacements(
                nodes[-1].inputs.structure,
                self.inputs.number_of_snapshots,
                self.inputs.temperature,
                random_seed=self.inputs.random_seed,
                **data)
        else:
            displacements = get_random_displacements(
                nodes[-1].inputs.structure,
                self.inputs.number_of_snapshots,
                self.inputs.temperature,
                **data)
        return displacements

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

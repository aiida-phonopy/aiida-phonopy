import numpy as np
import phonopy
from phonopy.harmonic.displacement import get_displacements_and_forces
from aiida.engine import WorkChain
from aiida.plugins import WorkflowFactory, DataFactory
from aiida.orm import Float, Bool, Str, Int
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
    @classmethod
    def define(cls, spec):
        super(IterHarmonicApprox, cls).define(spec)
        spec.expose_inputs(PhonopyWorkChain,
                           exclude=['immigrant_calculation_folders',
                                    'calculation_nodes', 'dry_run'])
        spec.input('max_iteration',
                   valid_type=Int, required=False, default=Int(5))
        spec.input('number_of_snapshots',
                   valid_type=Int, required=False, default=Int(100))
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
        self.report("initialize")
        self.ctx.iteration = 0
        self.ctx.max_iteration = self.inputs.max_iteration.value
        self.ctx.prev_nodes = []
        self.ctx.num_nodes_for_average = 4
        if 'random_seed' in self.inputs:
            self.ctx.random_seed = self.inputs.random_seed.value
        else:
            self.ctx.random_seed = None

    def is_loop_finished(self):
        self.ctx.iteration += 1
        return self.ctx.iteration <= self.ctx.max_iteration

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

        n_ave = self.ctx.num_nodes_for_average
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

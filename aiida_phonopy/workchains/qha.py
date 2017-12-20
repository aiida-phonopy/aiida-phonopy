from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction
from aiida.work.run import run, submit, async

from aiida.orm import load_node, DataFactory, WorkflowFactory
from aiida.orm.data.base import Str, Float, Bool, Int

import numpy as np

__testing__ = False

ForceConstantsData = DataFactory('phonopy.force_constants')
ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

BandStructureData = DataFactory('phonopy.band_structure')
GruneisenPhonopy = WorkflowFactory('phonopy.gruneisen')
PhononPhonopy = WorkflowFactory('phonopy.phonon')


@workfunction
def phonopy_qha(**kwargs):
    from aiida_phonopy.workchains.gruneisen import get_qha

    n_structures = len([key for key, value in kwargs.items() if 'structure_' in key.lower()])

    volumes = []
    energies = []
    stresses = []
    entropy = []
    free_energy = []
    cv = []
    temperature = None
    for i in range(n_structures):
        volumes.append(kwargs.pop('structure_{}'.format(i)).get_cell_volume())
        output_data = kwargs.pop('output_data_{}'.format(i))
        energies.append(output_data.dict.energy)
        stresses.append(np.average(np.diag(output_data.dict.stress)))

        thermal_properties = kwargs.pop('thermal_properties_{}'.format(i))
        entropy.append(thermal_properties.get_array('entropy'))
        free_energy.append(thermal_properties.get_array('free_energy'))
        cv.append(thermal_properties.get_array('heat_capacity'))
        temperature = thermal_properties.get_array('temperature')

    free_energy = np.array(free_energy)
    entropy = np.array(entropy).T
    cv = np.array(cv).T
    volumes = np.array(volumes)
    energies = np.array(energies)
    stresses = np.array(stresses)

    eos = ArrayData()
    eos.set_array('volumes', volumes)
    eos.set_array('energies', energies)
    eos.set_array('stresses', stresses)

    qha_output = get_qha(eos,
                         temperature,
                         free_energy,
                         cv,
                         entropy)

    qha_results = ArrayData()
    for key, value in qha_output.items():
        qha_results.set_array(key, np.array(value))

    return {'qha_results': qha_results}


class QHAPhonopy(WorkChain):
    """
    Workchain to calculate the mode Gruneisen parameters
    """
    @classmethod
    def define(cls, spec):
        super(QHAPhonopy, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("ph_settings", valid_type=ParameterData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional arguments
        spec.input("num_expansions", valid_type=Int, required=False, default=Int(10))
        spec.input("use_nac", valid_type=Bool, required=False, default=Bool(True))

        spec.outline(cls.get_gruneisen_prediction, cls.create_unit_cell_expansions, cls.calculate_qha)

    def get_gruneisen_prediction(self):
        print('start qha (pk={})'.format(self.pid))


        if __testing__:
            self.ctx._content['gruneisen'] = load_node(19945)
            return

        future = submit(GruneisenPhonopy,
                        structure=self.inputs.structure,
                        ph_settings=self.inputs.ph_settings,
                        es_settings=self.inputs.es_settings,
                        pressure=Float(0.0),
                        use_nac=self.inputs.use_nac
                        )

        print ('gruneisen workchain: {}'.format(future.pid))
        return ToContext(gruneisen=future)

    def create_unit_cell_expansions(self):

        print('start Gruneisen (pk={})'.format(self.pid))
        prediction = self.ctx.gruneisen.out.prediction
        stress_range = prediction.dict.stress_range
        stress_delta = stress_range[-1] - stress_range[0]

        stress_samples = np.linspace(stress_range[0] - stress_delta * 0.5,
                                     stress_range[-1] + stress_delta * 0.5,
                                     self.inputs.num_expansions)

        print prediction.dict.stress_range
        print prediction.dict.volume_range

        # For testing
        if __testing__:
            self.ctx._content['phonon_0'] = load_node(19218)
            self.ctx._content['phonon_1'] = load_node(19221)
            self.ctx._content['phonon_2'] = load_node(19224)
            self.ctx._content['phonon_3'] = load_node(19227)
            self.ctx._content['phonon_4'] = load_node(19230)
            self.ctx._content['phonon_5'] = load_node(19233)
            self.ctx._content['phonon_6'] = load_node(19236)
            self.ctx._content['phonon_7'] = load_node(19239)
            self.ctx._content['phonon_8'] = load_node(19242)
            self.ctx._content['phonon_9'] = load_node(19245)
            return

        calcs = {}
        for i, stress in enumerate(stress_samples):

            future = submit(PhononPhonopy,
                            structure=self.inputs.structure,
                            ph_settings=self.inputs.ph_settings,
                            es_settings=self.inputs.es_settings,
                            pressure=Float(stress),
                            optimize=Bool(True),
                            use_nac=self.inputs.use_nac
                            )

            calcs['phonon_{}'.format(i)] = future
            print ('phonon workchain: {} {}'.format('phonon_{}'.format(i), future.pid))

        return ToContext(**calcs)

    def calculate_qha(self):
        print('start Gruneisen (pk={})'.format(self.pid))
        input_qha = {}
        for i in range(int(self.inputs.num_expansions)):
            input_qha['structure_{}'.format(i)] = self.ctx.get('phonon_{}'.format(i)).out.final_structure
            input_qha['output_data_{}'.format(i)] = self.ctx.get('phonon_{}'.format(i)).out.optimized_data
            input_qha['thermal_properties_{}'.format(i)] = self.ctx.get('phonon_{}'.format(i)).out.thermal_properties

        qha_results = phonopy_qha(**input_qha)

        for key, value in qha_results.items():
            self.out(key, value)

        print('QHA WorkChain finished')

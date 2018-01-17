from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction
from aiida.work.run import run, submit, async

from aiida.orm import load_node, DataFactory, WorkflowFactory, CalculationFactory, Code
from aiida.orm.data.base import Str, Float, Bool, Int
from aiida.work.workchain import _If, _While

import numpy as np

__testing__ = True

ForceConstantsData = DataFactory('phonopy.force_constants')
ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

PhononPhonopy = WorkflowFactory('phonopy.phonon')
PhononPhono3py = WorkflowFactory('phonopy.phonon3')


def generate_phono3py_params(structure, ph_settings, force_sets, nac_data=None):
    """
    Generate inputs parameters needed to do a remote phonopy calculation

    :param code: Code object of phono3py
    :param structure: StructureData Object that constains the crystal structure unit cell
    :param ph_settings: ParametersData object containing a dictionary with the phonopy input data
    :param force_sets: ForceSetsData object containing the atomic forces and displacement information
    :param force_sets: NacData object containing the dielectric tensor and Born effective charges info
    :return: Calculation process object, input dictionary
    """


    try:
        code = Code.get_from_string(ph_settings.dict.code['fc3'])
    except :
        code = Code.get_from_string(ph_settings.dict.code)

    plugin = code.get_attr('input_plugin')
    PhonopyCalculation = CalculationFactory(plugin)

    # The inputs
    inputs = PhonopyCalculation.process().get_inputs_template()

    # code
    inputs.code = code

    # structure
    inputs.structure = structure

    # parameters
    inputs.parameters = ph_settings

    # resources
    inputs._options.resources = ph_settings.dict.machine['resources']
    inputs._options.max_wallclock_seconds = ph_settings.dict.machine['max_wallclock_seconds']

    # data_sets
    inputs.data_sets = force_sets

    # non-analytical corrections
    if nac_data is not None:
        inputs.nac_data = nac_data

    return PhonopyCalculation.process(), inputs


class ThermalPhono3py(WorkChain):

    @classmethod
    def define(cls, spec):
        super(ThermalPhono3py, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("ph_settings", valid_type=ParameterData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional arguments
        spec.input("optimize", valid_type=Bool, required=False, default=Bool(True))
        spec.input("pressure", valid_type=Float, required=False, default=Float(0.0))
        spec.input("use_nac", valid_type=Bool, required=False, default=Bool(False))  # false by default
        spec.input("chunks", valid_type=Int, required=False, default=Int(100))
        spec.input("step", valid_type=Float, required=False, default=Float(1.0))
        spec.input("initial_cutoff", valid_type=Float, required=False, default=Float(2.0))

        spec.outline(cls.harmonic_calculation,
                     _While(cls.not_converged)(cls.get_data_sets,
                                               cls.get_thermal_conductivity),
                     cls.collect_data)

    def not_converged(self):

        is_converged = False

        if 'thermal_conductivity' in self.ctx:
            print ('thermal test')
            self.ctx.iteration += 1
            self.ctx.cutoff += float(self.inputs.step)
            self.ctx.conductivity = self.ctx.thermal_conductivity.out.kappa.get_array('kappa')
            self.ctx.input_data_sets = self.ctx.anharmonic.out.force_sets
            print ('kappa=', self.ctx.thermal_conductivity.out.kappa)
            if 'old_conductivity' in self.ctx:
                print ('old_conductivity')
                is_converged = np.allclose(self.ctx.conductivity, self.ctx.old_conductivity, rtol=0.3, atol=0.1)

            self.ctx.old_conductivity = self.ctx.conductivity

        else:
            print ('initial')
            self.ctx.cutoff = float(self.inputs.initial_cutoff)
            self.ctx.iteration = 1
            self.ctx.input_data_sets = self.ctx.harmonic.out.force_sets
            self.ctx.final_structure = self.ctx.harmonic.out.final_structure
            self.out('final_structure', self.ctx.final_structure)

            if not self.ctx.harmonic.out.dos.is_stable(tol=1e-3):
                self.report('This structure is not dynamically stable. Workchain will stop')
                exit()

        print ('iteration', self.ctx.iteration)
        print ('is_converged', is_converged)
        return not is_converged

    def harmonic_calculation(self):

        print('start thermal conductivity (pk={})'.format(self.pid))

        if __testing__:
            self.ctx._content['harmonic'] = load_node(83559)
            return

        self.report('Harmonic calculation')

        future = submit(PhononPhonopy,
                        structure=self.inputs.structure,
                        ph_settings=self.inputs.ph_settings,
                        es_settings=self.inputs.es_settings,
                        pressure=self.inputs.pressure,
                        optimize=self.inputs.optimize,
                        use_nac=self.inputs.use_nac
                        )
        return ToContext(harmonic=future)

    def get_data_sets(self):

        print ('cutoff:', self.ctx.cutoff)

        future = submit(PhononPhono3py,
                        structure=self.ctx.final_structure,
                        ph_settings=self.inputs.ph_settings,
                        es_settings=self.inputs.es_settings,
                        optimize=Bool(False),
                        use_nac=Bool(False),
                        cutoff=Float(self.ctx.cutoff),
                        chunks=self.inputs.chunks,
                        data_sets=self.ctx.input_data_sets,
                        #data_sets=load_node(81481),
                        )

        print('start phonon3 (pk={})'.format(self.pid))
        return ToContext(anharmonic=future)

    def get_thermal_conductivity(self):

        JobCalculation, calculation_input = generate_phono3py_params(structure=self.ctx.final_structure,
                                                                     ph_settings=self.inputs.ph_settings,
                                                                     force_sets=self.ctx.anharmonic.out.force_sets)
        future = submit(JobCalculation, **calculation_input)
        print('start phono3py (pk={})'.format(self.pid))

        return ToContext(thermal_conductivity=future)

    def collect_data(self):

        self.out('kappa', self.ctx.thermal_conductivity.out.kappa)
        self.report('finish thermal conductivity')

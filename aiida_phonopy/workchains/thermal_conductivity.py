from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.run import run, submit

from aiida.orm import load_node, DataFactory, WorkflowFactory, CalculationFactory, Code
from aiida.orm.data.base import Str, Float, Bool, Int
from aiida.work.workchain import _If, _While

from aiida_phonopy.workchains.phono3py_dist import generate_phono3py_params

import numpy as np
__testing__ = False

ForceConstantsData = DataFactory('phonopy.force_constants')
ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

PhononPhonopy = WorkflowFactory('phonopy.phonon')
PhononPhono3py = WorkflowFactory('phonopy.phonon3')
Phono3pyDist = WorkflowFactory('phonopy.phono3py_dist')


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
        spec.input("use_nac", valid_type=Bool, required=False, default=Bool(False))
        spec.input("chunks", valid_type=Int, required=False, default=Int(100))
        spec.input("step", valid_type=Float, required=False, default=Float(1.0))
        spec.input("initial_cutoff", valid_type=Float, required=False, default=Float(2.0))
        spec.input("gp_chunks", valid_type=Int, required=False, default=Int(1))
        # Convergence criteria
        spec.input("rtol", valid_type=Float, required=False, default=Float(0.3))
        spec.input("atol", valid_type=Float, required=False, default=Float(0.1))

        spec.outline(cls.harmonic_calculation,
                     _While(cls.not_converged)(cls.get_data_sets,
                                               cls.get_thermal_conductivity),
                     cls.collect_data)

    def not_converged(self):

        is_converged = False

        if 'thermal_conductivity' in self.ctx:
            self.ctx.iteration += 1
            self.ctx.cutoff += float(self.inputs.step)
            self.ctx.conductivity = self.ctx.thermal_conductivity.out.kappa.get_array('kappa')
            self.ctx.input_data_sets = self.ctx.anharmonic.out.force_sets

            print ('kappa=', self.ctx.thermal_conductivity.out.kappa)
            if 'prev_conductivity' in self.ctx:
                is_converged = np.allclose(self.ctx.conductivity, self.ctx.prev_conductivity,
                                           rtol=float(self.inputs.rtol),
                                           atol=float(self.inputs.atol))

            self.ctx.prev_conductivity = self.ctx.conductivity

        else:
            print ('initial')
            self.ctx.cutoff = float(self.inputs.initial_cutoff)
            self.ctx.iteration = 1
            self.ctx.input_data_sets = self.ctx.harmonic.out.force_sets
            self.ctx.final_structure = self.ctx.harmonic.out.final_structure
            # Add all harmonic data outputs
            self.out('final_structure', self.ctx.final_structure)
            self.out('thermal_properties', self.ctx.harmonic.out.thermal_properties)
            self.out('dos', self.ctx.harmonic.out.dos)
            self.out('band_structure', self.ctx.harmonic.out.band_structure)

            if not self.ctx.harmonic.out.dos.is_stable(tol=1e-3):
                self.report('This structure is not dynamically stable. WorkChain will stop')
                exit()

        print ('iteration {}'.format(self.ctx.iteration))
        print ('is_converged {}'.format(is_converged))
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

        print ('cutoff: {}'.format(self.ctx.cutoff))

        future = submit(PhononPhono3py,
                        structure=self.ctx.final_structure,
                        ph_settings=self.inputs.ph_settings,
                        es_settings=self.inputs.es_settings,
                        optimize=Bool(False),
                        cutoff=Float(self.ctx.cutoff),
                        chunks=self.inputs.chunks,
                        data_sets=self.ctx.input_data_sets,
                        # data_sets=load_node(81481), #  Test purposes only
                        )

        print('start phonon3 (pk={})'.format(self.pid))
        return ToContext(anharmonic=future)

    def get_thermal_conductivity(self):

        inputs = {'structure': self.ctx.final_structure,
                  'parameters': self.inputs.ph_settings,
                  'force_sets': self.ctx.anharmonic.out.force_sets}

        if bool(self.inputs.use_nac):
            inputs.update({'nac_data': self.ctx.harmonic.out.nac_data})

        if int(self.inputs.gp_chunks) > 1:
            inputs.update({'gp_chunks' : self.inputs.gp_chunks})
            future = submit(Phono3pyDist, **inputs)
        else:
            JobCalculation, calculation_input = generate_phono3py_params(**inputs)
            future = submit(JobCalculation, **calculation_input)
            print ('phono3py (pk = {})'.format(future.pid))

        return ToContext(thermal_conductivity=future)

    def collect_data(self):

        self.out('kappa', self.ctx.thermal_conductivity.out.kappa)
        self.report('finish thermal conductivity')

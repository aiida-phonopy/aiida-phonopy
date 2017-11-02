# Works run by the daemon (using submit)

from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext

from aiida.orm.data.base import Str, Float, Bool, Int
from aiida.work.workchain import _If, _While

import numpy as np
from generate_inputs import *
from parse_interface import *

class OptimizeStructure(WorkChain):
    """
    Workchain to do crystal structure optimization and ensure proper convergence
    """

    @classmethod
    def define(cls, spec):
        super(OptimizeStructure, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("machine", valid_type=ParameterData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional
        spec.input("pressure", valid_type=Float, required=False, default=Float(0.0))
        spec.input("tolerance_forces", valid_type=Float, required=False, default=Float(1e-5))
        spec.input("tolerance_stress", valid_type=Float, required=False, default=Float(1e-2))
        spec.input("max_iterations", valid_type=Int, required=False, default=Int(10))

        spec.outline(cls.optimize_cycle, _While(cls.not_converged)(cls.optimize_cycle), cls.get_data)

    def not_converged(self):

        print ('Check convergence')

        #output_array = self.ctx.get('optimize').out.output_array
        #forces = output_array.get_array('forces')
        #stresses = output_array.get_array('stress')

        parsed_data = parse_optimize_calculation(self.ctx.get('optimize'))
        forces = parsed_data['forces']
        stresses = parsed_data['stresses']

        if len(stresses.shape) > 2:  # For quantum espresso
            stresses = stresses[-1] * 10

        not_converged_forces = len(np.where(abs(forces) > float(self.inputs.tolerance_forces))[0])

        stress_compare_matrix = stresses - np.diag([float(self.inputs.pressure)]*3)
        not_converged_stress = len(np.where(abs(stress_compare_matrix) > float(self.inputs.tolerance_stress))[0])

        not_converged = not_converged_forces + not_converged_stress

        if not_converged == 0:
            print ('Converged')
            self.report('converged')
            return False

        self.report('Not converged: F:{} S:{}'.format(not_converged_forces, not_converged_stress))

        # Check max optimizations
        if self.ctx.counter > self.inputs.max_iterations:
            self.report('Max optimization iterations reached')
            return False

        return True

    def optimize_cycle(self):

        if not 'counter' in self.ctx:
            self.ctx.counter = 0

        self.ctx.counter +=1

        # self.ctx._get_dict()
        print ('start optimization')

        if not 'optimize' in self.ctx:
            structure = self.inputs.structure
        else:
            structure = self.ctx.optimize.out.output_structure

        JobCalculation, calculation_input = generate_inputs(structure,
                                                            self.inputs.machine,
                                                            self.inputs.es_settings,
                                                            pressure=self.inputs.pressure,
                                                            type='optimize',
                                                            )

        # calculation_input._label = 'optimize'
        future = submit(JobCalculation, **calculation_input)
        self.report('optimize calculation pk = {}'.format(future.pid))

        return ToContext(optimize=future)

    def get_data(self):
        print ('get_job')

        # self.ctx.structure = self.ctx.get('optimize').out.output_structure

        self.out('optimized_structure', self.ctx.optimize.out.output_structure)
        self.out('optimized_structure_data', self.ctx.optimize.out.output_parameters)

        return


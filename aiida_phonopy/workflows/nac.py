from aiida.engine import WorkChain
from aiida.plugins import DataFactory
from aiida_phonopy.common.builders import get_calcjob_builder, get_calcjob_inputs

Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')


class NacParamsWorkChain(WorkChain):
    """Wrapper to computer non-analytical term correction parameters"""
    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.input('structure', valid_type=StructureData, required=True)
        spec.input('calculator_settings', valid_type=Dict, required=True)

        spec.outline(
            cls.run_calculation,
            cls.finalize
        )

        spec.output('dielectrics', valid_type=ArrayData, required=True)
        spec.output('born_charges', valid_type=ArrayData, required=True)

    def run_calculation(self):
        """Born charges and dielectric constant calculation"""
        self.report('calculate born charges and dielectric constant')
        builder_inputs = get_calcjob_inputs(
            self.inputs.calculator_settings, self.inputs.structure)
        builder = get_calcjob_builder(
            self.inputs.structure,
            self.inputs.calculator_settings['code_string'],
            builder_inputs,
            label='born_and_epsilon')
        future = self.submit(builder)
        self.report('born_and_epsilon: {}'.format(future.pk))
        self.to_context(**{'calc': future})

    def finalize(self):
        self.out('dielectrics', self.ctx.calc.outputs.dielectrics)
        self.out('born_charges', self.ctx.calc.outputs.born_charges)

"""Workflow to calculate supercell forces."""

from aiida.engine import WorkChain
from aiida.plugins import DataFactory
from aiida_phonopy.common.builders import (
    get_calcjob_builder, get_calcjob_inputs)

Float = DataFactory('float')
Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')


class ForcesWorkChain(WorkChain):
    """Wrapper to compute supercell forces."""

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)
        spec.input('structure', valid_type=StructureData, required=True)
        spec.input('calculator_settings', valid_type=Dict, required=True)

        spec.outline(
            cls.run_calculation,
        )

        spec.output('forces', valid_type=ArrayData, required=True)
        spec.output('energy', valid_type=ArrayData, required=False)

        spec.exit_code(
            1001, 'ERROR_NO_FORCES',
            message='forces could not be retrieved from calculaton.')
        spec.exit_code(
            1002, 'ERROR_NO_ENERGY',
            message='energy could not be retrieved from calculaton.')

    def run_calculation(self):
        """Run supercell force calculation."""
        self.report('calculate supercell forces')
        builder_inputs = get_calcjob_inputs(
            self.inputs.calculator_settings, self.inputs.structure)
        builder = get_calcjob_builder(
            self.inputs.structure,
            self.inputs.calculator_settings['code_string'],
            builder_inputs,
            label=self.label)
        future = self.submit(builder)
        self.report('{} pk = {}'.format(self.label, future.pk))
        self.to_context(**{'calc': future})

    def finalize(self):
        """Finalize calculation."""
        outputs = self.ctx.calc.outputs
        if 'forces' in outputs and 'final' in outputs.forces.get_arraynames():
            forces = ArrayData()
            forces.set_array('forces', outputs.forces.get_array('final'))
            self.out('forces', forces)
        else:
            return self.exit_codes.ERROR_NO_FORCES

        ekey = 'energy_extrapolated'
        if 'energies' in outputs and ekey in outputs.energies.get_arraynames():
            energy = ArrayData()
            forces.set_array('energy', outputs.forces.get_array(ekey))
            self.out('energy', energy)
        else:
            return self.exit_codes.ERROR_NO_ENERGY

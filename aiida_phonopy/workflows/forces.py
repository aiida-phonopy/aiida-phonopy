"""Workflow to calculate supercell forces."""

import numpy as np
from aiida.engine import WorkChain, calcfunction
from aiida.plugins import DataFactory
from aiida.orm import Code
from aiida_phonopy.common.builders import (
    get_calcjob_inputs, get_calculator_process)


Float = DataFactory('float')
Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')


def _get_forces(outputs, code_string):
    """Return supercell force ArrayData."""
    code = Code.get_from_string(code_string)
    plugin_name = code.get_input_plugin_name()
    if plugin_name == 'vasp.vasp':
        if 'forces' in outputs and 'final' in outputs.forces.get_arraynames():
            forces_data = get_vasp_forces(outputs.forces)
        else:
            return None
    elif plugin_name == 'quantumespresso.pw':
        if ('output_trajectory' in outputs and
            'forces' in outputs.output_trajectory.get_arraynames()):
            forces_data = get_qe_forces(outputs.output_trajectory)
        else:
            return None
    return forces_data


@calcfunction
def get_vasp_forces(forces):
    """Return VASP forces ArrayData."""
    forces_data = ArrayData()
    forces_data.set_array('forces', forces.get_array('final'))
    forces_data.label = 'forces'
    return forces_data


@calcfunction
def get_qe_forces(output_trajectory):
    """Return QE forces ArrayData."""
    forces_data = ArrayData()
    forces_data.set_array('forces', output_trajectory.get_array('forces')[-1])
    forces_data.label = 'forces'
    return forces_data


def _get_energy(outputs, code_string):
    """Return supercell energy ArrayData."""
    code = Code.get_from_string(code_string)
    plugin_name = code.get_input_plugin_name()
    if plugin_name == 'vasp.vasp':
        ekey = 'energy_extrapolated'
        if 'energies' in outputs and ekey in outputs.energies.get_arraynames():
            energy_data = get_vasp_energy(outputs.energies)
        else:
            return None
    elif plugin_name == 'quantumespresso.pw':
        if ('output_parameters' in outputs and
            'energy' in outputs.output_parameters.dict):
            energy_data = get_qe_energy(outputs.output_parameters)
    return energy_data


@calcfunction
def get_vasp_energy(energies):
    """Return VASP energy ArrayData."""
    energy_data = ArrayData()
    ekey = 'energy_extrapolated'
    energy_data.set_array('energy', np.array(
        [energies.get_array(ekey), ], dtype=float))
    energy_data.label = 'energy'
    return energy_data


@calcfunction
def get_qe_energy(output_parameters):
    """Return VASP energy ArrayData."""
    energy_data = ArrayData()
    energy_data.set_array('energy', np.array(
        [output_parameters.dict.energy, ], dtype=float))
    energy_data.label = 'energy'
    return energy_data


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
            cls.finalize
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
        process_inputs = get_calcjob_inputs(self.inputs.calculator_settings,
                                            self.inputs.structure,
                                            label=self.metadata.label)
        CalculatorProcess = get_calculator_process(
            self.inputs.calculator_settings['code_string'])
        future = self.submit(CalculatorProcess, **process_inputs)
        self.report('{} pk = {}'.format(self.metadata.label, future.pk))
        self.to_context(**{'calc': future})

    def finalize(self):
        """Finalize force calculation."""
        outputs = self.ctx.calc.outputs
        self.report('create forces ArrayData')
        forces = _get_forces(
            outputs,
            self.inputs.calculator_settings['code_string'])
        if forces is None:
            return self.exit_codes.ERROR_NO_FORCES
        else:
            self.out('forces', forces)

        self.report('create energy ArrayData')
        energy = _get_energy(
            outputs,
            self.inputs.calculator_settings['code_string'])
        if energy is None:
            return self.exit_codes.ERROR_NO_ENERGY
        else:
            self.out('energy', energy)

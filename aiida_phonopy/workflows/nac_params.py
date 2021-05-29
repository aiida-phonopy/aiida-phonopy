"""Workflow to calculate NAC params."""

from aiida.engine import WorkChain, calcfunction
from aiida.plugins import DataFactory
from aiida_phonopy.common.builders import get_calcjob_builder, get_calcjob_inputs
from aiida_phonopy.common.utils import phonopy_atoms_from_structure
from phonopy.structure.symmetry import symmetrize_borns_and_epsilon

Float = DataFactory('float')
Dict = DataFactory('dict')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')


@calcfunction
def get_nac_params(born_charges, epsilon, nac_structure, symmetry_tolerance,
                   primitive=None):
    """Obtain Born effective charges and dielectric constants in primitive cell.

    When Born effective charges and dielectric constants are calculated within
    phonopy workchain, those values are calculated in the primitive cell.
    However using immigrant, the cell may not be primitive cell and can be
    unit cell. In this case, conversion of data is necessary. This conversion
    needs information of the structure where those values were calcualted and
    the target primitive cell structure.

    Two kargs parameters
    primitive : StructureData
    symmetry_tolerance : Float

    """
    borns = born_charges.get_array('born_charges')
    eps = epsilon.get_array('epsilon')

    nac_cell = phonopy_atoms_from_structure(nac_structure)
    kwargs = {}
    kwargs['symprec'] = symmetry_tolerance.value
    if primitive is not None:
        kwargs['primitive'] = phonopy_atoms_from_structure(primitive)
    borns_, epsilon_ = symmetrize_borns_and_epsilon(
        borns, eps, nac_cell, **kwargs)

    nac_params = ArrayData()
    nac_params.set_array('born_charges', borns_)
    nac_params.set_array('epsilon', epsilon_)
    nac_params.label = 'born_charges & epsilon'

    return nac_params


class NacParamsWorkChain(WorkChain):
    """Wrapper to compute non-analytical term correction parameters."""

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)
        spec.input('structure', valid_type=StructureData, required=True)
        spec.input('calculator_settings', valid_type=Dict, required=True)
        spec.input('symmetry_tolerance', valid_type=Float,
                   default=lambda: Float(1e-5))

        spec.outline(
            cls.run_calculation,
            cls.create_nac_params
        )

        spec.output('nac_params', valid_type=ArrayData, required=True)

        spec.exit_code(
            1001, 'ERROR_NO_BORN_EFFECTIVE_CHARGES',
            message=('Born effecti charges could not be retrieved '
                     'from calculaton.'))
        spec.exit_code(
            1002, 'ERROR_NO_DIELECTRIC_CONSTANT',
            message=('dielectric constant could not be retrieved '
                     'from calculaton.'))

    def run_calculation(self):
        """Born charges and dielectric constant calculation."""
        self.report('Calculate born charges and dielectric constant')
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

    def create_nac_params(self):
        """Create NAC params ArrayData."""
        self.report('Create nac params data')

        calc = self.ctx.calc
        if type(calc) is dict:
            calc_dict = calc
            structure = calc['structure']
        else:
            calc_dict = calc.outputs
            structure = calc.inputs.structure

        if 'born_charges' not in calc_dict:
            return self.exit_codes.ERROR_NO_BORN_EFFECTIVE_CHARGES

        if 'dielectrics' not in calc_dict:
            return self.exit_codes.ERROR_NO_DIELECTRIC_CONSTANT

        nac_params = get_nac_params(calc_dict['born_charges'],
                                    calc_dict['dielectrics'],
                                    structure,
                                    self.inputs.symmetry_tolerance)
        self.out('nac_params', nac_params)

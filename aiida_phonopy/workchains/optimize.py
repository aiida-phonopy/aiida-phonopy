from aiida import load_dbenv, is_dbenv_loaded
if not is_dbenv_loaded():
    load_dbenv()

from aiida.work.workchain import WorkChain, ToContext
from aiida.work.workfunction import workfunction
from aiida.work.run import run, submit, async

from aiida.orm.data.base import Str, Float, Bool, Int
from aiida.work.workchain import _If, _While

from aiida.orm import DataFactory, load_node

ForceConstantsData = DataFactory('phonopy.force_constants')
ParameterData = DataFactory('parameter')
ArrayData = DataFactory('array')
StructureData = DataFactory('structure')

from aiida_phonopy.common.generate_inputs import generate_inputs
from aiida_phonopy.common.parse_interface import parse_optimize_calculation

import numpy as np


@workfunction
def standardize_cell(structure):
    import spglib
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
    from phonopy.structure.atoms import atom_data

    bulk = PhonopyAtoms(symbols=[site.kind_name for site in structure.sites],
                        positions=[site.position for site in structure.sites],
                        cell=structure.cell)

    structure_data = (structure.cell,
                      bulk.get_scaled_positions(),
                      bulk.get_atomic_numbers())

    #lattice, refined_positions, numbers = spglib.refine_cell(structure_data, symprec=1e-5)
    lattice, standardized_positions, numbers = spglib.standardize_cell(structure_data,
                                                                       symprec=1e-5,
                                                                       to_primitive=False,
                                                                       no_idealize=False)

    symbols = [atom_data[i][1] for i in numbers]

    # print lattice, standardized_positions, numbers
    # print [site.kind_name for site in structure.sites]
    standardized_bulk = PhonopyAtoms(symbols=symbols,
                                     scaled_positions=standardized_positions,
                                     cell=lattice)

    # create new aiida structure object
    standarized = StructureData(cell=standardized_bulk.get_cell())
    for position, symbol in zip(standardized_bulk.get_positions(), standardized_bulk.get_chemical_symbols()):
        standarized.append_atom(position=position,
                                      symbols=symbol)

    return {'standardized_structure': standarized}


class OptimizeStructure(WorkChain):
    """
    Workchain to do crystal structure optimization and ensure proper convergence
    """

    @classmethod
    def define(cls, spec):
        super(OptimizeStructure, cls).define(spec)
        spec.input("structure", valid_type=StructureData)
        spec.input("es_settings", valid_type=ParameterData)
        # Optional
        spec.input("pressure", valid_type=Float, required=False, default=Float(0.0))
        spec.input("tolerance_forces", valid_type=Float, required=False, default=Float(1e-5))
        spec.input("tolerance_stress", valid_type=Float, required=False, default=Float(1e-2))
        spec.input("max_iterations", valid_type=Int, required=False, default=Int(3))
        spec.input("standarize_cell", valid_type=Bool, required=False, default=Bool(True))

        spec.outline(cls.optimize_cycle, _While(cls.not_converged)(cls.optimize_cycle), cls.get_data)

    def not_converged(self):

        print ('Check convergence')

        parsed_data = parse_optimize_calculation(self.ctx.optimize)
        forces = parsed_data['forces']
        stresses = parsed_data['stresses']

        # Temporal way around. This has to be uniform for all codes
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

        self.ctx.counter += 1

        # self.ctx._get_dict()
        print ('start optimization')

        if not 'optimize' in self.ctx:
            structure = self.inputs.structure
        else:
            structure = parse_optimize_calculation(self.ctx.optimize)['output_structure']

        if self.inputs.standarize_cell:
            structure = standardize_cell(structure)['standardized_structure']

        JobCalculation, calculation_input = generate_inputs(structure,
                                                            self.inputs.es_settings,
                                                            pressure=self.inputs.pressure,
                                                            type='optimize',
                                                            )

        calculation_input._label = 'optimize'
        future = submit(JobCalculation, **calculation_input)
        self.report('optimize calculation pk = {}'.format(future.pid))

        return ToContext(optimize=future)
        # self.ctx.optimize = load_node(5487)  # for testing purposes

    def get_data(self):

        self.out('optimized_structure', parse_optimize_calculation(self.ctx.optimize)['output_structure'])
        self.out('optimized_structure_data', self.ctx.optimize.out.output_array)

        return


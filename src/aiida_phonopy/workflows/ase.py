# -*- coding: utf-8 -*-
"""Workflow for automatic frozen phonons calculations using Phonopy and any ASE calculator."""
from typing import Callable, Dict, List, Union

from aiida import orm
from aiida.common import AttributeDict
from aiida.engine import calcfunction, if_, while_
from aiida_pythonjob import PythonJob
from aiida_pythonjob.data.atoms import AtomsData
from ase.atoms import Atoms

from .phonopy import PhonopyWorkChain


@calcfunction
def get_forces_array(result: orm.Data) -> orm.ArrayData:
    """Convert the forces from a PythonJob result to ArrayData."""
    array = orm.ArrayData()
    array.set_array('forces', result.value)
    return array


class PhonopyAseWorkChain(PhonopyWorkChain):
    """Workflow for automated frozen phonons calculations using Phonopy and any ASE calculator.

    Phonopy is used to produce structures with displacements, while the forces are calculated
    with any calculator that can be interfaced with ASE. The submission of the calculation
    using ASE are handled via aiida-pythonjob.
    """

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        # yapf: disable
        super().define(spec)

        spec.expose_inputs(PythonJob, namespace='pythonjob')
        spec.input('settings.max_concurrent_base_workchains', valid_type=int, non_db=True,  default=1,
            help='Maximum number of concurrent running force calculation. If negative, run all at once.')

        spec.outline(
            cls.setup,
            cls.run_supercells,
            while_(cls.should_run_forces)(
                cls.run_forces,
            ),
            cls.inspect_forces,
            cls.run_results,
            if_(cls.should_run_phonopy)(
              cls.run_phonopy,
              cls.inspect_phonopy,
            ),
        )
        # yapf: enable

    @classmethod
    def get_populated_builder(
        cls,
        structure: Union[orm.StructureData, Atoms],
        calculator: Callable,
        pythonjob_inputs: Dict,
        max_number_of_atoms: int = None,
        supercell_matrix: List[List[int]] = None,
        primitive_matrix: List[List[int]] = None,
        symmetry_tolerance: float = 1.0e-5,
        is_symmetry: bool = True,
        displacement_generator: Dict = None,
        phonopy_inputs: Dict = None,
    ):
        """Return a builder with all required inputs."""
        from aiida.tools.data.structure import structure_to_spglib_tuple
        from aiida_pythonjob import prepare_pythonjob_inputs
        import numpy as np
        from phonopy.structure.cells import estimate_supercell_matrix
        import spglib

        # We need to define the function here to avoid validation error.
        # If would import this function from somewhere else, aiida-pythonjob
        # would serialize this function with a module path which raises an error.
        def calculate_forces(calculator, atoms):
            """Calculate the forces of an ASE Atoms structure given an ASE calculator."""
            atoms.calc = calculator
            return atoms.get_forces()

        builder = cls.get_builder()

        if isinstance(structure, orm.StructureData):
            builder.structure = structure
        elif isinstance(structure, Atoms):
            builder.structure = orm.StructureData(ase=structure)
        else:
            raise ValueError('structure must be either StructureData or ASE Atoms')

        if supercell_matrix is None and max_number_of_atoms is None:
            raise ValueError('Either `supercell_matrix` or `max_number_of_atoms` must be provided')
        if supercell_matrix is not None and max_number_of_atoms is not None:
            raise ValueError('Only one between `supercell_matrix` and `max_number_of_atoms` must be provided, not both')

        if supercell_matrix is not None:
            builder.supercell_matrix = orm.List(supercell_matrix)
        if max_number_of_atoms is not None:
            spglib_cell, _, _ = structure_to_spglib_tuple(builder.structure)
            spglib_dataset = spglib.get_symmetry_dataset(spglib_cell)
            to_conv = spglib_dataset.transformation_matrix  # to conventional cell
            matrix = estimate_supercell_matrix(spglib_dataset, max_number_of_atoms)
            supercell_matrix = np.dot(np.linalg.inv(to_conv), np.diag(matrix)).tolist()
            builder.supercell_matrix = orm.List(supercell_matrix)

        builder.symmetry_tolerance = orm.Float(symmetry_tolerance)
        builder.is_symmetry = orm.Bool(is_symmetry)
        if primitive_matrix is not None:
            builder.primitive_matrix = orm.List(primitive_matrix)
        if displacement_generator is not None:
            builder.displacement_generator = orm.Dict(displacement_generator)

        pythonjob_dict_inputs = prepare_pythonjob_inputs(
            function=calculate_forces,
            function_inputs={
                'calculator': calculator,
                'atoms': builder.structure.get_ase()
            },
            **pythonjob_inputs
        )

        if phonopy_inputs is not None:
            for key, value in phonopy_inputs.items():
                builder['phonopy'][key] = value
        else:
            builder.pop('phonopy', None)

        return {**builder._data, 'pythonjob': pythonjob_dict_inputs}

    def should_run_forces(self):
        """Whether to run or not forces."""
        return len(self.ctx.supercells) > 0

    def run_forces(self):
        """Run the forces for each supercell."""
        n_base_parallel = self.inputs.settings.max_concurrent_base_workchains
        if self.inputs.settings.max_concurrent_base_workchains < 0:
            n_base_parallel = len(self.ctx.supercells)

        for _ in self.ctx.supercells[:n_base_parallel]:
            key, supercell = self.ctx.supercells.pop(0)
            num = key.split('_')[-1]
            label = f'forces_{num}'

            inputs = AttributeDict(self.exposed_inputs(PythonJob, namespace='pythonjob'))
            inputs.function_inputs.atoms = AtomsData(supercell.get_ase())

            inputs.metadata.label = label
            inputs.metadata.call_link_label = label

            future = self.submit(PythonJob, **inputs)
            self.report(f'submitting `PythonJob` <PK={future.pk}> with supercell n.o {num}')
            self.to_context(**{label: future})

    def inspect_forces(self):
        """Inspect the calculation of forces for each supercell."""
        failed_runs = []

        for label, calculation in self.ctx.items():
            if label.startswith('forces_'):
                if calculation.is_finished_ok:
                    forces = calculation.outputs.result
                    self.out(f"supercells_forces.forces_{label.split('_')[-1]}", get_forces_array(forces))
                else:
                    self.report(f'PythonJob with <PK={calculation.pk}> failed')
                    failed_runs.append(calculation.pk)

        if failed_runs:
            self.report('one or more workchains did not finish succesfully')
            return self.exit_codes.ERROR_SUB_PROCESS_FAILED

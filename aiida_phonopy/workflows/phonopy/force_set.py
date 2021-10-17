"""Turn-key solution for automatic frozen phonons calculations."""
from abc import ABCMeta, abstractmethod
from aiida.engine import WorkChain, if_
from aiida import orm
from aiida_phonopy.calculations.functions.preprocess import phonopy_preprocess  
from aiida_phonopy.common.utils import get_force_sets
from copy import deepcopy

def validate_inputs(value, _):
    """Validate the entire input namespace."""
    disp = value['displacements']
 
    # Check the inputs has just one between ``dataset`` and ``random``
    if "dataset" in disp and "random" in disp:
        return 'too many inputs in the ``displacements`` namespace.'

    # Validate the ``displacements.random`` namespace
    if "random" in disp:
        rand = disp["random"]
        if "number_of_snapshots" in rand or "random_seed" in rand:
            if not("number_of_snapshots" in rand and "random_seed" in rand):
                return 'incomplete ``displacements.random`` inputs.'

def validate_supercell_matrix(value, _):
    """Validate the `validate_supercell_matrix` input."""
    if not len(value) == 3:
        return 'need exactly 3 diagonal elements or 3x3 arrays.'
    
    for row in value:
        if not len(row) == 3:
            return 'supercell matrix need to have 3x1 or 3x3 shape.'

        for element in row:
            if not type(element)==int:
                return f'type `{type(element)}` of {element} is not an accepted type in supercell matrix; only `int` is valid.'

def validate_positive_integer(value, _):
    """Validate that `value` is positive."""
    if not value.value > 0:
        return f'{value} is not positive.'


class ForceSetWorkChain(WorkChain, metaclass=ABCMeta):
    """
    Workflow to compute automatically the force set of a given structure
    using the frozen phonons approach.
    
    Phonopy is used to produce structures with displacements,
    while the forces are calculated with a quantum engine of choice.
    """

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)
        spec.input("structure", valid_type=orm.StructureData, help='The structure at equilibrium volume.')
        spec.input("supercell_matrix", valid_type=orm.List, required=True, validator=validate_supercell_matrix,
            help='The matrix used to generate the supercell from the input structure in the List format. Allowed shapes are 3x1 and 3x3 lists.')
        spec.input_namespace('displacements',
            help='The inputs regarding the type and size of the finite displacements. ´´dataset´´ and ´´random´´ are optional.')
        spec.input("displacements.distance", valid_type=orm.Float, default=lambda: orm.Float(0.01),
            help='The size of the finite displacement in angstrom.') # in angstrom
        spec.input("displacements.dataset", valid_type=orm.Dict, required=False,
            help='Dataset of displacements in a readable phonopy format.')
        spec.input("displacements.random.random_seed", valid_type=orm.Int, required=False, validator=validate_positive_integer)
        spec.input("displacements.random.number_of_snapshots", valid_type=orm.Int, required=False, validator=validate_positive_integer)
        spec.input("symmetry_tolerance", valid_type=orm.Float, default=lambda: orm.Float(1e-5),
            help='Symmetry tolerance for space group analysis on the input structure.')
        spec.input("subtract_residual_forces", valid_type=orm.Bool, default=lambda: orm.Bool(False),
            help='If `True` it performs a pristine supercell calculation.') # to explain better?
        spec.inputs.validator = validate_inputs

        spec.outline(
            cls.setup,
            cls.run_forces,
            cls.inspect_forces,
            cls.run_results,
        )

        spec.output_namespace("supercells", valid_type=orm.StructureData,
            help='The (super)cell structures containing the displacements.')
        spec.output_namespace("primitive", dynamic=True,
            help='The outputs regarding the primitive cell.')
        spec.output("displacements_dataset", valid_type=orm.ArrayData, # namespace?
            help='The set of displacements referred to the (super)cell.')
        spec.output_namespace("supercells_forces", valid_type=orm.ArrayData,
            help='The forces acting on the atoms of each (super)cell.') 

        spec.output("force_sets", valid_type=orm.ArrayData, required=False,
            help='The sets of forces needed for the phonopy post-process.')

    def setup(self):
        """Setup the workflow generating (super)cells and primitive cell."""

        preprocess = phonopy_preprocess(
            self.inputs.structure, self.inputs.symmetry_tolerance,
            self.inputs.displacements, self.inputs.supercell_matrix)

        for key in preprocess.keys():
            self.out(key, preprocess['key'])

        self.ctx.supercells = deep.copy(preprocess["supercells"])
        if not self.inputs.subtract_residual_forces:
            self.ctx.supercells.pop("supercell", None)
            # ^ - need to change the force set calculation
    
    @abstractmethod
    def run_forces(self):
        """Run supercells using the selected quantum engine code."""
        pass

    @abstractmethod
    def inspect_forces(self):
        """
        Inspect all scfs and return an error message if 
        one or more did not finished correctly
        """
        pass

    def create_force_sets(self):
        """Build datasets from forces of supercells with displacments."""
        # should move into the phonopy calculation ?
        force_set = get_force_sets(**self.ctx.forces)
        self.out("force_set", force_set)

    @abstractmethod
    def run_results(self):
        """Collect and expose the outputs."""
        pass


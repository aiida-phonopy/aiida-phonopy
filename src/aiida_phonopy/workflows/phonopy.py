# -*- coding: utf-8 -*-
"""Abstract workflow for automatic frozen phonons calculations."""

from abc import ABCMeta

from aiida import orm
from aiida.engine import WorkChain

from aiida_phonopy.data import ForceConstantsData, PhonopyData, PreProcessData


def validate_matrix(value, _):
    """Validate the `supercell_matrix` and `primitive_matrix` inputs."""
    import numpy as np

    if not isinstance(value, (list, orm.List, np.ndarray)):
        return 'value is not of the right type; only `list`, `aiida.orm.List` and `numpy.ndarray`'

    if isinstance(value, np.ndarray):
        value = value.tolist()

    if not len(value) == 3:
        return 'need exactly 3 diagonal elements or 3x3 arrays.'

    for row in value:
        if isinstance(row, list):
            if not len(row) in [0, 3]:
                return 'matrix need to have 3x1 or 3x3 shape.'
            for element in row:
                if not isinstance(element, (int, float)):
                    return (
                        f'type `{type(element)}` of {element} is not an accepted '
                        'type in matrix; only `int` and `float` are valid.'
                    )


def validate_positive_integer(value, _):
    """Validate that `value` is positive."""
    if not value.value > 0:
        return f'{value} is not positive.'


def validate_nac(value, _):
    """Validate that `value` is positive."""
    try:
        value.get_array('dielectric')
        value.get_array('born_charges')
    except KeyError:
        return 'data does not contain `dieletric` and/or `born_charges` arraynames.'


def validate_inputs(inputs, _):
    """Validate the entire inputs namespace."""
    ensamble_inputs = ['structure', 'supercell_matrix', 'primitive_matrix', 'symprec', 'is_symmetry']

    given_inputs = []

    for input_ in ensamble_inputs:
        if input_ in inputs:
            given_inputs.append(input_)

    if 'preprocess_data' in inputs and given_inputs:
        return 'too many inputs have been provided.'

    if given_inputs and 'structure' not in given_inputs:
        return 'a structure data is required'

    if not given_inputs and not 'preprocess_data' in inputs:
        return 'at least one between `preprocess_data` and `structure` must be provided in input'


class PhonopyWorkChain(WorkChain, metaclass=ABCMeta):
    """
    Abstract Workflow to compute automatically the force sets and force constants
    of a given structure using the frozen phonons approach. Phonopy is used to produce
    structures with displacements, while the forces are calculated with a quantum engine of choice.

    This workchain is meant to be used as a base for other specific force calculato plugin workchains,
    or as an example on how to set a possible workchain/workflow. For this reason, the outline of
    this class is not defined, while it provides the inputs and a `setup` method, which can be used
    in a specific workflow outline. Ideally, the workflow would look like:

    1. Setup the preprocess data.

        This is already provided in this class. It setups a `PreProcessData` node, from where
        supercell, primitive cell and supercells with displacements can be easily extracted using
        the methods of the nodes. This node can be taken from `self.ctx.preprocess_data`, and used
        during the outline of the workflow.

    2. Run supercells using the selected quantum engine/force calculator code.

        In specific code implementations, a force calculation on supercells needs to be run.
        To get these supercells, one need simply to run:

        ```self.ctx.preprocess_data.calcfunctions.get_supercells_with_displacements()```

        This will return a dictionary with all the supercells as StructureData to run for the phonon calculation.
        The keys of this dictionary are of the type `supercell_{number}`, where `number` is an integer.
        These numbers are essentials since the `phonopy` force sets is generated following these numbers,
        in order to make sure to refer to the correct displacement. Thus, it is required to keep track
        of them.
        Moreover,a calculation over the pristine supercell structure should be run before hand as reference.
        This structure can instead be gotten via:

        ```self.ctx.preprocess_data.calcfunctions.get_supercell()```

        This will return a StructureData without any label.

        For an example of implementation, refer to aiidateam/aiida-common-worfklows.

        * Note: some type of force calculation needs to map some variables from the unitcell to the supercell
        (and in certain case even the primitive cell), e.g. the atomic spin in VASP. Since this is code dependent,
        you will need to map these parameters before launching the force calculation of a certain supercell
        with displacement. This information can be gotten via:

        ```self.ctx.preprocess_data.get_cells_mappings()```

        Moreover, consider that cells in phonopy will always (re)fold the atoms in order to have positive coordinates.

    3. Inspect all runs and expose the forces and energies (not mandatory) outputs.

        * Suggested: when the calculation on each supercell has finished (correctly)
        expose the output forces (and energies) in the dynamical `supercells_forces(energies)` namespace(s).
        Provide each supercell forces as an `ArrayData` with the forces stored as `forces`
        (e.g. if your code plugin stores  the forces in `TrajectoryData`, extract them with a `calcfunction`).
        Expose each `ArrayData` choosing a **common prefix**, while as **suffix use
        _{number}**, with `{number}` referring to the correspective supercell label suffix (that you are supposed to
        keep track somewhere, e.g. in the label of the code calculation/workchain).
        Now you can gather all the information in one data noe, i.e. in a `PhonopyData` node.
        To do so, you can simple run:

        ```self.ctx.preprocess_data.calcfunctions.generate_phonopy_data(**self.outputs.supercells_forces)```

        and then expose it as output in the `output_phonopy_data` namespace.

        * Alternatively: instead of exposing the supercell forces as outputs, you can directly gather all the forces
        in a dictionary and run directly to the `generate_phonopy_data` method using this dictionary (always using
        the double *).

        See the implementation in aiidateam/aiida-common-workflows for an example.

    4. (optional) Run the non-analytical constants on the primitive cell.

        Non-analytical constants should be run for polar insulators. These require usually a linear response code
        or a finite difference approach (e.g. using the electric enthalpy). Since this is usually the most expensive
        part, you should run them on the primitive cell. To get it, use:

        ```self.ctx.preprocess_data.calcfunctions.get_primitive_cell()```

        If you compute also these, collect the dielectric tensor and the effectic born charges in an ArrayData,
        with the arraynames `dielectric` and `born_charges` (in Cartesian coordinates!).
        Then, gather all the information of nac and forces in a unique `PhonopyData` via:

        ```
        self.ctx.preprocess_data.calcfunctions.generate_phonopy_data(
            nac_parameters=nac_paramters,
            **self.outputs.supercells_forces
            )
        ```

        and expose the output.

        * Note: we require in the input for generating the full phonopy data, to give the nac in the primitive cell.
        The primitive cell of phonopy will just rotate the lattice vectors, thus mantaining the Cartasian coordinate
        system. It can happen, though, that the unitcell is not the primitive cell of the system, meaning that the
        primitive cell will contain less atoms. We expect in input the nac computed on this number of atoms. If you
        want, for some reason, compute the nac on the unitcell, you will need to get the reduced nac.
        To do so, you can consider using a built-in function in phonopy, namely:

        ```phonopy.structure.symmetry.elaborate_borns_and_epsilon```

    """

    _ENABLED_DISPLACEMENT_GENERATOR_FLAGS = {
        'distance': [float],
        'is_plusminus': ['auto', float],
        'is_diagonal': [bool],
        'is_trigonal': [bool],
        'number_of_snapshots': [int, None],
        'random_seed': [int, None],
        'temperature': [float, None],
        'cutoff_frequency': [float, None],
    }

    _ENABLED_FC_OPTIONS_FLAGS = {
        'calculate_full_force_constants': [bool],
        'fc_calculator': [str, None],
        'fc_calculator_options': [None],
    }

    @classmethod
    def define(cls, spec):
        """Define inputs, outputs, and outline."""
        super().define(spec)

        spec.input(
            'preprocess_data',
            valid_type=(PhonopyData, PreProcessData),
            required=False,
            help='The preprocess data for frozen phonon calcualtion.'
        )
        spec.input(
            'structure', valid_type=orm.StructureData, required=False, help='The structure at equilibrium volume.'
        )
        spec.input(
            'supercell_matrix',
            valid_type=orm.List,
            required=False,
            validator=validate_matrix,
            help=(
                'The matrix used to generate the supercell from the input '
                'structure in the List format. Allowed shapes are 3x1 and 3x3 lists.'
            ),
        )
        spec.input(
            'primitive_matrix',
            valid_type=orm.List,
            required=False,
            validator=validate_matrix,
            help=(
                'The matrix used to generate the primitive cell from the input '
                'structure in the List format. Allowed shapes are 3x1 and 3x3 lists.'
            ),
        )
        spec.input(
            'symmetry_tolerance',
            valid_type=orm.Float,
            required=False,
            help='Symmetry tolerance for space group analysis on the input structure.',
        )
        spec.input(
            'is_symmetry',
            valid_type=orm.Bool,
            required=False,
            help='Whether using or not the space group symmetries.',
        )
        spec.input(
            'displacement_generator',
            valid_type=orm.Dict,
            required=False,
            validator=cls._validate_displacements,
            help=(
                'Info for displacements generation. The following flags are allowed:\n ' +
                '\n '.join(f'{flag_name}' for flag_name in cls._ENABLED_DISPLACEMENT_GENERATOR_FLAGS)
            ),
        )
        spec.input(
            'fc_options',
            valid_type=orm.Dict,
            required=False,
            validator=cls._validate_fc_options,
            help=(
                'Options for force constants calculation (optional). The following flags are allowed:\n ' +
                '\n '.join(f'{flag_name}' for flag_name in cls._ENABLED_FC_OPTIONS_FLAGS)
            ),
        )
        spec.input(
            'nac_parameters',
            valid_type=orm.ArrayData,
            required=False,
            validator=validate_nac,
            help='Non-analytical parameters.',
        )
        spec.input_namespace(
            'options',
            help='Options for how to run the workflow.',
        )
        spec.input(
            'options.run_parallel_nac_forces',
            valid_type=bool,
            required=False,
            non_db=True,
            default=True,
            help='Whether running nac parameters and forces calculations in parallel.',
        )
        spec.inputs.validator = validate_inputs

        spec.output_namespace(
            'supercells',
            valid_type=orm.StructureData,
            dynamic=True,
            required=False,
            help='The supercells with displacements.'
        )
        spec.output_namespace(
            'supercells_forces',
            valid_type=orm.ArrayData,
            required=True,
            help='The forces acting on the atoms of each supercell.'
        )
        spec.output_namespace(
            'supercells_energies', valid_type=orm.Float, required=False, help='The total energy of each supercell.'
        )
        spec.output(
            'output_phonopy_data',
            valid_type=PhonopyData,
            help=(
                'The phonopy data with supercells displacements, forces and (optionally)'
                'nac parameters to use in the post-processing calculation.'
            )
        )
        spec.output(
            'output_force_constants',
            valid_type=ForceConstantsData,
            required=False,
            help='The matrix of force constants computed with finite displacements.'
        )

    @classmethod
    def _validate_displacements(cls, value, _):
        """Validate the ``displacements`` input namespace."""
        if value:
            value_dict = value.get_dict()
            enabled_dict = cls._ENABLED_DISPLACEMENT_GENERATOR_FLAGS
            unknown_flags = set(value_dict.keys()) - set(enabled_dict.keys())
            if unknown_flags:
                return (
                    f"Unknown flags in 'displacements': {unknown_flags}, "
                    f'allowed flags are {cls._ENABLED_DISPLACEMENT_GENERATOR_FLAGS.keys()}.'
                )
            invalid_values = [
                value_dict[key]
                for key in value_dict.keys()
                if not (type(value_dict[key]) in enabled_dict[key] or value_dict[key] in enabled_dict[key])
            ]
            if invalid_values:
                return f'Displacement options must be of the correct type; got invalid values {invalid_values}.'

    @classmethod
    def _validate_fc_options(cls, value, _):
        """Validate the ``fc_options`` input namespace."""
        if value:
            value_dict = value.get_dict()
            enabled_dict = cls._ENABLED_FC_OPTIONS_FLAGS
            unknown_flags = set(value_dict.keys()) - set(enabled_dict.keys())
            if unknown_flags:
                return (
                    f"Unknown flags in 'fc_options': {unknown_flags}, "
                    f'allowed flags are {cls._ENABLED_FC_OPTIONS_FLAGS.keys()}.'
                )
            invalid_values = [
                value_dict[key]
                for key in value_dict.keys()
                if not (type(value_dict[key]) in enabled_dict[key] or value_dict[key] in enabled_dict[key])
            ]
            if invalid_values:
                return f'Displacement options must be of the correct type; got invalid values {invalid_values}.'

    def setup(self):
        """Setup the workflow generating the PreProcessData."""
        if 'preprocess_data' in self.inputs:
            preprocess = self.inputs.preprocess_data
            if 'displacement_generator' in self.inputs:
                preprocess = preprocess.calcfunctions.get_preprocess_with_new_displacements(
                    self.inputs.displacement_generator
                )
        else:
            preprocess_inputs = {}
            for input_ in [
                'structure', 'supercell_matrix', 'primitive_matrix', 'symprec', 'is_symmetry', 'displacement_generator'
            ]:
                if input_ in self.inputs:
                    preprocess_inputs.update({input: self.inputs[input_]})
            preprocess = PreProcessData.generate_preprocess_data(**preprocess_inputs)

        self.ctx.preprocess_data = preprocess

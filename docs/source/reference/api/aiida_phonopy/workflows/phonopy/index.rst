:py:mod:`aiida_phonopy.workflows.phonopy`
=========================================

.. py:module:: aiida_phonopy.workflows.phonopy

.. autoapi-nested-parse::

   Abstract workflow for automatic frozen phonons calculations.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.workflows.phonopy.PhonopyWorkChain



Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.workflows.phonopy.validate_matrix
   aiida_phonopy.workflows.phonopy.validate_positive_integer
   aiida_phonopy.workflows.phonopy.validate_nac
   aiida_phonopy.workflows.phonopy.validate_inputs



.. py:function:: validate_matrix(value, _)

   Validate the `supercell_matrix` and `primitive_matrix` inputs.


.. py:function:: validate_positive_integer(value, _)

   Validate that `value` is positive.


.. py:function:: validate_nac(value, _)

   Validate that `value` is positive.


.. py:function:: validate_inputs(inputs, _)

   Validate the entire inputs namespace.


.. py:class:: PhonopyWorkChain(inputs: dict | None = None, logger: logging.Logger | None = None, runner: 'Runner' | None = None, enable_persistence: bool = True)

   Bases: :py:obj:`aiida.engine.WorkChain`

   Abstract workflow for automated frozen phonons.

   Phonopy is used to produce structures with displacements,
   while the forces are calculated with a quantum engine of choice.

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

       `self.ctx.preprocess_data.calcfunctions.get_supercells_with_displacements()`

       This will return a dictionary with all the supercells as StructureData to run for the phonon calculation.
       The keys of this dictionary are of the type `supercell_{number}`, where `number` is an integer.
       These numbers are essentials since the `phonopy` force sets is generated following these numbers,
       in order to make sure to refer to the correct displacement. Thus, it is required to keep track
       of them.
       Moreover,a calculation over the pristine supercell structure should be run before hand as reference.
       This structure can instead be gotten via:

       `self.ctx.preprocess_data.calcfunctions.get_supercell()`

       This will return a StructureData without any label.

       For an example of implementation, refer to aiidateam/aiida-common-worfklows.

       * Note: some type of force calculation needs to map some variables from the unitcell to the supercell
       (and in certain case even the primitive cell), e.g. the atomic spin in VASP. Since this is code dependent,
       you will need to map these parameters before launching the force calculation of a certain supercell
       with displacement. This information can be gotten via:

       `self.ctx.preprocess_data.get_cells_mappings()`

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

           `self.ctx.preprocess_data.calcfunctions.generate_phonopy_data(**self.outputs.supercells_forces)`

           and then expose it as output in the `output_phonopy_data` namespace.

       * Alternatively: instead of exposing the supercell forces as outputs, you can directly gather all the forces
           in a dictionary and run directly to the `generate_phonopy_data` method using this dictionary (always using
           the double *).

       See the implementation in aiidateam/aiida-common-workflows for an example.

   4. (optional) Run the non-analytical constants on the primitive cell.

       Non-analytical constants should be run for polar insulators. These require usually a linear response code
       or a finite difference approach (e.g. using the electric enthalpy). Since this is usually the most expensive
       part, you should run them on the primitive cell. To get it, use:

       `self.ctx.preprocess_data.calcfunctions.get_primitive_cell()`

       If you compute also these, collect the dielectric tensor and the effectic born charges in an ArrayData,
       with the arraynames `dielectric` and `born_charges` (in Cartesian coordinates!).
       Then, gather all the information of nac and forces in a unique `PhonopyData` via:

       `self.ctx.preprocess_data.calcfunctions.generate_phonopy_data(
           nac_parameters=nac_paramters,
           **self.outputs.supercells_forces
           )`

       and expose the output.

       .. note:: we require in the input for generating the full phonopy data, to give the nac in the primitive cell.
       The primitive cell of phonopy will just rotate the lattice vectors, thus mantaining the Cartasian coordinate
       system. It can happen, though, that the unitcell is not the primitive cell of the system, meaning that the
       primitive cell will contain less atoms. We expect in input the nac computed on this number of atoms. If you
       want, for some reason, compute the nac on the unitcell, you will need to get the reduced nac.
       To do so, you can consider using a built-in function in phonopy, namely:

       :py:func:`phonopy.structure.symmetry.elaborate_borns_and_epsilon`


   .. py:attribute:: _ENABLED_DISPLACEMENT_GENERATOR_FLAGS



   .. py:attribute:: _ENABLED_FC_OPTIONS_FLAGS



   .. py:method:: define(spec)
      :classmethod:

      Define inputs, outputs, and outline.


   .. py:method:: _validate_displacements(value, _)
      :classmethod:

      Validate the ``displacements`` input namespace.


   .. py:method:: _validate_fc_options(value, _)
      :classmethod:

      Validate the ``fc_options`` input namespace.


   .. py:method:: setup()

      Set up the workflow generating the PreProcessData.

:py:mod:`aiida_phonopy.calculations.phonopy`
============================================

.. py:module:: aiida_phonopy.calculations.phonopy

.. autoapi-nested-parse::

   CalcJob for phonopy post-processing.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.calculations.phonopy.PhonopyCalculation



Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.calculations.phonopy.get_default_metadata_options



.. py:function:: get_default_metadata_options()

   Get a default metadata option for Phonopy calculation.


.. py:class:: PhonopyCalculation(*args, **kwargs)

   Bases: :py:obj:`aiida.engine.CalcJob`

   Base `CalcJob` implementation for Phonopy post-processing.

   .. py:attribute:: _OUTPUTS



   .. py:attribute:: _INPUT_FORCE_CONSTANTS
      :value: 'force_constants.hdf5'



   .. py:attribute:: _DEFAULT_INPUT_FILE
      :value: 'aiida.in'



   .. py:attribute:: _DEFAULT_OUTPUT_FILE
      :value: 'aiida.out'



   .. py:attribute:: _DEFAULT_PHONOPY_FILE
      :value: 'phonopy.yaml'



   .. py:attribute:: _INPUT_SUBFOLDER
      :value: './'



   .. py:attribute:: _OUTPUT_SUBFOLDER
      :value: './'



   .. py:attribute:: _AVAILABLE_TAGS



   .. py:attribute:: _BLOCKED_TAGS
      :value: ['DIM', 'ATOM_NAME', 'MASS', 'MAGMOM', 'CREATE_DISPLACEMENTS', 'DISPLACEMENT_DISTANCE', 'DIAG',...



   .. py:method:: define(spec)
      :classmethod:

      Define inputs, outputs, and outline.


   .. py:method:: _validate_parameters(value, _)
      :classmethod:

      Validate the ``parameters`` input namespace.


   .. py:method:: prepare_for_submission(folder)

      Prepare the calculation job for submission by transforming input nodes into input files.

      In addition to the input files being written to the sandbox folder, a `CalcInfo` instance will be returned that
      contains lists of files that need to be copied to the remote machine before job submission, as well as file
      lists that are to be retrieved after job completion.

      :param folder: a sandbox folder to temporarily write files on disk.
      :return: :py:class:`~aiida.common.datastructures.CalcInfo` instance.


   .. py:method:: _get_p2s_map()

      Get the primitive to supercell map.


   .. py:method:: write_phonopy_info(folder)

      Write in `folder` the `phonopy.yaml` file.


   .. py:method:: write_force_constants(folder)

      Write in `folder` the force constants file.


   .. py:method:: write_calculation_input(folder, parameters: dict, filename: str)

      Write in `folder` the input file containing the information regarding the calculation.

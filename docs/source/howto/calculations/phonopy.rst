
.. _howto:calculations:phonopy:

``PhonopyCalculation``
----------------------

The ``phonopy`` of Phonopy computes phonon properties within the harmonic approximation, using either finite difference of forces or already computed
second-order force constants with other methods (e.g. density functional perturbation theory, DFPT).
Examples of these properties include phonon band structure, Helmoltz free energy, thermal capacity, and so on.
For a detailed list of such properties, please consult the `official documentation <https://phonopy.github.io/phonopy/>`_.

================== ===============================================================================
Plugin class       :class:`~aiida_phonopy.calculations.phonopy.PhonopyCalculation`
Plugin entry point ``phonopy.phonopy``
================== ===============================================================================


How to launch a ``phonopy`` calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Below is a script with a basic example of how to run a ``phonopy`` calculation through the ``PhonopyCalculation``
plugin that computes the phonon band structure of an fcc silicon crystal:

.. literalinclude:: ../../tutorials/include/run_phonopy_basic.py
    :language: python

Note that you may have to change the name of the code that is loaded using ``load_code``.


How to define input file parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``phonopy`` code supports many parameters that can be defined through the input file,
as shown on the `official documentation <https://phonopy.github.io/phonopy/>`_.
The parameters are simple tags that activate a particular post-processing.

.. code-block:: python

    parameters = {'bands':'auto'}

The parameters dictionary should be wrapped in a :class:`~aiida.orm.Dict` node
and assigned to the ``parameters`` input of the process builder:

.. code-block:: python

    from aiida.orm import Dict, load_code
    builder = load_code('phonopy').get_builder()
    parameters = {
        ...
    }
    builder.parameters = Dict(parameters)

.. warning::

    There are a number of input parameters that *cannot* be set, as they will be automatically set by
    the plugin based on other inputs. Defining them anyway will result in an exception when launching the calculation.


How to add extra settings
^^^^^^^^^^^^^^^^^^^^^^^^^

You can specify extra settings to set the fine details of the calculations. These are:

    * ``symmetrize_nac``: :class:`bool`; it symmetrizes the dielectric and Born effective charges, and applies sum rules
    * ``factor_nac``: :class:`float`; it will use this conversion factor in the non-analytical correction
    * ``subtract_residual_forces``: :class:`bool`; it will subtract the forces of the pristine supercell

For example:

.. code-block:: python

    builder.settings = Dict({'symmetrize_nac': True})


How to retrieve animation files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``PhonopyCalculation`` plugin will retrieve the most important and common output files by default.
To retrieve animation output files, specify the ``keep_animation_files`` key in the ``settings`` input:

.. code-block:: python

    builder.settings = Dict({'keep_animation_files': True})

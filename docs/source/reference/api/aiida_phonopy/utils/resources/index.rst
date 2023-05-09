:py:mod:`aiida_phonopy.utils.resources`
=======================================

.. py:module:: aiida_phonopy.utils.resources

.. autoapi-nested-parse::

   Utilities for CalcJob resources.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.utils.resources.get_default_options



.. py:function:: get_default_options(max_num_machines=1, max_wallclock_seconds=300, with_mpi=False)

   Return an instance of the options dictionary with the minimally required parameters for a `CalcJob`.

   :param max_num_machines: set the number of nodes, default=1
   :param max_wallclock_seconds: set the maximum number of wallclock seconds, default=1800
   :param with_mpi: whether to run the calculation with MPI enabled

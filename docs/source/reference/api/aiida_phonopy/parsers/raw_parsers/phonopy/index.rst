:py:mod:`aiida_phonopy.parsers.raw_parsers.phonopy`
===================================================

.. py:module:: aiida_phonopy.parsers.raw_parsers.phonopy

.. autoapi-nested-parse::

   Raw parsers of the phonopy output files.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.parsers.raw_parsers.phonopy.parse_stdout



.. py:function:: parse_stdout(stdout)

   Raw parser of the phonopy std output.

   :param stdout: str of the std output of phonopy
   :returns: tuple of two dictionaries, with the parsed data and log messages, respectively

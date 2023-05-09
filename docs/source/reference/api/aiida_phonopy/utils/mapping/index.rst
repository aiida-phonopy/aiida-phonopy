:py:mod:`aiida_phonopy.utils.mapping`
=====================================

.. py:module:: aiida_phonopy.utils.mapping

.. autoapi-nested-parse::

   Utilities to deal with various mapping data structures.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.utils.mapping.get_logging_container
   aiida_phonopy.utils.mapping._case_transform_dict
   aiida_phonopy.utils.mapping._lowercase_dict
   aiida_phonopy.utils.mapping._uppercase_dict



.. py:function:: get_logging_container()

   Return an `AttributeDict` that can be used to map logging messages to certain log levels.

   This datastructure is useful to add log messages in a function that does not have access to the right logger. Once
   returned, the caller who does have access to the logger can then easily loop over the contents and pipe the messages
   through the actual logger.

   :return: :py:class:`~aiida.common.extendeddicts.AttributeDict`


.. py:function:: _case_transform_dict(dictionary: dict, dict_name: str, func_name, transform)

   Transform a dictionary to have the first keys capitalized or decapitalized.


.. py:function:: _lowercase_dict(dictionary: dict, dict_name: str)

   Transform a dictionary to have the first keys decapitalized.


.. py:function:: _uppercase_dict(dictionary: dict, dict_name: str)

   Transform a dictionary to have the first keys capitalized.

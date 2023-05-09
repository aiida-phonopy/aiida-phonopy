:py:mod:`aiida_phonopy.parsers.base`
====================================

.. py:module:: aiida_phonopy.parsers.base

.. autoapi-nested-parse::

   Defines a `Parser` base class for `aiida-phonopy`.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.parsers.base.Parser




.. py:class:: Parser(node: aiida.orm.CalcJobNode)

   Bases: :py:obj:`aiida.parsers.Parser`

   Custom `Parser` class for `aiida-phonopy` parser implementations.

   .. py:method:: emit_logs(logging_dictionaries, ignore=None)

      Emit the messages in one or multiple "log dictionaries" through the logger of the parser.

      A log dictionary is expected to have the following structure: each key must correspond to a log level of the
      python logging module, e.g. `error` or `warning` and its values must be a list of string messages. The method
      will loop over all log dictionaries and emit the messages it contains with the log level indicated by the key.

      Example log dictionary structure::

          logs = {
              'warning': ['Could not parse the `etot_threshold` variable from the stdout.'],
              'error': ['Self-consistency was not achieved']
          }

      :param logging_dictionaries: log dictionaries
      :param ignore: list of log messages to ignore


   .. py:method:: exit(exit_code)

      Log the exit message of the give exit code with level `ERROR` and return the exit code.

      This is a utility function if one wants to return from the parse method and automically add the exit message
      associated to the exit code as a log message to the node: e.g. `return self.exit(self.exit_codes.LABEL))`

      :param exit_code: an `ExitCode`
      :return: the exit code

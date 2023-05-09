:py:mod:`aiida_phonopy.parsers.phonopy`
=======================================

.. py:module:: aiida_phonopy.parsers.phonopy

.. autoapi-nested-parse::

   Parsers of `PhonopyCalculation` output files.



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   aiida_phonopy.parsers.phonopy.PhonopyParser



Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.parsers.phonopy.file_opener



.. py:function:: file_opener(folder_path: str, filename: str)

   Return a function to open different expected output format.


.. py:class:: PhonopyParser(node: aiida.orm.CalcJobNode)

   Bases: :py:obj:`aiida_phonopy.parsers.base.Parser`

   Parser the files produced by a phonopy post processing calculation.

   .. py:method:: parse(**kwargs)

      Parse retrieved files from remote folder.


   .. py:method:: parse_stdout() -> tuple

      Parse the stdout output file.

      :param parameters: the input parameters dictionary
      :param parsed_phonopy: the collected parsed data from the yaml phonopy output
      :return: raw parsed data


   .. py:method:: get_expected_filenames_keys() -> set

      Return the retrieve file keys (that map to the filenames outputs) depending on the tags in `parameters`.


   .. py:method:: parse_force_constants(filepath: str) -> aiida.orm.ArrayData

      Parse the `force_constants.hdf5` output file.


   .. py:method:: load_with_numpy(file: str) -> numpy.ndarray

      Load a txt file using numpy.


   .. py:method:: load_with_yaml(file: str) -> dict

      Load a yaml file using.


   .. py:method:: parse_yaml(file: str) -> aiida.orm.Dict

      Parse a `.yaml` file and return it as a Dict.


   .. py:method:: parse_total_dos(file: str) -> aiida.orm.XyData

      Parse `total_dos.dat` output file.


   .. py:method:: parse_projected_dos(file: str) -> aiida.orm.XyData

      Parse `projected_dos.dat` output file.


   .. py:method:: parse_thermal_properties(file: str) -> aiida.orm.XyData

      Parse the `thermal_properties.yaml`` output file.


   .. py:method:: parse_band_structure(file: str, freqs_units: str = 'THz') -> aiida.orm.BandsData

      Parse the `band.hdf5`` output file.

      Expected keys are:
          * **nqpoint**: array with total number of q-points in the band structure, array(1,)
          * **frequency**: array of frequencies at each q-point; array(npath, segment_nqpoint, nband)
          * **label**: array of labels, two per each path; array(npath, 2)
          * **path**: number of q-points per path; array(npath, segment_nqpoint, qpoint/qposition(3,))
          * **distance**: distance between consecutive q-points in a segment of path; array(npath, segment_nqpoint)
          * **segment_nqpoint**: number of q-points per segment of the entire path; array(npath,)
          * **eigenvector**: (optional) eigenvector of each phonon mode, i.e. each frequency of each q-point;
              array(npath, segment_nqpoint, nband, eigenvector(6,))
          * **group_velocity**: (optional) group velocity of each phonon mode (as for eigenvector);
              array(npath, segment_nqpoint, nband, group_velocity(3,))


   .. py:method:: parse_qpoints(file: str, freqs_units: str = 'THz') -> aiida.orm.BandsData

      Parse the `mesh.hdf5`` and `qpoints.hdf5`` output files.

      Expected keys are:
          * **frequency**: array of frequencies at each q-point; array(nqpoint, nband)
          * **mesh**: qpoint mesh; array(3,)
          * **qpoint**: qpoints; array(nqpoint, 3)
          * **weight**: weight of qpoints; array(nqpoint,)
          * **eigenvector**: (optional) eigenvector of each phonon mode, i.e. each frequency of each q-point;
              array(npath, segment_nqpoint, nband, eigenvector(6,))
          * **group_velocity**: (optional) group velocity of each phonon mode (as for eigenvector);
              array(npath, segment_nqpoint, nband, group_velocity(3,))


   .. py:method:: _get_p2s_map()

      Get the primitive to supercell map.

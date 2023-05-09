:py:mod:`aiida_phonopy.calculations.functions.link_structures`
==============================================================

.. py:module:: aiida_phonopy.calculations.functions.link_structures

.. autoapi-nested-parse::

   Functions for linking PhonopyAtoms and StructureData.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   aiida_phonopy.calculations.functions.link_structures.phonopy_atoms_to_structure
   aiida_phonopy.calculations.functions.link_structures.phonopy_atoms_from_structure
   aiida_phonopy.calculations.functions.link_structures.if_to_map



.. py:function:: phonopy_atoms_to_structure(cell: phonopy.structure.cells.PhonopyAtoms, mapping: dict | None = None, pbc: tuple[bool, bool, bool] = (True, True, True)) -> aiida.orm.StructureData

   Return a StructureData from a PhonopyAtoms instance.

   :param cell: a PhonopyAtoms instance
   :param mapping: a number to kinds and symbols map, defaults to None
   :param pbc: periodic boundary conditions in the three lattice directions


.. py:function:: phonopy_atoms_from_structure(structure: aiida.orm.StructureData) -> phonopy.structure.cells.PhonopyAtoms

   Return a tuple containg the PhonopyAtoms and the mapping from a StructureData.

   :param structure: StructureData instance
   :return: tuple with element:
       * a :class:`~phonopy.structure.cells.PhonopyAtoms` istance
       * mapping dictionary, with key:pair of the type int:str, string
           referring to the custom atomic name


.. py:function:: if_to_map(structure: aiida.orm.StructureData) -> bool

   Return a bool according whether kind names and symbols are the same.

   :param structure: StructureData instance

   :return: True if kind names different from symbols are used, False otherwise

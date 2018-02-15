from aiida.orm.data.array import ArrayData
from aiida.orm.data.structure import StructureData

import numpy


class NacData(ArrayData):
    """
    Store the force constants on disk as a numpy array. It requires numpy to be installed.
    """

    def __init__(self, *args, **kwargs):
        super(NacData, self).__init__(*args, **kwargs)
        self._cached_arrays = {}

    def set_structure(self, structure):
        """
        Sets the structure used to calculate the Born effective charges and the dielectric tensor (NAC parameters).
        Setting this structure is necessary for get_born_parameters_phonopy function to work.

        :param structure: StructureData object that constains the unitcell used to calculate the NAC parameters
        :return:
        """

        self._set_attr('cell', structure.cell)
        self._set_attr('positions', [site.position for site in structure.sites])
        self._set_attr('symbols', [site.kind_name for site in structure.sites])

    def get_epsilon(self):
        """
        Return dielectric tensor as a numpy array
        """

        return self.get_array('epsilon')

    def get_born_charges(self):
        """
        Return born charges  as a numpy array
        """

        return self.get_array('born_charges')

    def set_born_charges(self, born_charges):
        """
        Store Born charges as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.

        :param born_charges: The numpy array to store.
        """

        self.set_array('born_charges', numpy.array(born_charges))

    def set_epsilon(self, epsilon):
        """
        Store the dielectric tensor as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.

        :param epsilon: The numpy array to store.
        """

        self.set_array('epsilon', numpy.array(epsilon))

    def get_structure(self):

        """
        Returns the structure used to calculate the NAC parameters

        :return: StructureData
        """

        structure = StructureData(cell=self.get_attr('cell'))

        symbols = self.get_attr('symbols')
        positions = self.get_attr('positions')
        for symbol, position in zip(symbols, positions):
            structure.append_atom(position=position,
                                  symbols=symbol)

        return structure

    def get_born_parameters_phonopy(self, primitive_cell=None, symprec=1e-5):
        """
        Returns a dictionary in phonopy format that contains the non-analytical corrections.
        By default the born charges for all the atoms in the unit cell used to calculate them.
        Using primitive_cell allow to choose to use a custom unit cell.

        :param primitive_cell (optional): (NumpyArray) lattice vectors matrix that define the unit cell in which the born_charges are returned. (this should be the primitive cell used in phonopy)
        :return:
        """

        import numpy as np
        from phonopy.structure.cells import get_primitive, get_supercell
        from phonopy.structure.atoms import Atoms as PhonopyAtoms
        from phonopy.units import Hartree, Bohr

        born_charges = self.get_array('born_charges')
        epsilon = self.get_array('epsilon')
        structure_born = self.get_structure()

        ucell = PhonopyAtoms(symbols=[site.kind_name for site in structure_born.sites],
                             positions=[site.position for site in structure_born.sites],
                             cell=structure_born.cell)

        if primitive_cell is None:
            target_mat = np.identity(3)
        else:
            inv_target = np.linalg.inv(primitive_cell)
            target_mat = np.dot(inv_target, structure_born.cell)

        if np.linalg.det(structure_born.cell) < np.linalg.det(primitive_cell):

            # inv_pmat = np.identity(3)
            # pcell = get_primitive(ucell, inv_pmat, symprec=symprec)

            inv_smat = np.linalg.inv(target_mat)

            scell = get_supercell(ucell, inv_smat, symprec=symprec)

            s2p = scell.get_supercell_to_unitcell_map()
            map_primitive = scell.get_unitcell_to_unitcell_map()

            reduced_born = [born_charges[map_primitive[i]] for i in s2p]

        else:

            # inv_smat = np.identity(3)
            # scell = get_supercell(ucell, inv_smat, symprec=symprec)

            inv_pmat = np.linalg.inv(target_mat)

            pcell = get_primitive(ucell, inv_pmat, symprec=symprec)

            s2p = pcell.get_supercell_to_primitive_map()
            map_primitive = pcell.get_primitive_to_primitive_map()

            reduced_born = [born_charges[map_primitive[i]] for i in s2p]

        factor = Hartree * Bohr
        non_anal = {'born': reduced_born,
                    'factor': factor,
                    'dielectric': epsilon}

        return non_anal

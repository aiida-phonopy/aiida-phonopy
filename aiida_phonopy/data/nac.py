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

        self._set_attr('cell', structure.cell)
        self._set_attr('positions', [site.position for site in structure.sites])
        self._set_attr('symbols', [site.kind_name for site in structure.sites])

    def get_epsilon(self):
        """
        Return dielectric tensor stored in the node as a numpy array
        """

        return self.get_array('epsilon')

    def get_born_charges(self):
        """
        Return born charges stored in the node as a numpy array
        """

        return self.get_array('born_charges')

    def set_born_charges(self, born_charges):
        """
        Store Born charges as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.

        :param array: The numpy array to store.
        """

        self.set_array('born_charges', numpy.array(born_charges))

    def set_epsilon(self, epsilon):
        """
        Store the dielectric tensor as a numpy array. Possibly overwrite the array
        if it already existed.
        Internally, it is stored as a force_constants.npy file in numpy format.

        :param array: The numpy array to store.
        """

        self.set_array('epsilon', numpy.array(epsilon))

    def get_structure(self):

        structure = StructureData(cell=self.get_attr('cell'))

        symbols = self.get_attr('symbols')
        positions = self.get_attr('positions')
        for symbol, position in zip(symbols, positions):
            structure.append_atom(position=position,
                                  symbols=symbol)

        return structure

    def get_born_parameters_phonopy(self, symprec=1e-5):

        from phonopy.units import Hartree, Bohr
        from phonopy.structure.atoms import Atoms as PhonopyAtoms
        from phonopy.interface.vasp import symmetrize_borns_and_epsilon
        from phonopy.structure.cells import get_supercell, get_primitive

        import numpy as np

        pmat = self.get_array('primtive_matrix')
        smat = self.get_array('supercell_matrix')

        cell = PhonopyAtoms(symbols=self.get_attr('symbols'),
                            positions=self.get_attr('positions'),
                            cell=self.get_attr('cell'))

        borns_, epsilon_ = symmetrize_borns_and_epsilon(self.get_born_charges(),
                                                        self.get_epsilon(),
                                                        cell,
                                                        symprec=symprec)

        scell = get_supercell(cell, smat, symprec=symprec)
        u2u_map = scell.get_unitcell_to_unitcell_map()

        pcell = get_primitive(scell, np.dot(np.linalg.inv(smat), pmat), symprec=symprec)

        p2s_map = pcell.get_primitive_to_supercell_map()
        born_primitive = borns_[[u2u_map[i] for i in p2s_map]]

        factor = Hartree * Bohr

        non_anal = {'born': born_primitive,
                    'factor': factor,
                    'dielectric': epsilon_}
        return non_anal
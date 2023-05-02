# -*- coding: utf-8 -*-
"""Calcfunctions Utils for aiida-phonopy DataTypes."""
from aiida import orm
from aiida.engine import calcfunction
from aiida.plugins import DataFactory

PreProcessData = DataFactory('phonopy.preprocess')
PhonopyData = DataFactory('phonopy.phonopy')

__all__ = (
    'get_unitcell', 'get_primitive', 'get_supercell', 'get_supercells_with_displacements', 'get_displacements',
    'get_preprocess_with_new_displacements', 'generate_preprocess_data', 'generate_phonopy_data', 'CalcfunctionMixin'
)


@calcfunction
def get_unitcell(preprocess_data: PreProcessData):
    """Get the unitcell of a PreProcessData as a StructureData."""
    structure_data = preprocess_data.get_unitcell()
    return structure_data


@calcfunction
def get_primitive(preprocess_data: PreProcessData):
    """Get the primitive cell of a PreProcessData as a StructureData."""
    structure_data = preprocess_data.get_primitive_cell()
    return structure_data


@calcfunction
def get_supercell(preprocess_data: PreProcessData):
    """Get the supercell (pristine) of a PreProcessData as a StructureData."""
    structure_data = preprocess_data.get_supercell()
    return structure_data


@calcfunction
def get_supercells_with_displacements(preprocess_data: PreProcessData):
    """Get the supercells with displacements of a PreProcessData as a StructureData."""
    structures_data = preprocess_data.get_supercells_with_displacements()
    return structures_data


@calcfunction
def get_displacements(preprocess_data: PreProcessData):
    """Get the displacements of a PreProcessData as an ArrayData with array name `displacements`."""
    displacements = preprocess_data.get_displacements()
    the_displacements = orm.ArrayData()
    the_displacements.set_array('displacements', displacements)

    return the_displacements


@calcfunction
def generate_preprocess_data(
    structure: orm.StructureData,
    displacement_generator=None,
    supercell_matrix=None,
    primitive_matrix=None,
    symprec=None,
    is_symmetry=None,
    distinguish_kinds=None,
):
    """Returns a complete stored PreProcessData node.

    :param structure: structure data node representing the unitcell
    :type structure: orm.StructureData
    :param displacement_generator: dictionary containing the info for generating the displacements
    :type displacement_generator: orm.Dict
    :param supercell_matrix: supercell matrix, defaults to diag(1,1,1)
    :type supercell_matrix: orm.List, optional
    :param primitive_matrix: primitive matrix, defaults to "auto"
    :type primitive_matrix: orm.List or orm.List, optional
    :param symprec: symmetry precision on atoms, defaults to 1e-5
    :type symprec: orm.Float, optional
    :param is_symmetry: if using space group symmetry, defaults to True
    :type is_symmetry: orm.Bool, optional
    :param distinguish_kinds: if distinguish names of same specie by symmetry, defaults to True
    :type distinguish_kinds: orm.Bool, optional

    :return: PreProcessData node
    """
    kwargs = {}

    kwargs['structure'] = structure

    if displacement_generator is not None:
        displ_generator = displacement_generator.get_dict()
    else:
        displ_generator = {
            'distance': 0.01,
            'is_plusminus': 'auto',
            'is_diagonal': True,
            'is_trigonal': False,
            'number_of_snapshots': None,
            'random_seed': None,
            'temperature': None,
            'cutoff_frequency': None,
        }

    if primitive_matrix is not None:
        if isinstance(primitive_matrix, orm.List):
            kwargs['primitive_matrix'] = primitive_matrix.get_list()
        if isinstance(primitive_matrix, orm.Str):
            kwargs['primitive_matrix'] = primitive_matrix.value
    else:
        kwargs['primitive_matrix'] = 'auto'

    if supercell_matrix is not None:
        kwargs['supercell_matrix'] = supercell_matrix.get_list()

    if symprec is not None:
        kwargs['symprec'] = symprec.value

    if is_symmetry is not None:
        kwargs['is_symmetry'] = is_symmetry.value

    if distinguish_kinds is not None:
        kwargs['distinguish_kinds'] = distinguish_kinds.value

    preprocess = PreProcessData(**kwargs)

    preprocess.set_displacements(**displ_generator)

    return preprocess


@calcfunction
def get_preprocess_with_new_displacements(preprocess_data: PreProcessData, displacement_generator: orm.Dict):
    """Get a new PreProcessData from an old one from new displacement generator settings."""
    displacement_dataset = preprocess_data.generate_displacement_dataset(**displacement_generator.get_dict())

    preprocess = PreProcessData(
        structure=preprocess_data.get_unitcell(),
        supercell_matrix=preprocess_data.supercell_matrix,
        primitive_matrix=preprocess_data.primitive_matrix,
        symprec=preprocess_data.symprec,
        is_symmetry=preprocess_data.is_symmetry,
        distinguish_kinds=preprocess_data.distinguish_kinds,
    )

    preprocess.set_displacements_from_dataset(displacement_dataset)

    return preprocess


@calcfunction
def generate_phonopy_data(preprocess_data: PreProcessData, nac_parameters=None, forces_index=None, **forces_dict):
    """Create a PhonopyData node from a PreProcess(Phonopy)Data node, storing forces and (optionally)
        non-analytical constants.

        `Forces` must be passed as **kwargs**, since we are calling a calcfunction with a variable
        number of supercells forces.

        :param nac_parameters: ArrayData containing 'dielectric' and 'born_charges' as arrays
            with their correct shape
        :param forces_index: Int if a TrajectoryData is given, in order to get the correct slice of the array.
        :param forces_dict: dictionary of supercells forces as ArrayData stored as `forces`, each Data
            labelled in the dictionary in the format `forces_{suffix}`.
            The prefix is common and the suffix corresponds to the suffix number of the supercell with
            displacement label given from the `get_supercells_with_displacements` method.

            For example:
                {'forces_1':ArrayData, 'forces_2':ArrayData}
                <==>
                {'supercell_1':StructureData, 'supercell_2':StructureData}
                and forces in each ArrayData stored as 'forces',
                i.e. ArrayData.get_array('forces') must not raise error

            .. note: if residual forces would be stored, label it with 0 as suffix.
    """
    prefix = 'forces'

    forces_0 = forces_dict.pop(f'{prefix}_0', None)
    # Setting the dictionary of forces
    dict_of_forces = {}

    for key, value in forces_dict.items():
        if key.startswith(prefix):
            dict_of_forces[key] = value.get_array('forces')

    if forces_index is not None:
        forces_index = forces_index.value

    # Setting data on a new PhonopyData
    new_phonopy_data = PhonopyData(preprocess_data=preprocess_data)
    new_phonopy_data.set_forces(dict_of_forces=dict_of_forces, forces_index=forces_index)

    if forces_0 is not None:
        new_phonopy_data.set_residual_forces(forces=forces_0.get_array('forces'))

    if nac_parameters is not None:
        new_phonopy_data.set_dielectric(nac_parameters.get_array('dielectric'))
        new_phonopy_data.set_born_charges(nac_parameters.get_array('born_charges'))

    return new_phonopy_data


class CalcfunctionMixin:
    """Set of calcfunctions to be called from the aiida-phonopy DataTypes."""

    def __init__(self, data_node):
        self._data_node = data_node

    def get_unitcell(self):
        """Get the unitcell as a StructureData through a calfunction."""
        return get_unitcell(preprocess_data=self._data_node)

    def get_primitive_cell(self):
        """Get the primitive cell as a StructureData through a calfunction."""
        return get_primitive(preprocess_data=self._data_node)

    def get_supercell(self):
        """Get the supercell (pristine) as a StructureData through a calfunction."""
        return get_supercell(preprocess_data=self._data_node)

    def get_supercells_with_displacements(self):
        """Get the supercells with displacements as a StructureData through a calfunction."""
        return get_supercells_with_displacements(preprocess_data=self._data_node)

    def get_displacements(self):
        """Get the displacements as an ArrayData through a calfunction."""
        return get_displacements(preprocess_data=self._data_node)

    def get_preprocess_with_new_displacements(self, displacement_generator: orm.Dict):
        """Create a PreProcessData node from a PreProcess/PhonopyData with a new set of displacements.

        :param displacement_generator: a `storable` dictionary  """
        return get_preprocess_with_new_displacements(
            preprocess_data=self._data_node, displacement_generator=displacement_generator
        )

    def generate_phonopy_data(self, nac_parameters=None, forces_index=None, **forces_dict):
        """Create a PhonopyData node from a PreProcess(Phonopy)Data node, storing forces and (optionally)
        non-analytical constants.

        `Forces` must be passed as **kwargs**, since we are calling a calcfunction with a variable
        number of supercells forces.

        :param nac_parameters: ArrayData containing 'dielectric' and 'born_charges' as arrays
            with their correct shape
        :param forces_index: Int if a TrajectoryData is given, in order to get the correct slice of the array.
        :param forces_dict: dictionary of supercells forces as ArrayData stored as `forces`, each Data
            labelled in the dictionary in the format `{prefix}_{suffix}`.
            The prefix is common and the suffix corresponds to the suffix number of the supercell with
            displacement label given from the `get_supercells_with_displacements` method

            For example:
                {'forces_1':ArrayData, 'forces_2':ArrayData} goes along with
                {'supercell_1':StructureData, 'supercell_2':StructureData}
                and forces in each ArrayData stored as 'forces',
                i.e. ArrayData.get_array('forces') must not raise error
        """
        return generate_phonopy_data(
            preprocess_data=self._data_node, nac_parameters=nac_parameters, forces_index=forces_index, **forces_dict
        )

# -*- coding: utf-8 -*-
"""Calcfunctions to gather computed forces into force sets."""

# from aiida import orm
# from aiida.engine import calcfunction
# import numpy as np

# @calcfunction
# def get_dataset(displacement_dataset, **forces_dict):
#     """Return force sets from supercell forces.

#     :param forces_dict: must contain keys ending in `_{int}`;
#     each key must refer to an `ArrayData` with the forces array stored with the label `forces`.
#     Forces key related to the pristine supercell must finish in `_0`, to subtract the
#     residual forces.
#     :return: phonopy dataset dictionary, ready postprocess use
#     :raises: `ValueError` if the keys are not in the correct format, `RuntimeError` if
#     the number of displacements (i.e. supercells with displacements) in the `displacement_dataset` is different from
#     the number of forces in `forces_dict` ()
#     """
#     for key in forces_dict.keys():
#         try:
#             int(key.split('_')[-1])
#             # dumb way of getting the prefix including cases of multiple `_`, e.g. `force_calculation_001`
#             prefix = key[:-(len(key.split('_')[-1]) + 1)]
#         except ValueError:
#             raise ValueError(f'{key} is not an acceptable key, must finish with an int number.')

#     forces_0 = forces_dict.pop(f'{prefix}_0', None)

#     # Setting force sets array
#     force_sets = [0 for _ in forces_dict.keys()]

#     # Filling arrays in numeric order determined by the key label
#     for key, value in forces_dict.items():
#         index = int(key.split('_')[-1])
#         force_sets[index - 1] = value.get_array('forces')

#     # Finilizing force sets array
#     force_sets = np.array(force_sets)
#     if forces_0 is not None:
#         force_sets = force_sets - forces_0.get_array('forces')

#     displacement_data = displacement_dataset.get_dict()
#     # Check if the given `displacement_dataset` and `supercells_forces` match.
#     # If not, we raise an error.
#     if not len(displacement_data['first_atoms']) == len(force_sets):
#         raise RuntimeError('`displacement_dataset` and `supercells_forces` do not match')

#     for displacement, force in zip(displacement_data['first_atoms'], force_sets):
#         displacement.update({'forces': np.array(force)})

#     return orm.Dict(dict=displacement_data)

# @calcfunction
# def get_force_constants(preprocess_info, dataset: orm.Dict):
#     """Compute and return the force constants matrix.

#     :param structure: StructureData of the unitcell
#     :param preprocess_info: the preprocess info containing supercell matrix and symmetry tolerance
#     :param dataset: phonopy dataset with displacements and forces
#     :return: ArrayData with force constants
#     """
#     from aiida_phonopy.utils.mapping import get_phonopy_instance

#     ph = get_phonopy_instance(preprocess_info)

#     ph.dataset = dataset.get_dict()
#     ph.produce_force_constants(**preprocess_info.get_attribute('fc_options'))

#     fc = orm.ArrayData()
#     fc.set_array('force_constants', ph.force_constants)
#     fc.set_array('p2s_map', preprocess_info.get_attribute('symmetry')['cells_maps']['p2s_map'])

#     return fc

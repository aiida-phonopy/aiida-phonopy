# -*- coding: utf-8 -*-
"""
This module defines a subclass of `BandsData` to include other quantities.
"""

# from aiida.orm.nodes.data.array import BandsData
# import numpy as np

# __all__ = ('ResolvedBandsData')

# class ResolvedBandsData(BandsData):
#     """
#     This class extends the base `BandsData` to include also the eigenvectors and
#     group velocities (e.g. for phonons). The former is used to estimate the band
#     connections to resolve the phonon band structure.
#     """

#     @property
#     def eigenvectors(self):
#         """Get the the phonon eigenvectors.

#         :return: (nkpoints, nbands, 6) shape numpy.ndarray"""
#         try:
#             the_eigenvectors = np.array(self.get_array('eigenvectors'))
#         except AttributeError:
#             the_eigenvectors = None
#         return the_eigenvectors

#     @eigenvectors.setter
#     def eigenvectors(self, value):
#         """Set the phonon eigenvectors.

#         :param value: (nkpoints, nbands, 6) shape numpy.ndarray
#         """
#         self.set_eigenvectors(value)

#     def set_eigenvectors(self, value):
#         """Set the phonon eigenvectors.

#         :param value: (nkpoints, nbands, 6) shape numpy.ndarray
#         """
#         self._if_can_modify()

#         try:
#             kpoints = self.get_kpoints()
#         except AttributeError:
#             raise AttributeError('Must first set the kpoints, then the bands, and then the eigentvectors')

#         try:
#             bands = self.get_bands()
#             if len(bands) != 1:
#                 bands = bands[0]
#         except AttributeError:
#             raise AttributeError('Must first set the bands, then the eigentvectors')

#         the_eigenvectors = np.array(value)

#         if len(the_eigenvectors.shape) != 3:
#             raise ValueError(f'Eigenvectors lenght must be 3, not {len(the_eigenvectors.shape)}')

#         if the_eigenvectors.shape[0] != len(kpoints):
#             raise ValueError('There must be eigenvectors for every kpoint')

#         if the_eigenvectors.shape[1] != len(bands):
#             raise ValueError('There must be eigenvectors for every energy in the bands')

#         self.set_array('eigenvectors', the_eigenvectors)

#     def _if_can_modify(self):
#         """Check if the object is stored and raise an error if so. To use in every setter."""
#         from aiida.common.exceptions import ModificationNotAllowed

#         if self.is_stored:
#             raise ModificationNotAllowed('The object cannot be modified, it has already been stored')

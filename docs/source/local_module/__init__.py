# -*- coding: utf-8 -*-
"""A module of code related to the tutorial."""
import os
import pathlib
import warnings

os.environ['AIIDA_PATH'] = str(pathlib.Path(__file__).parent / '_aiida_path')
# load the configuration without emitting a warning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', UserWarning)
    from aiida.manage.configuration import settings

    settings

from .temp_profile import load_temp_profile  # noqa: F401

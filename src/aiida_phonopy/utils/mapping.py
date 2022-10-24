# -*- coding: utf-8 -*-
"""Utilities to deal with various mapping data structures."""


def get_logging_container():
    """Return an `AttributeDict` that can be used to map logging messages to certain log levels.

    This datastructure is useful to add log messages in a function that does not have access to the right logger. Once
    returned, the caller who does have access to the logger can then easily loop over the contents and pipe the messages
    through the actual logger.

    :return: :py:class:`~aiida.common.extendeddicts.AttributeDict`
    """
    from aiida.common import AttributeDict

    return AttributeDict({
        'debug': [],
        'info': [],
        'warning': [],
        'error': [],
        'critical': [],
    })


def _case_transform_dict(dictionary, dict_name, func_name, transform):
    from collections import Counter

    from aiida.common import InputValidationError

    if not isinstance(dictionary, dict):
        raise TypeError(f'{func_name} accepts only dictionaries as argument, got {type(dictionary)}')
    new_dict = dict((transform(str(k)), v) for k, v in dictionary.items())
    if len(new_dict) != len(dictionary):
        num_items = Counter(transform(str(k)) for k in dictionary.keys())
        double_keys = ','.join([k for k, v in num_items if v > 1])
        message = (
            f'Inside the dictionary `{dict_name}` there are the following keys that '
            f'are repeated more than once when compared case-insensitively: {double_keys}.'
            'This is not allowed.'
        )
        raise InputValidationError(message)
    return new_dict


def _lowercase_dict(dictionary, dict_name):
    return _case_transform_dict(dictionary, dict_name, '_lowercase_dict', str.lower)


def _uppercase_dict(dictionary, dict_name):
    return _case_transform_dict(dictionary, dict_name, '_uppercase_dict', str.upper)

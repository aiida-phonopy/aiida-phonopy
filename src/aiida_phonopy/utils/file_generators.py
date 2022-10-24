# -*- coding: utf-8 -*-
"""Methods for writing needed `CalcJob` input files."""
from phonopy.file_IO import get_FORCE_CONSTANTS_lines


def get_FORCE_CONSTANTS_txt(force_constants):
    """Return the FORCE_CONSTANTS file in str."""
    fc = force_constants.get_array('force_constants')
    p2s_map = force_constants.get_array('p2s_map')
    lines = get_FORCE_CONSTANTS_lines(fc, p2s_map=p2s_map)

    return '\n'.join(lines)


# def get_phonopy_yaml_txt(preprocess_info, nac_parameters=None, nac_settings=None, dataset=None):
#     """Generate the `phonopy.yaml` input file."""
#     from phonopy.interface.phonopy_yaml import PhonopyYaml

#     ph = get_phonopy_instance(preprocess_info)

#     # Setting the phonopy yaml obtject to produce yaml lines
#     phpy_yaml = PhonopyYaml()
#     phpy_yaml.set_phonon_info(ph)

#     return str(phpy_yaml)

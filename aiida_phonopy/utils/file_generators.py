# -*- coding: utf-8 -*-
"""Methods for writing needed `CalcJob` input files."""

from io import StringIO
from phonopy import Phonopy
from phonopy.structure.cells import get_supercell
from phonopy.structure.dataset import forces_in_dataset
from phonopy.interface.phonopy_yaml import PhonopyYaml
from phonopy.file_IO import get_FORCE_SETS_lines, get_BORN_lines, get_FORCE_CONSTANTS_lines

from aiida_phonopy.calculations.functions.preprocess import phonopy_atoms_from_structure
from aiida import orm


def get_BORN_txt(nac_data: orm.ArrayData, structure: orm.StructureData, symmetry_tolerance: orm.Float):
    """Returns a string of BORN file. Structure is assumed as the primitive one."""

    born_charges = nac_data.get_array("born_charges")
    epsilon = nac_data.get_array("epsilon")
    # in the new implementation we expect to pass the already mapped
    # StructureData, do we don't have to have the mapping, so we
    # leave it to the default value, i.e. `False`.
    pcell = phonopy_atoms_from_structure(structure)
    lines = get_BORN_lines(pcell, born_charges, epsilon, symprec=symmetry_tolerance.value)

    return "\n".join(lines)


def get_FORCE_SETS_txt(dataset, force_sets=None):
    if dataset is None:
        return None
    if force_sets is None and not forces_in_dataset(dataset):
        return None
    if dataset is not None and force_sets is not None:
        forces = force_sets.get_array("force_sets")
        lines = get_FORCE_SETS_lines(dataset, forces=forces)
    elif dataset is not None and force_sets is None:
        lines = get_FORCE_SETS_lines(dataset)

    return "\n".join(lines)


def get_FORCE_CONSTANTS_txt(force_constants_object):
    force_constants = force_constants_object.get_array("force_constants")
    p2s_map = force_constants_object.get_array("p2s_map")
    lines = get_FORCE_CONSTANTS_lines(force_constants, p2s_map=p2s_map)

    return "\n".join(lines)


def get_phonopy_yaml_txt(structure, supercell_matrix=None, primitive_matrix="auto"):
    from aiida_phonopy.calculations.functions.preprocess import phonopy_atoms_from_structure

    unitcell, _ = phonopy_atoms_from_structure(structure)
    ph = Phonopy(
        unitcell,
        supercell_matrix=supercell_matrix,
        primitive_matrix=primitive_matrix,
    )
    phpy_yaml = PhonopyYaml()
    phpy_yaml.set_phonon_info(ph)

    return str(phpy_yaml)


def get_disp_fc3_txt(structure, parameters_data, force_sets):
    from phono3py.file_IO import write_cell_yaml

    dataset = force_sets.get_datasets3()
    supercell = get_supercell(
        phonopy_atoms_from_structure(structure),
        parameters_data.dict.supercell,
        symprec=parameters_data.dict.symmetry_tolerance,
    )

    w = StringIO.StringIO()
    w.write("natom: %d\n" % dataset["natom"])

    num_first = len(dataset["first_atoms"])
    w.write("num_first_displacements: %d\n" % num_first)
    if "cutoff_distance" in dataset:
        w.write("cutoff_distance: %f\n" % dataset["cutoff_distance"])

    num_second = 0
    num_disp_files = 0
    for d1 in dataset["first_atoms"]:
        num_disp_files += 1
        num_second += len(d1["second_atoms"])
        for d2 in d1["second_atoms"]:
            if "included" in d2:
                if d2["included"]:
                    num_disp_files += 1
            else:
                num_disp_files += 1

    w.write("num_second_displacements: %d\n" % num_second)
    w.write("num_displacements_created: %d\n" % num_disp_files)

    w.write("first_atoms:\n")
    count1 = 1
    count2 = num_first + 1
    for disp1 in dataset["first_atoms"]:
        disp_cart1 = disp1["displacement"]
        w.write("- number: %5d\n" % (disp1["number"] + 1))
        w.write("  displacement:\n")
        w.write("    [%20.16f,%20.16f,%20.16f ] # %05d\n" % (disp_cart1[0], disp_cart1[1], disp_cart1[2], count1))
        w.write("  second_atoms:\n")
        count1 += 1

        included = None
        atom2 = -1
        for disp2 in disp1["second_atoms"]:
            if atom2 != disp2["number"]:
                atom2 = disp2["number"]
                if "included" in disp2:
                    included = disp2["included"]
                pair_distance = disp2["pair_distance"]
                w.write("  - number: %5d\n" % (atom2 + 1))
                w.write("    distance: %f\n" % pair_distance)
                if included is not None:
                    if included:
                        w.write("    included: %s\n" % "true")
                    else:
                        w.write("    included: %s\n" % "false")
                w.write("    displacements:\n")

            disp_cart2 = disp2["displacement"]
            w.write("    - [%20.16f,%20.16f,%20.16f ] # %05d\n" % (disp_cart2[0], disp_cart2[1], disp_cart2[2], count2))
            count2 += 1

    write_cell_yaml(w, supercell)
    w.seek(0)
    lines = w.read()
    w.close()
    return lines


def get_forces_txt(force_sets):
    w = StringIO.StringIO()
    from phono3py.file_IO import write_FORCES_FC3

    write_FORCES_FC3(force_sets.get_datasets3(), force_sets.get_forces3(), fp=w)
    w.seek(0)
    lines = w.read()
    w.close()
    return lines


def write_fc2_to_hdf5_file(force_constants, filename):
    from phono3py.file_IO import write_fc2_to_hdf5

    write_fc2_to_hdf5(force_constants.get_data(), filename)


def write_fc3_to_hdf5_file(force_constants, filename):
    from phono3py.file_IO import write_fc3_to_hdf5

    write_fc3_to_hdf5(force_constants.get_data(), filename)


def write_kappa_to_hdf5_file(gp, filename="kappa"):
    import h5py

    with h5py.File(filename, "w") as w:
        for key in gp:
            w.create_dataset(key, data=gp[key])

import sys
from phonopy import Phonopy
from aiida_phonopy.common.utils import phonopy_atoms_from_structure
from aiida.orm import load_node
from aiida import load_profile

load_profile()


def dump_phonopy(pk):
    n = load_node(pk)
    unitcell = phonopy_atoms_from_structure(n.inputs.structure)
    smat = n.outputs.phonon_setting_info["supercell_matrix"]
    ph = Phonopy(unitcell, smat, primitive_matrix="auto")
    force_sets = n.outputs.force_sets.get_array("force_sets")
    dataset = n.outputs.phonon_setting_info["displacement_dataset"]
    ph.dataset = dataset
    ph.forces = force_sets
    if "nac_params" in n.outputs:
        borns = n.outputs.nac_params.get_array("born_charges")
        epsilon = n.outputs.nac_params.get_array("epsilon")
        nac_params = {"born": borns, "factor": 14.399652, "dielectric": epsilon}
        ph.nac_params = nac_params

    # phonopy-params.yaml is written out.
    ph.save()
    print("phonopy_params.yaml was made for PK=%d" % pk)


if __name__ == "__main__":
    # PK as the first argument
    dump_phonopy(int(sys.argv[1]))

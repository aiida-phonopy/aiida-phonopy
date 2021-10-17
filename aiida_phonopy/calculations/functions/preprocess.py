"""Functions for phonopy pre-process."""

import numpy as np
from numpy.lib.function_base import disp
from aiida.engine import calcfunction
from aiida import orm
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import get_default_displacement_distance, get_default_physical_units


@calcfunction
def phonopy_preprocess(structure, symmetry_tolerance, displacements, supercell_matrix):
    """
    Set up the pre-process phonopy calculation. It returns the (super)cells, the primitive cell,
    the displacements dataset, the primitive_matrix and the code version. 
    """ # are the last two really necessary?

    ph = Phonopy(phonopy_atoms_from_structure(structure), 
        supercell_matrix=get_supercell_matrix(supercell_matrix),
        primitive_matrix="auto", symprec=symmetry_tolerance.value, )

    if "dataset" in displacements:
        ph.dataset = dataset.get_dict()
    else:
        kwargs = {"distance":displacements["distance"]}
        if "random" in displacements:
            rand = displacements["random"]
            for key in rand.items():
                kwargs.append({key:rand[key]})
        ph.generate_displacements(**kwargs) # phonopy instance method

    structures = get_phonon_structures(ph)
    version = ph.version # phonopy version - needed in output? Where?
    displacements_dataset = orm.Dict(dict=ph.dataset)
    primitive_matrix = orm.List(list=ph.primitive_matrix)

    return {'supercells':structures['supercells'], 
        'primitive':{'structure':structures['primitive'], 'matrix':primitive_matrix},
        'displacements_dataset':displacements_dataset}

def get_supercell_matrix(supercell_matrix):
    """Return a phonopy readible supercell matrix."""
    if len(np.ravel(supercell_matrix)) == 3:
        smat = np.diag(supercell_matrix)
    else:
        smat = np.array(supercell_matrix)

    return smat.tolist()

def phonopy_atoms_to_structure(cell):
    """Convert PhonopyAtoms to StructureData."""
    symbols = cell.symbols
    positions = cell.positions
    structure = StructureData(cell=cell.cell)
    for symbol, position in zip(symbols, positions):
        structure.append_atom(position=position, symbols=symbol)
    return structure

def phonopy_atoms_from_structure(structure):
    """Convert StructureData to PhonopyAtoms."""
    # need improvement for kind_name issue
    cell = PhonopyAtoms(
        symbols=[site.kind_name for site in structure.sites],
        positions=[site.position for site in structure.sites],
        cell=structure.cell,
    )
    return cell

def get_phonon_structures(ph):
    """Generate AiiDA structures of phonon related cells."""
    structures_dict = {'supercells':{}}

    digits = len(str(len(ph.supercells_with_displacements)))
    for i, scell in enumerate(ph.supercells_with_displacements):
        structure = phonopy_atoms_to_structure(scell)
        label = "supercell_%s" % str(i + 1).zfill(digits) # is it really necessary? Maybe simpler?
        structure.label = "%s %s" % (structure.get_formula(mode="hill_compact"), label) # is it really necessary?
        structures_dict["supercells"].update({label:structure})

    supercell_structure = phonopy_atoms_to_structure(ph.supercell)

    supercell_structure.label = "%s %s" % (supercell_structure.get_formula(mode="hill_compact"), "supercell")
    structures_dict["supercells"].update({"supercell":supercell_structure})

    primitive_structure = phonopy_atoms_to_structure(ph.primitive)
    # is the label really necessary?
    primitive_structure.label = "%s %s" % (primitive_structure.get_formula(mode="hill_compact"), "primitive cell")
    structures_dict["primitive"] = primitive_structure

    return structures_dict
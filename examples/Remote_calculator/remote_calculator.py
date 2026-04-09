from __future__ import annotations

from pathlib import Path

from aiida import load_profile
from aiida.engine import run_get_node
from aiida.orm import Dict, load_code
from aiida.plugins import DataFactory, WorkflowFactory
from ase.io import read

load_profile()

PhonopyAseWorkChain = WorkflowFactory("phonopy.ase")
StructureData = DataFactory("core.structure")

potential_fname = "nep89_20250409.txt"
HERE = Path(__file__).resolve().parent
atoms_path = HERE / "model.xyz"
potential_path = HERE / potential_fname


atoms = read(atoms_path, format="extxyz")
structure = StructureData(ase=atoms)

# LAMMPS (via ASE LAMMPSRun)
parameters = {
    "pair_style": f"matpl/nep/kk {potential_path.name}",
    "pair_coeff": ["* * Si"],
    "lammps_options": "-echo log -screen none -log /dev/stdout -k on g 1 -sf kk -pk kokkos neigh half comm device newton on",
}

files = [potential_path.name]


def make_lammps_calculator():
    """Factory executed on the remote worker to avoid pickling tmp_dir paths."""
    from ase.calculators.lammpsrun import LAMMPS

    return LAMMPS(files=files, **parameters)


inputs = PhonopyAseWorkChain.get_populated_builder(
    structure=structure,
    calculator=make_lammps_calculator,
    max_number_of_atoms=200,
    pythonjob_inputs={
        "computer": "sugon-4",
        "upload_files": {
            potential_path.name: str(potential_path),
        },
        "metadata": {
            "options": {
                "resources": {
                    "num_machines": 1,
                    "num_mpiprocs_per_machine": 1,
                },
                "queue_name": "pg_g4J4",
                "custom_scheduler_commands": "#SBATCH --gres=gpu:1",
            }
        },
    },
    phonopy_inputs={
        "code": load_code("phonopy@localhost"),
        "parameters": Dict({"band": "auto"}),
    },
)

results, node = run_get_node(PhonopyAseWorkChain, **inputs)

node = results["output_phonopy"]["phonon_bands"]
node.show_mpl()

with open("band.png",'wb') as f:
    content, _ = node._prepare_mpl_png()
    f.write(content)

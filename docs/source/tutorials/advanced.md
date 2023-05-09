---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.14.5
kernelspec:
  display_name: base
  language: python
  name: python3
---


# Advanced Tutorial

This tutorial shows some more advanced features of the Data types presented in the [previous tutorial](./intermidiate.ipynb).

## Working with custom atomic names

Many quantum codes, and AiiDA itself, allow for defining custom atomic names. This can be usefull when defining some extra features, such as magnetic ordering or on-site Hubbard values. Nevertheless, `Phonopy` is not handling such cases, and this could be a bottleneck for the usage of the package in AiiDA. In `aiida-phonopy`, we manage to overcome this problem! Let's see how does it work!

### A structure with kinds

Let's define a very simple structure that contains two atoms of the same species. Let's take cubic silicon!

```{code-cell} ipython3
from aiida import load_profile
from aiida.plugins import DataFactory

load_profile()
StructureData = DataFactory("core.structure")

a = 2.716
cell = [[0,a,a],[a,0,a],[a,a,0]]

structure = StructureData(cell=cell)
structure.append_atom(position=(a,a,a), symbols="Si", name="Si1")
structure.append_atom(position=(1.5*a,1.5*a,1.5*a), symbols="Si", name="Si2")
```

Now we can pass load it in the {py:class}`aiida_phonopy.data.preprocess.PreProcessData` with a supercell matrix of (2,2,2). Let's see if the supercell has the correct kinds.

```{code-cell} ipython3
PreProcessData = DataFactory("phonopy.preprocess")

preprocess_data =  PreProcessData(structure=structure, supercell_matrix=[2,2,2])

supercell = preprocess_data.get_supercell()
supercell.sites
```

Great! The supercell structure has the correct _kind names_.

If we have a look at the number of structure with displacements, we will notice they will be higher than the silicon structure with only chemical symbols.

```{code-cell} ipython3
supercells = preprocess_data.get_supercells_with_displacements()
print(f"Number of displacements: {len(supercells)}")
```

This may or may not be a wanted behaviour. If you want `Phonopy` to not distinguish atoms on their _name_, you can initialize the `PreProcessData` with the flag `distinguish_kinds = False`, as follows:

```{code-cell} ipython3
preprocess_data = PreProcessData(structure, supercell_matrix=[2,2,2], distinguish_kinds=False)
supercells = preprocess_data.get_supercells_with_displacements()
print(f"Now the number of displacements are: {len(supercells)}")
```

```{admonition} Exercise
:class: tip
Verify that the supercell still have the same _kind names_.
```

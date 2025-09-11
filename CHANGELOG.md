## v1.3.0

Minor release to support some breaking compatibilities with new versions of `phonopy` and `aiida-pythonjob`.
We also address in passing some deprecation warnings from `phonopy`.

### ⬆️ Update dependencies

- Drop support for phonopy versions < 2.39.
- Drop support for aiida-pythonjob versions < 0.3.


## v1.2.1

Patch release to avoid aiida-core v2.7.0 that ships with a mysterious behaviour that makes
some PhonopyCalculation to crash. Tests are also updated due to phonopy new default
calculation for force constants that makes use of an other package that is still under active developments and does not manage to reproduce/interpolate the old default behaviour.

### ⬆️ Update dependencies

- Do not support aiida-core v2.7.0 due to mysterious bug that makes PhonopyCalculation to crash.



# Changelog

## v1.2.0 - 2024-01-27

[full changelog](https://github.com/aiida-phonopy/aiida-phonopy/compare/v1.2.0...v1.1.4)

### ✨ New features

- `PhonopyAseWorkChain`: automated phonons using _any_ ASE calculator (ideal for ML potentials).

### 🗑️ Deprecations

- Deprecate python v3.8 and lower.

### 📚 Documentation

- Added tutorials

### ⬆️ Update dependencies

- Add `aiida-pythonjob` as dependency to run ASE on remote computers (e.g. ML calculators).
- Update python dependencies up to v3.12.


## v1.1.4 - 2023-12-09

[full changelog](https://github.com/aiida-phonopy/aiida-phonopy/compare/v1.1.4...v1.1.3)

Minor release that solves a small bug, improves documentation and typing.

- Fixes:
  - `RawData` was not settings the `primitive_matrix` when `is_symmetry = False`, since `Phonopy` doesn't specify it in this case; the solution is to set it to `diag(1,1,1)` when the symmetry is turned off (https://github.com/aiida-phonopy/aiida-phonopy/commit/6ef9e8e10e5d6ea6a27f75b4c38b255f3c8e9d84)
  - Typo in the "Get Started" documentation (https://github.com/aiida-phonopy/aiida-phonopy/commit/57327f534a32bb1beda8596711d148f50f68e1db)

- Docs:
  - Improvements to the "Get Started" page (https://github.com/aiida-phonopy/aiida-phonopy/commit/57327f534a32bb1beda8596711d148f50f68e1db)

- Tests:
  - Add more tests for the data types (https://github.com/aiida-phonopy/aiida-phonopy/commit/6ef9e8e10e5d6ea6a27f75b4c38b255f3c8e9d84)

- New contributors:
  - @mbercx made his first contribution to the repo :tada:

## v1.1.3 - 2023-05-24

[full changelog](https://github.com/aiida-phonopy/aiida-phonopy/compare/v1.1.3...v1.1.2)

Minor release that solves a small bug.

- Fixes:
  - `generate_preprocess_data` had a typo in the typing which made validation to crash (https://github.com/aiida-phonopy/aiida-phonopy/commit/d2c32cf25ebe469797feab46d7e50dde4e7baba3)

## v1.1.2 - 2023-05-22

[full changelog](https://github.com/aiida-phonopy/aiida-phonopy/compare/v1.1.2...v1.1.1)

Minor release that solves a small bug.

- Fixes:
  - `PhonopyCalculation` crashed when `force_constants` was used as input (https://github.com/aiida-phonopy/aiida-phonopy/commit/449447aa19a495edc6a28b0333a5e6d8e1f868bb)

## v1.1.1 - 2023-05-10

[full changelog](https://github.com/aiida-phonopy/aiida-phonopy/compare/v1.1.1...v1.1.0)

This release comes with some bug fixes of the previous version for data storage,
it improves the docstrings and typing, and finally adds the documentation.

## v1.1.0 - 2023-05-02

[full changelog](https://github.com/aiida-phonopy/aiida-phonopy/compare/v1.1.0...v1.0.0)

In this release we apply some improvements in the way the data is stored, and refactor the calculation utils to make names shorter and more comprehensible.

## v1.0.0 - 2022-11-0X

[full changelog](https://github.com/aiida-phonopy/aiida-phonopy/compare/v1.0.0...v0.7.0)

In this release the new version of AiiDA v2.x is fully supported.


## v0.7.0 - 2022-11-0X

[full changelog](https://github.com/aiida-phonopy/aiida-phonopy/compare/v0.7.0...v0.6.0)

This new release contains radical changes from the previous versions. The new changes have been made to encounter the request and the possibility
to interface with more (possibly any) force calculator.

### New design

The design of this new version was thought in terms of ``DataTypes``. The new data types will allow to load
the StructureData and save the information regarding the displacements dataset. This allows for an easy, direct
and fast access to the supercells on top of which one might want to calculate the forces.
The design of these data types allows for a consistent load of information, from the structure to the forces on supercells.
A specific workflows which run all the structures is thus no more necessary, but only convenient. An example of workflows is still
provided in the `workflows`. Since a workflow of this type can can be optimized for the specific force calculator (plugin),
we are not here addressing to the request of a general workflos. On the contrary, we have designed a general interface.
More technically, we have "wrapped" the ``phonopy`` code itself in a way that it can interface with aiida in the
more easy and accessible way.

### New Data

We have designed new data types. We can distinguished two main data types: the *pre-processing* data types, and the *post-processing* data types.
As the names suggest, the formers one is meant to store the pre-process information, e.g. the displacement dataset, structure information, ... .
The latters will wrap all or part of the pre-process information, and is meant to store the forces/force constants information.
The post-processing data types can then be used to get information through the **new CalcJob** that we designed.
The new DataTypes are:
- **RawData**: this is meant to be a base class for the other ones; it already contains the main benefits of the new design.
  With this data one can store the information of a StructureData in a way that we can convert it to phonopy format.
  Structures with *kind names* are handled in a way that can be read and distinguished when interfacing with phonopy.
  It contains many methods, among which a set of **calcfunctions** utilities, accessible via the ``.calcfunctions`` namespace.
- **PreProcessData**: based on RawData, allows to set and store the information regarding the displacement dataset.
  It contains many methods. From the ``.calcfunctions`` namespace one can store and get the supercells on which to calculate the forces,
  and many other information.
- **PhonopyData**: probably the main DataTypes, it is base on top of the PreProcessData, and here one can store the information
  regarding the forces related to the displacements set (a cross verification is made). This can be used as the main data for
  getting phonon properties.
- **ForceConstantsData**: base on the RawData, this is meant to store force constants data computed using other methods
  (e.g. linear responce) and structure information, in order to be able to exploit the post processing capabilities of phonopy
  and our new phonopy CalcJob.

### New CalcJob

The new CalcJob exploits all the phonopy post processing capabilities. In respect to the older versions, now
any post processing property can be calculated and it is **parsed** in the appropriate format. The hdf5 format
is always used for having fast calculation and fast retrieval of files. Moreover, the inputs/outputs interface
is now very user-friendly!

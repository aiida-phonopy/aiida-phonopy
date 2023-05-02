# Changelog

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

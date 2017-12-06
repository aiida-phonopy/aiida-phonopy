Gruneisen
=========

This WorkChain performs a mode Gruneisen parameters calculation using Phonopy. This WorkChain is designed to
keep the compatibility with phonon WorkChain input structure. For this reason the input parameters in the WorkChain
have the same structure. Please, check phonon Workchain documentation for detailed information.
This WorkChain calculates the mode Gruneisen parameters by performing 3 different phonon calculations.
One with the crystal structure optimized at a given external stress (default: 0 kB) and the other two
with the unit cell optimized with a slightly higher and lower stress (defined by stress_displacement)
obtaining a slightly smaller and larger unit cell respectively.
Stress_displacement can be set as an optional argument, by default its value is 1e-2 kB.

.. function:: GruneisenPhonopy(structure, machine, ph_settings, es_settings [, stress_displacement=1e-2])

   :param structure: AiiDA StructureData object that contains the crystal unit cell structure.
   :param ph_settings: AiiDA ParametersData data  object that contains the phonopy input parameters.
   :param es_settings: AiiDA ParameterData object that contains the calculator input parameters.
   :param pressure (optional): AiiDA FloatData object. This determines the absolute stress (in kBar) at which the reference crystal structure is optimized (default 0).
   :param stress_displacement (optional): AiiDA FloatData object. This determines the stress difference between the 3 phonon calculations (default 1e-2 kB).

The outputs of this WorkChain are:

* **band_structure**: BandStructure object that contains the phonon band structure and the mode Gruneisen parameters.
* **mesh**: ArrayData object that contains the wave vectors sampling mesh and the mode Gruneises parameters at each wave vector.

Each one of this objects has its own methods for extracting the information. Check the individual object documentation
for more details. **workchains/tools/plot_gruneisen.py** contains a complete example script showing how to extract the information from these outputs.

example of use
--------------
::

    GruneisenPhonopy = WorkflowFactory('phonopy.gruneisen')

Run in interactive ::

    result = run(GruneisenPhonopy,
                 structure=structure,
                 es_settings=es_settings,
                 ph_settings=ph_settings,
                 Optional settings
                 pressure=Float(0),
                 stress_displacement=Float(1e-2),
                 )


Run by the deamon ::

    future = submit(GruneisenPhonopy,
                    structure=structure,
                    es_settings=es_settings,
                    ph_settings=ph_settings,
                    pressure=Float(0),
                    stress_displacement=Float(1e-2)
                    )

    print('Running WorkChain with pk={}'.format(future.pid))

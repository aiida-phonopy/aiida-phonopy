Phonon
======

This WorkChain performs a phonon calculation using Phonopy. This WorkChain requires at least one of the following
AiiDA plugins to work: QuantumESPRESSO (https://github.com/aiidateam/aiida-quantumespresso),
VASP (https://github.com/DropD/aiida-vasp) or LAMMPS (https://github.com/abelcarreras/aiida-lammps).
At the present time Born effective charges are only calculated using VASP plugin.
Phonopy can be used either locally and remotely. To use phonopy remotely phonopy code must be setup as described
in AiiDA documentation (https://aiida-core.readthedocs.io/en/latest/get_started/index.html#code-setup-and-configuration).
using the phonopy plugin provided in this package.

.. function:: PhononPhonopy(structure, ph_settings, es_settings [, optimize=True, pressure=0.0, use_nac=False])

   :param structure: StructureData object that contains the crystal unit cell structure.
   :param ph_settings: ParametersData data  object that contains the phonopy input parameters.
   :param es_settings: ParameterData object that contains the calculator input parameters. These parameters depends on the code used (see workchains/launcher examples)
   :param optimize: (optional) BooleanData object. Determines if a crystal unit cell optimization is performed or not before the phonon calculation. By default this option is True.
   :param pressure: (optional) FloatData object. If optimize is True, this sets the external pressure (in kB) at which the unit cell optimization is preformed. By default this option takes value 0 kB.
   :param use_nac: (optional) BooleanData object. Determines if non-analytical corrections will be included in the phonon calculation. By default this option is False.

- ph_settings: This object contains a dictionary with all input parameters for phonopy. See plugins section for more information. Additional dictionary entries can be added to request a remote phonopy calculation. See example in examples/workchains/launh_phonon_gan ::

    code: phonopy@cluster
    machine: machine_dict

machine_dict dictionary should contain the following entries. resources_dict may change depending on the scheduler ::

    machine_dict = {'resources': resources_dict
                    'max_wallclock_seconds': 3600 * 10 # in seconds
                    }

    resources_dic = {'num_machines': 1,
                     'parallel_env': 'mpi*',
                     'tot_num_mpiprocs': 16
                     }

- es_settings: This object contains the parameters for each specific software used as calculator (QE, VASP or LAMMPS). Each calculator uses a different dictionary structure (See examples at examples/workchains folder for the details). The common basic structure in all calculators is ::

    settings_dict = {'code': {'optimize': 'vasp5@boston',
                              'forces': 'vasp4@boston',
                              'born': 'vasp4@boston'},

                     'parameters': parameters, # this depends on calculator (see examples)
                     'machine': machine_dict, # same dictionary defined above
                     ...
                     }

    es_settings = ParameterData(dict=settings_dict)

If the code used in all calculations types (optimize, forces and born) is the same, the dictionary can be written as ::

    settings_dict = {'code': vasp@boston',
                     'parameters': parameters,
                     'machine': machine_dict,
                     ...
                     }


The results outputs of this WorkChain are the following :

* **force_constants**: ForceConstantsData object that contains the harmonic force constants. If Born effective charges are calculated this object also contains the dielectric tensor and the born effective charges of each atom in the unit cell
* **thermal_properties**: ArrayData object that contains the thermal properties calculated using phonopy. These include the entropy, free energy and heat capacity at constant volume.
* **dos**: PhononDosData object that contains the phonon full and partial density of states.
* **band_structure**: BandStructureData object that contains the harmonic phonon band structure.
* **final_structure**: StructureData object that contains the optimized unit cell. If no optimization is performed this is the same StructureData object provided as a input.

Each one of this objects has its own methods for extracting the information. Check the individual object documentation
for details. **workchains/tools/plot_phonon.py** contains a complete example script showing how to extract the information from these outputs.

example of use
--------------
::

    PhononPhonopy = WorkflowFactory('phonopy.phonon')

Run in interactive ::

    result = run(PhononPhonopy,
                 structure=structure,
                 es_settings=es_settings,
                 ph_settings=ph_settings,
                 pressure=Float(0.0),
                 optimize=Bool(True),
                 use_nac=Bool(True)
                 )

Run by the deamon ::

    future = submit(PhononPhonopy,
                    structure=structure,
                    es_settings=es_settings,
                    ph_settings=ph_settings,
                    pressure=Float(0),
                    optimize=Bool(True),
                    use_nac=Bool(True)
                    )
    print('Running WorkChain with pk={}'.format(future.pid))

from setuptools import setup, find_packages

setup(
    name='aiida-phonopy',
    version='0.1',
    description='AiiDA plugin for running phonon calculations using phonopy',
    author='Abel Carreras',
    author_email='abelcarreras83@gmail.com',
    license='MIT license',
    packages=find_packages(exclude=['aiida']),
    requires=['phonopy', 'numpy'],
    setup_requires=['reentry'],
    reentry_register=True,
    include_package_data=True,
    entry_points={
        'aiida.data': [
            'band_structure = aiida_phonopy.data.band_structure: BandStructureData',
            'force_constants = aiida_phonopy.data.force_constants: ForceConstantsData',
            'force_sets = aiida_phonopy.data.force_sets: ForceSetsData',
            'phonon_dos = aiida_phonopy.data.phonon_dos: PhononDosData'
        ],
        'aiida.calculations': [
            'phonopy = aiida_phonopy.plugins.jobs.phonopy: PhonopyCalculation'
        ],
        'aiida.parsers': [
            'phonopy = aiida_phonopy.plugins.parsers.phonopy: PhonopyParser'
        ],
        'aiida.workflows': [
             'phonopy.optimize = aiida_phonopy.workchains.optimize: OptimizeStructure',
             'phonopy.phonon = aiida_phonopy.workchains.phonon: PhononPhonopy',
             'phonopy.gruneisen = aiida_phonopy.workchains.gruneisen: GruneisenPhonopy'
         #    'phonopy.optimize = workchains.wf_optimize: OptimizeStructure',
         #   'phonopy.phonon = workchains.wf_phonon: PhononPhonopy',
         #   'phonopy.gruneisen = workchains.wf_gruneisen: GruneisenPhonopy',
        ]
    }
)
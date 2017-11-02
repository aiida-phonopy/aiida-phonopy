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
    entry_points={
        'aiida.calculations': [
            'phonopy = plugins.jobs.phonopy: PhonopyCalculation'
        ],
        'aiida.data': [
            'band_structure = data.band_structure: BandStructureData'
            'force_constants = data.force_constants: ForceConstantsData'
            'force_sets = data.force_sets: ForceSetsData'
            'phonon_dos = data.phonon_dos: PhononDosData'
        ],
        'aiida.parsers': [
            'phonopy = plugins.parsers.phonopy: PhonopyParser'
#        ],
#        'aiida.workflows': [
#            'wc_phonon = workflows.wf_phonon:PhononPhonopy',
#            'wc_gruneisen = workflows.wf_gruneisen:GruneisenPhonopy',
#            'generate_inputs = workflows.generate_inputs:generate_inputs',
#            'parse_optimize_calculation = workflows.parse_interface:parse_optimize_calculation',
        ]
    }
)
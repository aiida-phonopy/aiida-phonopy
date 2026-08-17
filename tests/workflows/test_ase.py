"""Test for the :mod:`~aiida_phonopy.workflows.ase` module."""

import shutil
import sys

import pytest


@pytest.fixture
def generate_workchain_phonopy_ase(fixture_localhost, fixture_code, generate_workchain, generate_structure):
    """Generate an instance of a `PhonopyAseWorkChain`."""

    def _generate_workchain_phonopy_ase(append_inputs=None, phonon_inputs=None, return_inputs=False, calculator=None):
        from aiida.orm import Dict, InstalledCode
        import ase.calculators.lj

        from aiida_phonopy.workflows.ase import PhonopyAseWorkChain

        entry_point = 'phonopy.ase'

        if not shutil.which('phonopy'):
            phonopy_inputs = None
        else:
            code = InstalledCode(
                computer=fixture_localhost,
                filepath_executable=shutil.which('phonopy'),
                input_plugin_name='phonopy.phonopy',
                label='phonopy',
                default_calc_job_plugin='phonopy.phonopy',
            )
            phonopy_inputs = {
                'code': code,
                'parameters': Dict({'band': 'auto'}),
            }

        if calculator is None:
            calculator = ase.calculators.lj.LennardJones()

        inputs = PhonopyAseWorkChain.get_populated_builder(
            structure=generate_structure(),
            calculator=calculator,
            max_number_of_atoms=40,
            pythonjob_inputs={
                'computer': fixture_localhost.hostname,
                # Run the job with the same interpreter as the test session: the pickled
                # function must be unpickled by a compatible python having e.g. `cloudpickle`,
                # which a bare `python3` resolved through the computer shell may not provide.
                'command_info': {'filepath_executable': sys.executable},
            },
            phonopy_inputs=phonopy_inputs,
        )

        process = generate_workchain(entry_point, inputs)

        return process

    return _generate_workchain_phonopy_ase


def test_validation(generate_workchain_phonopy_ase):
    """Test if the validation and the population of inputs works."""
    generate_workchain_phonopy_ase()


def test_run(generate_workchain_phonopy_ase):
    """Test if the validation and the population of inputs works.

    ..note:: this is really running the workchain, so it might take a while.
    """
    from aiida.engine import run_get_node

    results, node = run_get_node(generate_workchain_phonopy_ase())

    print('*' * 50)
    print('*' * 50)
    print(node.called[-1].outputs.retrieved.get_object_content('aiida.out'))
    print('*' * 50)
    print('*' * 50)
    assert node.is_finished_ok

    # phonopy_data = results['phonopy_data'].get_phonopy_instance()
    # phonopy_data.produce_force_constants()

    # assert 'phonon_bands' in results['output_phonopy']


def test_run_with_calculator_factory(generate_workchain_phonopy_ase):
    """Run workchain with calculator as a factory callable (for remote PythonJob / tmp-dir calculators).

    ..note:: this is really running the workchain, so it might take a while.
    """
    from aiida.engine import run_get_node
    import ase.calculators.lj

    results, node = run_get_node(generate_workchain_phonopy_ase(calculator=lambda: ase.calculators.lj.LennardJones()))
    assert node.is_finished_ok

    phonopy_data = results['phonopy_data'].get_phonopy_instance()
    phonopy_data.produce_force_constants()

    assert 'phonon_bands' in results['output_phonopy']

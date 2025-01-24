# -*- coding: utf-8 -*-
"""Test for the :mod:`~aiida_phonopy.workflows.ase` module."""
import pytest


@pytest.fixture
def generate_workchain_phonopy_ase(fixture_localhost, generate_workchain, generate_structure):
    """Generate an instance of a `PhonopyAseWorkChain`."""

    def _generate_workchain_phonopy_ase(append_inputs=None, phonon_inputs=None, return_inputs=False):
        from ase.calculators.lj import LennardJones

        from aiida_phonopy.workflows.ase import PhonopyAseWorkChain

        entry_point = 'phonopy.ase'

        inputs = PhonopyAseWorkChain.get_populated_builder(
            structure=generate_structure(),
            calculator=LennardJones(),
            max_number_of_atoms=40,
            pythonjob_inputs={'computer': fixture_localhost.hostname}
        )

        process = generate_workchain(entry_point, inputs)

        return process

    return _generate_workchain_phonopy_ase


def test_validation(generate_workchain_phonopy_ase):
    """Test if the validation and the population of inputs works."""
    generate_workchain_phonopy_ase()


def test_run(generate_workchain_phonopy_ase):
    """Test if the validation and the population of inputs works."""
    from aiida.engine import run_get_node

    results, node = run_get_node(generate_workchain_phonopy_ase())
    assert node.is_finished_ok

    phonopy_data = results['phonopy_data'].get_phonopy_instance()
    phonopy_data.produce_force_constants()

# -*- coding: utf-8 -*-
"""Test for the :mod:`~aiida_phonopy.workflows.ase` module."""
import shutil

import pytest


@pytest.fixture
def generate_workchain_phonopy_ase(fixture_localhost, fixture_code, generate_workchain, generate_structure):
    """Generate an instance of a `PhonopyAseWorkChain`."""

    def _generate_workchain_phonopy_ase(append_inputs=None, phonon_inputs=None, return_inputs=False):
        from aiida.orm import Dict, InstalledCode
        from ase.calculators.lj import LennardJones

        from aiida_phonopy.workflows.ase import PhonopyAseWorkChain

        entry_point = 'phonopy.ase'

        if not shutil.which('phonopy'):
            phonopy_inputs = None
        else:
            code = InstalledCode(
                label='phonopy',
                computer=fixture_localhost,
                filepath_executable=shutil.which('phonopy'),
                default_calc_job_plugin='phonopy.phonopy',
            )
            phonopy_inputs = {
                'code': code,
                'parameters': Dict({'band': 'auto'}),
            }

        inputs = PhonopyAseWorkChain.get_populated_builder(
            structure=generate_structure(),
            calculator=LennardJones(),
            max_number_of_atoms=40,
            pythonjob_inputs={'computer': fixture_localhost.hostname},
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

    .. note:: this is really running the workchain, so it might take a while.
    """
    from aiida.engine import run_get_node

    results, node = run_get_node(generate_workchain_phonopy_ase())
    assert node.is_finished_ok

    phonopy_data = results['phonopy_data'].get_phonopy_instance()
    phonopy_data.produce_force_constants()

    assert 'phonon_bands' in results['output_phonopy']

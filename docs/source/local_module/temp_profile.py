# -*- coding: utf-8 -*-
"""Load and populate a temporary profile with a computer and code."""
from __future__ import annotations

from dataclasses import dataclass
import os
import pathlib
import shutil

from aiida import get_profile, load_ipython_extension, load_profile, manage, orm
from aiida.storage.sqlite_temp import SqliteTempBackend
import psutil


@dataclass
class AiiDALoaded:
    profile: manage.Profile
    computer: orm.Computer | None
    code: orm.Code | None
    structure: orm.StructureData | None
    cpu_count: int
    workdir: pathlib.Path
    phonopy_path: pathlib.Path


def load_temp_profile(
    name='temp_profile',
    add_computer=False,
    add_phonopy_code=False,
    add_structure_si=False,
    debug=False,
    wipe_previous=True,
    cpu_count: int | None = None,
):
    """Load a temporary profile with a computer and code.

    This function is idempotent, so it can be called multiple times without
    creating duplicate computers and codes.

    :param name: The name of the profile to load.
    :param add_computer: Whether to add a computer to the profile.
    :param add_phonopy_code: Whether to add a phonopy code to the profile.
    :param debug: Whether to enable debug mode (printing all SQL queries).
    :param wipe_previous: Whether to wipe any previous data
    """
    # load the ipython extension, if possible
    try:
        load_ipython_extension(get_ipython())
    except NameError:
        pass

    workdir_path = pathlib.Path(__file__).parent / '_aiida_workdir' / name
    repo_path = pathlib.Path(os.environ['AIIDA_PATH']) / '.aiida' / 'repository' / name

    profile = get_profile()

    if not (profile and profile.name == name):

        if wipe_previous and repo_path.exists():
            shutil.rmtree(repo_path)
        if wipe_previous and workdir_path.exists():
            shutil.rmtree(workdir_path)

        profile = SqliteTempBackend.create_profile(
            name,
            options={'runner.poll.interval': 1},
            debug=debug,
        )
        load_profile(profile, allow_switch=True)
        config = manage.get_config()
        config.add_profile(profile)

    cpu_count = cpu_count or min(2, psutil.cpu_count(logical=False))
    if not shutil.which('phonopy'):
        raise RuntimeError('phonopy not found in PATH')
    phonopy_path = pathlib.Path(shutil.which('phonopy'))

    computer = load_computer(workdir_path, cpu_count) if add_computer else None
    phonopy_code = load_phonopy_code(computer, phonopy_path) if (computer and add_phonopy_code) else None
    structure = create_si_structure() if add_structure_si else None

    return AiiDALoaded(
        profile,
        computer,
        phonopy_code,
        structure,
        cpu_count,
        workdir_path,
        phonopy_path,
    )


def load_computer(work_directory: pathlib.Path, cpu_count: int):
    """Idempotent function to add the computer to the database."""
    created, computer = orm.Computer.collection.get_or_create(
        label='local_direct',
        description='local computer with direct scheduler',
        hostname='localhost',
        workdir=str(work_directory.absolute()),
        transport_type='core.local',
        scheduler_type='core.direct',
    )
    if created:
        computer.store()
        computer.set_minimum_job_poll_interval(0.0)
        computer.set_default_mpiprocs_per_machine(cpu_count)
        computer.configure()
    return computer


def load_phonopy_code(computer, exec_path: pathlib.Path):
    """Idempotent function to add the code to the database."""
    try:
        code = orm.load_code('phonopy@local_direct')
    except:
        code = orm.Code(
            input_plugin_name='phonopy.phonopy',
            remote_computer_exec=[computer, str(exec_path)],
        )
        code.label = 'phonopy'
        code.description = 'phonopy code on local computer'
        code.store()
    return code


def create_si_structure():
    cell = [
        [3.7881476451529, 0.0, 0.0],
        [1.8940738225764, 3.2806320939886, 0.0],
        [1.8940738225764, 1.0935440313296, 3.0930096003167],
    ]
    structure = orm.StructureData(cell=cell)
    structure.append_atom(position=(0.0, 0.0, 0.0), symbols='Si')
    structure.append_atom(position=(1.8940738225764, 1.0935440313296, 0.77325240007918), symbols='Si')
    structure.store()
    return structure

# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name,too-many-statements
"""Initialise a text database and profile for pytest."""
from collections.abc import Mapping
import os
import shutil

import pytest

pytest_plugins = ['aiida.manage.tests.pytest_fixtures']  # pylint: disable=invalid-name


@pytest.fixture(scope='session')
def filepath_tests():
    """Return the absolute filepath of the `tests` folder.

    .. warning:: if this file moves with respect to the `tests` folder, the implementation should change.

    :return: absolute filepath of `tests` folder which is the basepath for all test resources.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def filepath_fixtures(filepath_tests):
    """Return the absolute filepath to the directory containing the file `fixtures`."""
    return os.path.join(filepath_tests, 'fixtures')


@pytest.fixture(scope='function')
def fixture_sandbox():
    """Return a `SandboxFolder`."""
    from aiida.common.folders import SandboxFolder

    with SandboxFolder() as folder:
        yield folder


@pytest.fixture
def fixture_localhost(aiida_localhost):
    """Return a localhost `Computer`."""
    localhost = aiida_localhost
    localhost.set_default_mpiprocs_per_machine(1)
    return localhost


@pytest.fixture
def fixture_code(fixture_localhost):
    """Return a `Code` instance configured to run calculations of given entry point on localhost `Computer`."""

    def _fixture_code(entry_point_name):
        from aiida.common import exceptions
        from aiida.orm import InstalledCode, load_code

        label = f'test.{entry_point_name}'

        try:
            return load_code(label)
        except exceptions.NotExistent:
            return InstalledCode(
                label=label,
                computer=fixture_localhost,
                filepath_executable='/bin/true',
                default_calc_job_plugin=entry_point_name,
            )

    return _fixture_code


@pytest.fixture
def serialize_builder():
    """Serialize the given process builder into a dictionary with nodes turned into their value representation.

    :param builder: the process builder to serialize
    :return: dictionary
    """

    def serialize_data(data):
        # pylint: disable=too-many-return-statements
        from aiida.orm import AbstractCode, BaseType, Data, Dict, KpointsData, RemoteData
        from aiida.plugins import DataFactory

        StructureData = DataFactory('core.structure')
        UpfData = DataFactory('pseudo.upf')

        if isinstance(data, dict):
            return {key: serialize_data(value) for key, value in data.items()}

        if isinstance(data, BaseType):
            return data.value

        if isinstance(data, AbstractCode):
            return data.full_label

        if isinstance(data, Dict):
            return data.get_dict()

        if isinstance(data, StructureData):
            return data.get_formula()

        if isinstance(data, UpfData):
            return f'{data.element}<md5={data.md5}>'

        if isinstance(data, RemoteData):
            # For `RemoteData` we compute the hash of the repository. The value returned by `Node._get_hash` is not
            # useful since it includes the hash of the absolute filepath and the computer UUID which vary between tests
            return data.base.repository.hash()

        if isinstance(data, KpointsData):
            try:
                return data.get_kpoints()
            except AttributeError:
                return data.get_kpoints_mesh()

        if isinstance(data, Data):
            return data.base.caching._get_hash()  # pylint: disable=protected-access

        return data

    def _serialize_builder(builder):
        return serialize_data(builder._inputs(prune=True))  # pylint: disable=protected-access

    return _serialize_builder


@pytest.fixture
def generate_calc_job():
    """Fixture to construct a new `CalcJob` instance and call `prepare_for_submission` for testing `CalcJob` classes.

    The fixture will return the `CalcInfo` returned by `prepare_for_submission` and the temporary folder that was passed
    to it, into which the raw input files will have been written.
    """

    def _generate_calc_job(folder, entry_point_name, inputs=None):
        """Fixture to generate a mock `CalcInfo` for testing calculation jobs."""
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager
        from aiida.plugins import CalculationFactory

        manager = get_manager()
        runner = manager.get_runner()

        process_class = CalculationFactory(entry_point_name)
        process = instantiate_process(runner, process_class, **inputs)

        calc_info = process.prepare_for_submission(folder)

        return calc_info

    return _generate_calc_job


@pytest.fixture
def generate_calc_job_node(fixture_localhost):
    """Fixture to generate a mock `CalcJobNode` for testing parsers."""

    def flatten_inputs(inputs, prefix=''):
        """Flatten inputs recursively like :meth:`aiida.engine.processes.process::Process._flatten_inputs`."""
        flat_inputs = []
        for key, value in inputs.items():
            if isinstance(value, Mapping):
                flat_inputs.extend(flatten_inputs(value, prefix=prefix + key + '__'))
            else:
                flat_inputs.append((prefix + key, value))
        return flat_inputs

    def _generate_calc_job_node(
        entry_point_name='phonopy.phonopy',
        computer=None,
        test_name=None,
        inputs=None,
        attributes=None,
        retrieve_temporary=None,
    ):
        """Fixture to generate a mock `CalcJobNode` for testing parsers.

        :param entry_point_name: entry point name of the calculation class
        :param computer: a `Computer` instance
        :param test_name: relative path of directory with test output files in the `fixtures/{entry_point_name}` folder.
        :param inputs: any optional nodes to add as input links to the corrent CalcJobNode
        :param attributes: any optional attributes to set on the node
        :param retrieve_temporary: optional tuple of an absolute filepath of a temporary directory and a list of
            filenames that should be written to this directory, which will serve as the `retrieved_temporary_folder`.
            For now this only works with top-level files and does not support files nested in directories.
        :return: `CalcJobNode` instance with an attached `FolderData` as the `retrieved` node.
        """
        from aiida import orm
        from aiida.common import LinkType
        from aiida.plugins.entry_point import format_entry_point_string

        if computer is None:
            computer = fixture_localhost

        filepath_folder = None

        if test_name is not None:
            basepath = os.path.dirname(os.path.abspath(__file__))
            filename = os.path.join(entry_point_name[len('phonopy.'):], test_name)
            filepath_folder = os.path.join(basepath, 'parsers', 'fixtures', filename)
            # filepath_input = os.path.join(filepath_folder, "aiida.in")

        entry_point = format_entry_point_string('aiida.calculations', entry_point_name)

        node = orm.CalcJobNode(computer=computer, process_type=entry_point)
        node.base.attributes.set('input_filename', 'aiida.in')
        node.base.attributes.set('output_filename', 'aiida.out')
        node.base.attributes.set('error_filename', 'aiida.err')
        node.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        node.set_option('max_wallclock_seconds', 1800)

        if attributes:
            node.base.attributes.set_many(attributes)  # here if you would specify temp folder etc

        if inputs:
            metadata = inputs.pop('metadata', {})
            options = metadata.get('options', {})

            for name, option in options.items():
                node.set_option(name, option)

            # here we link the inputs (which determines what to parse) to the CalcJobNode
            for link_label, input_node in flatten_inputs(inputs):
                input_node.store()
                node.base.links.add_incoming(input_node, link_type=LinkType.INPUT_CALC, link_label=link_label)

        node.store()

        if retrieve_temporary:
            dirpath, filenames = retrieve_temporary
            for filename in filenames:
                try:
                    shutil.copy(os.path.join(filepath_folder, filename), os.path.join(dirpath, filename))
                except FileNotFoundError:
                    pass  # To test the absence of files in the retrieve_temporary folder

        if filepath_folder:
            retrieved = orm.FolderData()
            retrieved.base.repository.put_object_from_tree(filepath_folder)

            # Remove files that are supposed to be only present in the retrieved temporary folder
            if retrieve_temporary:
                for filename in filenames:
                    try:
                        retrieved.base.repository.delete_object(filename)
                    except OSError:
                        pass  # To test the absence of files in the retrieve_temporary folder

            retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
            retrieved.store()

            remote_folder = orm.RemoteData(computer=computer, remote_path='/tmp')
            remote_folder.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='remote_folder')
            remote_folder.store()

        return node

    return _generate_calc_job_node


@pytest.fixture(scope='session')
def generate_parser():
    """Fixture to load a parser class for testing parsers."""

    def _generate_parser(entry_point_name):
        """Fixture to load a parser class for testing parsers.

        :param entry_point_name: entry point name of the parser class
        :return: the `Parser` sub class
        """
        from aiida.plugins import ParserFactory

        return ParserFactory(entry_point_name)

    return _generate_parser


@pytest.fixture
def generate_structure():
    """Return a `StructureData` representing bulk silicon."""

    def _generate_structure(structure_id='silicon'):
        """Return a `StructureData` representing bulk silicon."""
        from aiida.orm import StructureData

        param = 5.43
        cell = [[0.0, param / 2.0, param / 2.0], [param / 2.0, 0, param / 2.0], [param / 2.0, param / 2.0, 0.0]]
        structure = StructureData(cell=cell)

        if structure_id == 'silicon':
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols='Si', name='Si')
            structure.append_atom(position=(param / 4.0, param / 4.0, param / 4.0), symbols='Si', name='Si')
        elif structure_id == 'silicon-with-names':
            structure.append_atom(position=(0.0, 0.0, 0.0), symbols='Si', name='A')
            structure.append_atom(position=(param / 4.0, param / 4.0, param / 4.0), symbols='Si', name='B')
        else:
            raise KeyError(f"Unknown structure_id='{structure_id}'")
        return structure

    return _generate_structure


@pytest.fixture
def generate_raw_data(generate_structure):
    """Return a `RawData`."""

    def _generate_raw_data(structure_id='silicon', inputs=None):
        """Return a `RawData`."""
        from aiida_phonopy.data.raw import RawData
        structure = generate_structure(structure_id=structure_id)

        if inputs is None:
            inputs = {}

        try:
            structure = inputs.pop('structure')
        except (KeyError, AttributeError):
            pass

        raw_data = RawData(structure=structure, **inputs)

        return raw_data

    return _generate_raw_data


@pytest.fixture
def generate_preprocess_data(generate_structure):
    """Return a `PreProcessData`."""

    def _generate_preprocess_data(structure_id='silicon', inputs=None):
        """Return a `PreProcessData`."""
        from aiida_phonopy.data.preprocess import PreProcessData
        structure = generate_structure(structure_id=structure_id)

        if inputs is None:
            inputs = {}

        try:
            structure = inputs.pop('structure')
        except (KeyError, AttributeError):
            pass

        preprocess_data = PreProcessData(structure=structure, **inputs)

        return preprocess_data

    return _generate_preprocess_data


@pytest.fixture
def generate_phonopy_data(generate_preprocess_data):
    """Return a `PhonopyData`."""

    def _generate_phonopy_data(preprocess_data=None, forces=None, dielectric=None, born_charges=None):
        """Return a `PhonopyData`."""
        from aiida_phonopy.data.phonopy import PhonopyData

        if preprocess_data is None:
            preprocess_data = generate_preprocess_data()
            preprocess_data.set_displacements()

        phonopy_data = PhonopyData(preprocess_data=preprocess_data)

        if forces is None:
            phonopy_data.set_forces([[[1., 0., 0.], [-1., 0., 0.]]  # 1st displacement
                                     ])
        else:
            phonopy_data.set_forces(forces)

        if dielectric is not None:
            phonopy_data.set_dielectric(dielectric)

        if born_charges is not None:
            phonopy_data.set_born_charges(born_charges)

        return phonopy_data

    return _generate_phonopy_data


@pytest.fixture
def generate_example_phonopy_data():
    """Return BTO PhonopyData."""

    def _generate_example_phonopy_data():
        """Return BTO PhonopyData."""
        import os

        import phonopy

        from aiida_phonopy.data import PhonopyData, PreProcessData

        basepath = os.path.dirname(os.path.abspath(__file__))
        phyaml = os.path.join(basepath, 'fixtures', 'phonopy', 'phonopy.yaml')
        fsets = os.path.join(basepath, 'fixtures', 'phonopy', 'FORCE_SETS')

        ph = phonopy.load(phyaml, force_sets_filename=fsets)

        preprocess_data = PreProcessData(phonopy_atoms=ph.unitcell, supercell_matrix=[3, 3, 3])
        preprocess_data.set_displacements_from_dataset(dataset=ph.dataset)

        phonopy_data = PhonopyData(preprocess_data=preprocess_data)
        phonopy_data.set_forces(sets_of_forces=ph.forces)

        return phonopy_data

    return _generate_example_phonopy_data


@pytest.fixture
def generate_remote_data():
    """Return a `RemoteData` node."""

    def _generate_remote_data(computer, remote_path, entry_point_name=None):
        """Return a `RemoteData` node."""
        from aiida.common.links import LinkType
        from aiida.orm import CalcJobNode, RemoteData
        from aiida.plugins.entry_point import format_entry_point_string

        entry_point = format_entry_point_string('aiida.calculations', entry_point_name)

        remote = RemoteData(remote_path=remote_path)
        remote.computer = computer

        if entry_point_name is not None:
            creator = CalcJobNode(computer=computer, process_type=entry_point)
            creator.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
            remote.base.links.add_incoming(creator, link_type=LinkType.CREATE, link_label='remote_folder')
            creator.store()

        return remote

    return _generate_remote_data


@pytest.fixture
def generate_workchain():
    """Generate an instance of a `WorkChain`."""

    def _generate_workchain(entry_point, inputs):
        """Generate an instance of a `WorkChain` with the given entry point and inputs.

        :param entry_point: entry point name of the work chain subclass.
        :param inputs: inputs to be passed to process construction.
        :return: a `WorkChain` instance.
        """
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager
        from aiida.plugins import WorkflowFactory

        process_class = WorkflowFactory(entry_point)
        runner = get_manager().get_runner()
        process = instantiate_process(runner, process_class, **inputs)

        return process

    return _generate_workchain


@pytest.fixture
def generate_workchain_cls():
    """Generate an instance of a `WorkChain` with abstract methods."""

    def _generate_workchain_cls(cls, inputs):
        """Generate an instance of a `WorkChain` with the given entry point and inputs.

        :param cls: work chain test class.
        :param inputs: inputs to be passed to process construction.
        :return: a `WorkChain` instance.
        """
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager

        runner = get_manager().get_runner()
        process = instantiate_process(runner, cls, **inputs)

        return process

    return _generate_workchain_cls


@pytest.fixture
def generate_workchain_force_sets(generate_workchain_cls, generate_force_sets_cls, generate_structure):
    """Generate an instance of a `ForceSetsWorkChain`."""

    def _generate_workchain_force_sets(append_inputs=None, structure_id='silicon', return_inputs=False):
        from aiida.orm import List

        structure = generate_structure(structure_id=structure_id)
        supercell_matrix = List([1, 1, 1])

        inputs = {'structure': structure, 'supercell_matrix': supercell_matrix}

        if return_inputs:
            return inputs

        if append_inputs is not None:
            inputs.update(append_inputs)

        process = generate_workchain_cls(generate_force_sets_cls(), inputs)

        return process

    return _generate_workchain_force_sets

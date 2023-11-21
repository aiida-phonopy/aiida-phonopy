---
myst:
    substitutions:
        pip: '[`pip`](https://pip.pypa.io/en/stable/index.html)'
        PyPI: '[PyPI](https://pypi.org/)'
---

# Get started

(installation-requirements)=

## Requirements

To work with `aiida-phonopy`, you should have:

- installed `aiida-core`.
- configured an AiiDA profile.

Please refer to the
[documentation](https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html)
of `aiida-core` for detailed instructions.

(installation-installation)=

## Installation

The Python package can be installed from the Python Package Index ({{ PyPI }}) or directly from the source:

::::{tab-set}

:::{tab-item} PyPI

The recommended method of installation is to use the Python package manager {{ pip }}:

```console
$ pip install aiida-phonopy
```

This will install the latest stable version that was released to the {{ PyPI }}.
:::

:::{tab-item} Source

To install the package from source, first clone the repository and then install using {{ pip }}:

```console
$ git clone https://github.com/aiida-phonopy/aiida-phonopy
$ pip install -e aiida-phonopy
```

The ``-e`` flag will install the package in editable mode, meaning that changes to the source code will be automatically picked up.
:::

::::

(installation-setup)=

## Setup

(installation-setup-computer)=

### Computer

To run Phonopy calculations on a compute resource, the computer should first be set up in AiiDA.
This can be done from the command line interface (CLI) or the Python application programming interface (API).
In this example, we will set up the `localhost`, the computer where AiiDA itself is running:

::::{tab-set}

:::{tab-item} CLI

To set up a computer, use the ``verdi`` CLI of ``aiida-core``.

```console
$ verdi computer setup -n -L localhost -H localhost -T core.local -S core.direct -w $AIIDA_PATH/workdir
```

After creating the localhost computer, configure the `core.local` transport using:

```console
$ verdi computer configure core.local localhost -n --safe-interval 0
```

Verify that the computer was properly setup by running:

```console
$ verdi computer test localhost
```
:::

:::{tab-item} API

To setup a computer using the Python API, run the following code in a Python script with `verdi run` or in the `verdi` shell:

```python
from aiida.orm import Computer
from pathlib import Path

computer = Computer(
    label='localhost',
    hostname='localhost',
    transport_type='core.local',
    scheduler_type='core.direct',
    workdir=Path('$AIIDA_PATH/workdir').absolute().as_posix()
).store()
computer.configure()
```
:::

::::

For more detailed information, please refer to the documentation [on setting up compute resources](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-set-up-a-computer).

:::{note}
The computer `workdir` is where AiiDA will run the calculations for codes set up for this computer.
The commands above will set the `workdir` of the `localhost` computer in `$AIIDA_PATH/workdir`.
:::

(installation-setup-code)=

### Code

To run a Phonopy code, it should first be setup in AiiDA.
This can be done from the command line interface (CLI) or the Python application programming interface (API).
In this example, we will setup the `phonopy` code that is installed on the computer where AiiDA is running:

::::{tab-set}

:::{tab-item} CLI

To setup a particular Phonopy code, use the ``verdi`` CLI of ``aiida-core``.

```console
$ verdi code create core.code.installed -n --computer localhost --label phonopy --default-calc-job-plugin phonopy.phonopy --filepath-executable $(which phonopy)
```
:::

:::{tab-item} API

To set up the Phonopy code using the Python API, run the following code in a Python script with `verdi run` or in the `verdi` shell:

```python
from aiida.orm import InstalledCode
import shutil

computer = load_computer('localhost')
code = InstalledCode(
    label='phonopy',
    computer=computer,
    filepath_executable=shutil.which('phonopy'),
    default_calc_job_plugin='phonopy.phonopy',
).store()
```
:::

::::

:::{important}
The CLI command and Python script above will automatically find the first `phonopy` binary in your `PATH` directories using the `which` command:

```console
which phonopy
```

If you have installed `phonopy` in a different python environment, make sure to pass the correct filepath to the executable and activate the environment correctly.
The latter can be done by adding the required commands to the `preprend_text` of the code.
:::

For more detailed information, please refer to the documentation [on setting up codes](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-setup-a-code).

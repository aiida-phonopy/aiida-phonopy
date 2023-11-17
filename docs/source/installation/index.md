---
myst:
    substitutions:
        pip: '`pip`'
---

# Get started

(installation-requirements)=

## Requirements

To work with `aiida-phonopy`, you should have:

- installed `aiida-core`
- configured an AiiDA profile.

Please refer to the
[documentation](https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html)
of `aiida-core` for detailed instructions.

(installation-installation)=

## Installation

The Python package can be installed from the Python Package index [PyPI](https://pypi.org/) or directly from the source:

::::{tab-set}

:::{tab-item} PyPI

The recommended method of installation is to use the Python package manager {{ pip }}:

```console
$ pip install aiida-phonopy
```

This will install the latest stable version that was released to PyPI.
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
$ verdi computer setup -n -L localhost -H localhost -T core.local -S core.direct
```

After creating the localhost computer, configure it using:

```console
$ verdi computer configure core.local localhost -n --safe-interval 0
```

Verify that the computer was properly setup by running:

```console
$ verdi computer test localhost
```
:::

:::{tab-item} API

To setup a computer using the Python API, run the following code in a Python script or interactive shell:

```python
from aiida.orm import Computer

computer = Computer(
    label='localhost',
    hostname='localhost',
    transport_type='core.local',
    scheduler_type='core.direct'
).store()
computer.configure()
```
:::

::::

For more detailed information, please refer to the documentation [on setting up compute resources](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-set-up-a-computer).

(installation-setup-code)=

### Code

To run a Phonopy code, it should first be setup in AiiDA.
This can be done from the command line interface (CLI) or the Python application programming interface (API).
In this example, we will setup the `phonopy` code that is installed on the computer where AiiDA is running:

::::{tab-set}

:::{tab-item} CLI

To setup a particular Phonopy code, use the ``verdi`` CLI of ``aiida-core``.

```console
$ verdi code create core.code.installed -n --computer localhost --label phonopy --default-calc-job-plugin phonopy.phonopy --filepath-executable /path/to/phonopy
```
:::

:::{tab-item} API

To setup particular Phonopy code using the Python API, run the following code in a Python script or interactive shell:

```python

from aiida.orm import InstalledCode

computer = load_computer('localhost')
code = InstalledCode(
    label='phonopy',
    computer=computer,
    filepath_executable='/path/to/phonopy',
    default_calc_job_plugin='phonopy.phonopy',
).store()
```
:::

::::

:::{important}
Make sure to replace `/path/to/phonopy` with the actual absolute path to the `phonopy` binary.

If you have it in a different python environment, make also sure to source the environment correctly.
This can be done specifying the *preprend_text* for the code.
:::

:::{hint}
You can find it by simply typing:

```console
$ which phonopy
```
:::

For more detailed information, please refer to the documentation [on setting up codes](https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html#how-to-setup-a-code).

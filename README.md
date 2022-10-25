# `aiida-phonopy`
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![PyPI version](https://badge.fury.io/py/aiida-phonopy.svg)](https://badge.fury.io/py/aiida-phonopy)
[![Build Status](https://github.com/aiida-phonopy/aiida-phonopy/workflows/aiida-phonopy/badge.svg?branch=develop&event=push)](https://github.com/aiida-phonopy/aiida-phonopy/actions)
[![Docs status](https://readthedocs.org/projects/aiida-phonopy/badge)](http://aiida-phonopy.readthedocs.io/)

This is the official AiiDA plugin for [Phonopy](https://phonopy.github.io/phonopy/index.html).

## Compatibility

From `v0.7.0` this plugin does not support retro-compatibility with previous versions,
due to a restructure of the package.

| Plugin | AiiDA | Phonopy |
|-|-|-|
| >=`v1.0.0` < v`2.0.0` | >=`v2.0.0` <`v3.0.0` |  >=`v2.14.0` <`v3.0.0` |
| >=`v0.7.0` < v`1.0.0` | >=`v1.6.0` <`v2.0.0` |  >=`v2.14.0` <`v3.0.0` |

## Installation

To install from PyPI, simply execute:

    pip install aiida-phonopy

or when installing from source:

    git clone https://github.com/aiida-phonopy/aiida-phonopy
    pip install .

## License

The `aiida-phonopy` plugin package is released under the MIT license.
See the `LICENSE.txt` file for more details.


## Acknowlegements

We acknowledge support from:
* the [U Bremen Excellence Chairs](https://www.uni-bremen.de/u-bremen-excellence-chairs) program funded within the scope of the [Excellence Strategy of Germany’s federal and state governments](https://www.dfg.de/en/research_funding/excellence_strategy/index.html);
* the [MAPEX](https://www.uni-bremen.de/en/mapex) Center for Materials and Processes.

<img src="https://raw.githubusercontent.com/aiida-phonopy/aiida-phonopy/develop/docs/source/images/UBREMEN.png" width="300px" height="54px"/>
<img src="https://raw.githubusercontent.com/aiida-phonopy/aiida-phonopy/develop/docs/source/images/MAPEX.jpg" width="300px" height="99px"/>

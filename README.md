# `aiida-phonopy`

This is the official AiiDA plugin for [Phonopy](https://phonopy.github.io/phonopy/index.html).

|    | |
|-----|----------------------------------------------------------------------------|
| Reference | [![DOI](https://img.shields.io/badge/DOI-10.1038/s41524024012363-purple.svg)](https://doi.org/10.1038/s41524-024-01236-3) |
|Latest release| [![PyPI version](https://badge.fury.io/py/aiida-phonopy.svg)](https://badge.fury.io/py/aiida-phonopy)[![PyPI pyversions](https://img.shields.io/pypi/pyversions/aiida-phonopy.svg)](https://pypi.python.org/pypi/aiida-phonopy) |
|Getting help| [![Docs status](https://readthedocs.org/projects/aiida-phonopy/badge)](http://aiida-phonopy.readthedocs.io/) [![Discourse status](https://img.shields.io/discourse/status?server=https%3A%2F%2Faiida.discourse.group%2F)](https://aiida.discourse.group/)
|Build status| [![Build Status](https://github.com/aiida-phonopy/aiida-phonopy/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/aiida-phonopy/aiida-phonopy/actions) [![Coverage Status](https://codecov.io/gh/aiida-phonopy/aiida-phonopy/branch/main/graph/badge.svg)](https://codecov.io/gh/aiida-phonopy/aiida-phonopy) |
|Activity| [![PyPI-downloads](https://img.shields.io/pypi/dm/aiida-phonopy.svg?style=flat)](https://pypistats.org/packages/aiida-phonopy) [![Commit Activity](https://img.shields.io/github/commit-activity/m/aiida-phonopy/aiida-phonopy.svg)](https://github.com/aiida-phonopy/aiida-phonopy/pulse)
|Community|  [![Discourse](https://img.shields.io/discourse/topics?server=https%3A%2F%2Faiida.discourse.group%2F&logo=discourse)](https://aiida.discourse.group/)

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

## How to cite

If you use this plugin for your research, please cite the following works:

* L. Bastonero and N. Marzari, [*Automated all-functionals infrared and Raman spectra*](https://doi.org/10.1038/s41524-024-01236-3), npj Computational Materials **10**, 55 (2024)

* A. Togo and I. Tanaka, [*First principles phonon calculations in materials science*](http://dx.doi.org/10.1016/j.scriptamat.2015.07.021), Scripta Materialia **108**, 1-5 (2015)

* S. P. Huber _et al._, [*AiiDA 1.0, a scalable computational infrastructure for automated reproducible workflows and data provenance*](https://doi.org/10.1038/s41597-020-00638-4), Scientific Data **7**, 300 (2020)

* M. Uhrin _et al._, [*Workflows in AiiDA: Engineering a high-throughput, event-based engine for robust and modular computational workflows*](https://www.sciencedirect.com/science/article/pii/S0927025620305772), Computational Materials Science **187**, 110086 (2021)

Please, also cite the underlying Quantum ESPRESSO and Phonopy codes references.

## License

The `aiida-phonopy` plugin package is released under the MIT license.
See the `LICENSE.txt` file for more details.


## Acknowlegements

We acknowledge support from:
* the [U Bremen Excellence Chairs](https://www.uni-bremen.de/u-bremen-excellence-chairs) program funded within the scope of the [Excellence Strategy of Germany’s federal and state governments](https://www.dfg.de/en/research_funding/excellence_strategy/index.html);
* the [MAPEX](https://www.uni-bremen.de/en/mapex) Center for Materials and Processes.

<img src="https://raw.githubusercontent.com/aiida-phonopy/aiida-phonopy/main/docs/source/images/UBREMEN.png" width="300px" height="108px"/>
<img src="https://raw.githubusercontent.com/aiida-phonopy/aiida-phonopy/main/docs/source/images/MAPEX.jpg" width="300px" height="99px"/>

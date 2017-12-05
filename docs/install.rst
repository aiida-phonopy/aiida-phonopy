
Installation
============

Use setup.py to install the plugin::

   python setup.py install --user

or PIP::

   pip install -e aiida-phonopy

Requirements
------------

* AiiDA v0.9
* phonopy v12.0
* seekpath

Additional optional requirements
--------------------------------

* aiida-vasp (to use VASP as a calculator)
* aiida-quantumespresso (to use QuantumESPRESSO as a calculator)
* aiida-lammps (to use LAMMPS as a calculator)

Setup workchains
----------------
The following files have to be copied (or linked) to aiida/workflows directory:

* wc_optimize.py
* wc_gruneisen.py
* wc_phonon.py
* generate_inputs.py
* parse_interface.py


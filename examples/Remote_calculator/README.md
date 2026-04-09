# Remote calculator example (branch fix)

This example demonstrates the branch fix for remote calculator serialization:
use a **calculator factory function** so the ASE `LAMMPSRun` calculator is
constructed on the remote worker, instead of trying to pickle a local temporary
directory path.

## What this fixes

When a pre-built `LAMMPSRun` calculator is serialized and sent to a remote
PythonJob, internal temporary paths can be captured and become invalid on the
remote machine. The fix is to pass a callable:

- `make_lammps_calculator()`

The callable is executed remotely and creates a fresh calculator instance with
the uploaded potential file available in the remote working directory.

## Files in this folder

- `remote_calculator.py`: runnable example using `phonopy.ase`
- `model.xyz`: input structure (`extxyz`)
- `nep89_20250409.txt`: potential file uploaded through PythonJob inputs

## Run

From this folder:

```bash
python remote_calculator.py
```

The script submits/runs `PhonopyAseWorkChain`, displays the band structure, and
writes `band.png`.

## Notes

- Update `pythonjob_inputs["computer"]`, queue name, and scheduler options for
  your environment.
- Ensure `phonopy@localhost` exists in your AiiDA profile, or replace it with
  your configured code label.

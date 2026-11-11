# suncat

SUNA + CTD toolkit: readers, quality control, calibration and gridded sections for Sea-Bird
SBE 19plus V2 CTD and SUNA V2 nitrate data. Extracted from the `ctd_pipeline` processing work.

## Two layers

- **Portable core** (`suncat.io`, `calc`, `qc`, `nitrate`, `productivity`, `sections`, `viz`):
  works on already-produced `.cnv` and L2/L3 products. Installs anywhere.
- **Optional SBE backend** (`suncat.sbe`): drives the proprietary, Windows-only SBEDataProcessing
  executables through the `.psa` files. Used on the acquisition machine only, never a hard
  dependency of the core.

## Install (development)

    pip install -e ".[dev]"
    pytest -q

Optional extras: `.[viz]` (figures), `.[netcdf]` (NetCDF output).

## Configuration

No paths are hardcoded. Point the package at your pipeline root with the `SUNCAT_CTD_ROOT`
environment variable, or pass `root=` to `suncat.config.for_cruise(...)`.

## Regression baseline (do this before any refactor)

    python tools/make_baseline.py     # freezes current outputs for chosen casts, locally

`tests/baseline/local/` is git-ignored because cruise products are confidential and this repo is
public. The baseline test skips when that folder is empty (including in CI). A small,
permission-cleared example cast will later live in `tests/baseline/example/` to run in CI too.

## Data policy

Raw cruise data, confidential metadata, credentials, and vendor material stay out of the repo.
The existing `.gitignore` and nbstripout discipline carry over unchanged.

## License

TODO: choose before first public release. MIT or BSD-3-Clause are the usual choices for
research Python.

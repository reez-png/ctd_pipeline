# Changelog

All notable changes to this project are documented here. This project follows
[Semantic Versioning](https://semver.org).

## [0.1.0] - 2026-10-27

Initial public release.

### Added
- Cruise-centric L0/L1/L2/L3 processing of Sea-Bird CTD casts (step01): the SBE chain DatCnv, Filter,
  Align, Cell Thermal Mass, Section (Python soak and equilibration-dip removal), Loop Edit, Derive, and
  Bin Average (1 m downcast and 1 s full cast), with bottle (auto-fire) mode.
- SUNA nitrate path: spectral QA/QC front-end (SUNA step00 inventory and step01 spectral QC producing
  the per-frame `suna_L1_clean` flag), CTD/SUNA time co-location and merge (step02), and the
  Zheng (2024) bias-corrected nitrate QC (step03) producing the L3 product.
- Productivity calibration (step07): in-situ chlorophyll fluorescence against extracted chlorophyll-a
  with non-photochemical-quenching handling, and a coded but dormant SBE 43 versus Winkler oxygen path;
  NPQ-sensitivity variants (step07_noremoval, step07_npq_comparison).
- Visualization and diagnostics: profile plots (step04), profile diagnostics (step05), and cross-shelf
  hydrographic sections (step06), with interactive explorers.
- Shared path and calibration resolver (`project_paths.py`); cruise `_TEMPLATE` skeleton and input
  templates for nitrate, chlorophyll, Winkler oxygen, and station coordinates.
- Synthetic `_EXAMPLE` cruise for end-to-end runs without real data, and a validation record
  (`docs/VALIDATION.md`).
- pytest suite (`ctd_lib` parsers and seawater functions, bottle table) and GitHub Actions CI across
  Python 3.10 to 3.12.

### Notes
- Nitrate and chlorophyll calibrations are provisional: small bottle counts (7 to 10) from a single
  cruise (P45_06). See `docs/VALIDATION.md` for the numbers and limitations.
- Cruise data is not included; it is available on reasonable request (see the README).

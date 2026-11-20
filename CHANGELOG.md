# Changelog

All notable changes to this project are documented here. Releases up to 0.1.0 used semantic versioning;
from November 2026 versions are calendar versions, `YYYY.M.D` (README, section 9).

## [Unreleased] - Nov 2026

### Changed
- Data levels: calibrate and flag at native resolution (L2), grid once (L3); every folder defined in
  `project_paths.py` and pinned by `tests/test_layout.py`; cruises moved with `tools/migrate_layout.py`.
- SUNA Stage 3 calibration stated per cruise in `metadata/suna_calibration.json` (`fit` or `adopted`), clock
  offsets in `suna_time_offsets.csv`, slow casts in `suna_response.csv`, each written by its tool.
- Homogeneous-layer SUNA calibration needs 3 clean frames on each leg (P45_06: gain 0.877 to 0.862).
- step07 applies the chlorophyll calibration to every scan; its printed inference no longer reads the gain as a
  factory error (the ECO V2 sent Legacy counts the `.xmlcon` coefficients do not match).
- step04's readable tables moved to `L3/CTD/profiles_readable`; step04 plots the production nitrate and the
  calibrated chlorophyll from step08's grid (it had plotted the cross-check nitrate), and its tables carry both.
- step02 finds SUNA files whatever the case of `.csv` (none were found on a case-sensitive file system).
- step08 counts the daylight-rule NPQ bins on casts without a PAR channel too.
- `cruises/_EXAMPLE` runs the whole path after step01 (Derived files, bin-type headers, `suna_calibration.json`,
  `stations.csv`), checked by `tests/test_example_cruise.py`; `station_coords.csv` retired for `stations.csv`.

### Added
- step02b (SUNA response time and the lag-aware bottle test), SUNA step04 (clock lag), step08 (L3 gridding).
- Per-bin down/up nitrate check in step08 (`nitrate_updown_uM`, `nitrate_hyst_flag`).
- SUNA spectral QC: surface breaks inside the CTD's wet interval, out-of-water windows per cast.
- Bottle pairing with hand flags and a temperature check; stop detection in the CTD record (`suncat.stops`).
- `chl_above_cal` and `chl_cal_max`; NPQ by PAR beside the daylight rule (`chl_npq_risk_par`), and step07's
  comparison fit under it.
- Tools: `make_stations_table.py`, `inspect_eco_channels.py`, `make_synthetic_cruise.py`, `run_notebook.py`,
  the `set_suna_*` scripts, response and lag diagnostics.
- Docs: the instrument suite and its configuration, adapting the pipeline to other instruments, the two cruises'
  results and open questions (`docs/VALIDATION.md`).

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

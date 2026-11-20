# SUNCAT: SUNA and CTD Analysis Toolkit

![tests](https://github.com/reez-png/ctd_pipeline/actions/workflows/tests.yml/badge.svg)

**SUNCAT** is the cruise-centric processing of Sea-Bird CTD casts and a rosette-mounted SUNA nitrate sensor for
the Ocean Margins Initiative (OMI) on the Ghana shelf. It follows the June 2026 OMI specification (Marrec, Lucas),
the SUNA quality control of Zheng et al. (2024), and a data-level scheme agreed in November 2026: calibrate and
flag at native resolution (L2), then grid once (L3).

The repository holds two things:

- **`suncat`**, an installable Python package (`src/suncat/`) with the tested core: readers, sectioning, seawater
  functions, SUNA spectral QC, bottle pairing, response and calibration fits, gridding, and an optional
  Windows-only backend that drives SBE Data Processing.
- **The cruise notebooks and scripts** (`processing_scripts/notebook/`, `tools/`) that use it to turn raw
  instrument files into calibrated, flagged products, organised by cruise under L0 to L3.

Built on, and checked against, two cruises: P45_06 (April 2026, 11 casts) and P45_07 (31 casts). What the pipeline
recovers on them, and what it cannot yet settle, is in [`docs/VALIDATION.md`](docs/VALIDATION.md).

Contents: [instruments](#1-the-instruments-it-is-built-for) ·
[install](#2-install) · [data levels](#3-data-levels-and-layout) · [running a cruise](#4-running-a-cruise) ·
[the steps](#5-the-steps) · [results and findings](#6-results-and-findings) ·
[other instruments](#7-adapting-it-to-another-ctd-or-sensor-suite) · [tools](#8-helper-tools) ·
[tests](#9-tests-baseline-and-ci) · [known gaps](#10-known-gaps-and-planned-work) · [references](#references)

---

## 1. The instruments it is built for

One instrument suite, the same on P45_06 and P45_07, under one configuration file (`19-8460.xmlcon`):

| Instrument | Serial | Calibration | Reaches the pipeline as |
|---|---|---|---|
| Sea-Bird SBE 19plus V2 SeaCAT, profiling mode, 4 Hz, strain-gauge pressure | 8460 | temperature and conductivity 3 Aug 2025; pressure 4 Sep 2025 | `tv290C`, `c0S/m`, `prdM` |
| Sea-Bird SBE 43 dissolved oxygen | 4704 | 8 Jul 2025 | `sbeox0V`, `sbeox0Mm/L`, `sbox0Mm/Kg` |
| Satlantic log PAR | 2419 | 28 Jul 2025 | `par/sat/log` |
| WET Labs C-Star transmissometer, 25 cm path | 2655 | 10 Sep 2025 | `CStarAt0` |
| Sea-Bird ECO V2 (chlorophyll at 470 and 435 nm excitation, backscatter at 700 nm), on the CTD's RS-232 port | 10043 | characterization 8 Sep 2025 | `flECO-AFL` (470 nm chlorophyll); `flECO-AFL1` (most likely 435 nm, not used); backscatter not exported |
| Sea-Bird SUNA V2 nitrate, open path, own logger and clock | 2503 | `SNA2503C.cal`; saltwater verification 20 Aug 2025 | raw SATSLF frames (`.CSV`), merged by time |
| Sea-Bird SBE 55 ECO water sampler (Niskins) | not in the `.xmlcon` | | `.afm` auto-fire log per cast (P45_07 only) |

Two facts about this suite shape the processing (details in
[`docs/PROCESSING_PARAMETERS.md`](docs/PROCESSING_PARAMETERS.md)):

- **The ECO V2 sent Legacy (ECO Classic) counts**, which saturate at 4130, while the `.xmlcon` holds the
  characterization sheet's coefficients rather than the Legacy ones the sensor's `GetCC` reports. The bottle
  calibration in step07 absorbs the difference; the `.cnv` fluorescence itself cannot be read as the factory scale.
- **The SUNA logs on its own clock** and responds with a lag (a median response time of 7 s on P45_06, much
  slower on a few casts).
  Its frames are matched to the CTD by UTC time, its clock offset is measured, and its response is assessed per
  cast before any bottle calibration.

---

## 2. Install

```bash
git clone https://github.com/reez-png/ctd_pipeline.git
cd ctd_pipeline
python -m venv .venv
.venv\Scripts\activate           # Windows (macOS/Linux: source .venv/bin/activate)
pip install -r requirements.txt
pip install -e ".[dev]"          # the suncat package (editable) and pytest
pytest -q
```

Developed on Python 3.14; CI runs the tests on 3.10 to 3.12. `requirements.txt` gives minimum versions,
`requirements-lock.txt` the exact tested environment. The notebook kernel must be the environment `suncat` is
installed in. step01 also needs Sea-Bird **SBE Data Processing** (Windows); everything after it is pure Python.

**Try it without cruise data.** `cruises/_EXAMPLE` (three synthetic casts) runs the whole path after step01, from
the SUNA QA to step08 and the figures, and recovers what was planted in it; the order and the expected numbers are
in [`cruises/_EXAMPLE/README_EXAMPLE.md`](cruises/_EXAMPLE/README_EXAMPLE.md), and `tests/test_example_cruise.py`
checks them. `tools/make_synthetic_cruise.py <root>` writes cruises built to differ from P45 (another hemisphere,
other column names, cast numbering). The calibration files must be reachable: either run from the repo, or copy
`metadata/calibration/` next to `cruises/` under your data root.

```bash
python tools/run_notebook.py processing_scripts/notebook/suna/step01_suna_spectral_qc.ipynb --set CRUISE_ID=_EXAMPLE
```

To keep data outside the repository, point the notebooks and scripts at a data root holding `cruises/` with the
environment variable `CTD_PIPELINE_ROOT` (`suncat.config` and `tools/make_baseline.py` read `SUNCAT_CTD_ROOT`).

**Data availability.** This repository holds the processing code only. The cruise measurements (raw `.hex`, SUNA
frames, L0 to L3 products, discrete nutrient, chlorophyll and oxygen tables, station positions) are excluded by
`.gitignore`; they are part of the OMI P45 cruise dataset and available on reasonable request. The Sea-Bird PSA
files are not tracked either; [`docs/PROCESSING_PARAMETERS.md`](docs/PROCESSING_PARAMETERS.md) is the tracked
record of every processing parameter.

---

## 3. Data levels and layout

| Level | Meaning | Resolution |
|---|---|---|
| **L0** | raw, as logged: CTD `.hex` (and `.afm`), SUNA frames | instrument |
| **L1** | cleaned per cast: SBE stages DatCnv to Derive, the CTD 1 s time base, SUNA spectral-QC tables | native |
| **L2** | merged, calibrated and flagged: CTD and SUNA on one time base, nitrate on the bottle scale, chlorophyll calibrated, flags carried, down and up legs identified | native |
| **L3** | gridded: one table per cast on the SBE 1 dbar down-cast bins, plus the SBE bin files and readable tables | 1 dbar |

Every folder is defined once, in `project_paths.py` (`PRODUCTS`, `product_dir()`), and pinned by
`tests/test_layout.py`. The full map, and how a cruise was moved from the old layout, is in
[`docs/DATA_LEVELS.md`](docs/DATA_LEVELS.md). Never build `L1`/`L2`/`L3` paths by hand in a step.

```text
cruises/<CRUISE_ID>/
├── L0/CTD/                    *.hex (+ *.afm beside each, same name)
├── L0/SUNA/                   raw SUNA frames
├── L1/CTD/                    DatCnv *.cnv; Filter/ Align_CTD/ Cell_tm/ Section/ Loop_Edit/ Derived/ CTD_1s/
├── L1/SUNA/                   per-frame spectral QC (suna_L1_clean)
├── L2/SUNA/                   <cast>_SUNA_1s.csv (merged), <CRUISE>__suna_L2_colocated.csv
├── L2/SUNA_calibrated/        step03_suna_bottle_calibration's corrected profiles (a cross-check, not plotted)
├── L2/nitrate/                <cast>_nitrate_QC.csv: the production nitrate (step03_nitrate_qc)
├── L2/productivity/           calibrated chlorophyll, every scan (step07); productivity_noNPQ/ beside it
├── L2/CTD/Bottle/             .btl (bottle mode)
├── L3/CTD/CTD_1m_down/        SBE Bin Average, 1 dbar down-cast
├── L3/CTD/profiles_readable/  <cast>_readable.csv (step04); combined_all_casts_*.csv in L3/CTD/
├── L3/gridded/                <cast>_L3_1m_down.csv: CTD, nitrate and chlorophyll with their flags (step08)
├── metadata/                  per-cruise settings and discrete samples (below)
├── _audit/                    logs, QC tables, fit tables
└── _figures/                  qc/ and profiles/ (figures are not data)
```

All of `L0` to `L3`, `_audit` and `_figures` are git-ignored. Of `metadata/`, only `suna_cast_map.csv` and the
`*_TEMPLATE.csv` files are tracked.

---

## 4. Running a cruise

Copy `cruises/_TEMPLATE/` to `cruises/<CRUISE_ID>/`, put the raw files in `L0/CTD` and `L0/SUNA` (see
[`docs/CRUISE_DATA_COLLECTION.md`](docs/CRUISE_DATA_COLLECTION.md)), and run the steps in this order. Notebooks can
be run from Jupyter or the command line:

```bash
python tools/run_notebook.py processing_scripts/notebook/step08_grid_l3.ipynb --set CRUISE_ID=P45_07
```

| # | Step | Writes | Then |
|---|---|---|---|
| 1 | `step01_run_sbe_processing_v2.py` (PyCharm; dry run first) | L1/CTD stages, L3/CTD/CTD_1m_down, L2/CTD/Bottle | |
| 2 | `suna/step00_suna_inventory`, `suna/step01_suna_spectral_qc` | L1/SUNA flags | |
| 3 | `step02_merge_ctd_suna`, `suna/step02_suna_colocation` | L2/SUNA | |
| 4 | `suna/step04_suna_time_lag` | a proposed clock offset | if the SUNA clock is off: `tools/set_suna_time_offset.py`, rerun row 3 |
| 5 | `suna/step02b_suna_response_correction` | response time per cast, bottle test | `tools/set_suna_slow_casts.py` |
| 6 | `build_bottle_nitrate_table.py` (with `.btl`) or `build_bottle_table.py` (without) | `metadata/bottle_nitrate/` | |
| 7 | `suna/step03_suna_bottle_calibration` | homogeneous-layer fit, L2/SUNA_calibrated | `tools/set_suna_calibration.py` (`fit` or `adopted`) |
| 8 | `step03_nitrate_qc` | L2/nitrate | |
| 9 | `productivity/step07_productivity_calibration` (and `_noremoval`, `_npq_comparison`) | L2/productivity | |
| 10 | `tools/make_stations_table.py <CRUISE> --write` | `metadata/stations.csv` | |
| 11 | `step08_grid_l3` | L3/gridded | |
| 12 | `step04_plot_profiles`, `step05_profile_diagnostics`, `step06_hydrographic_sections`; `combine_readable_profiles.py` | figures, readable tables | |

Per-cruise settings live in `cruises/<CRUISE_ID>/metadata/`, written by the tools rather than edited by hand, so a
run command or a notebook default cannot change them silently:

| File | Says | Written by | Read by |
|---|---|---|---|
| `suna_cast_map.csv` | which SUNA file belongs to which cast | you (or step02's proposal) | step02, step03 |
| `suna_time_offsets.csv` | the SUNA clock offset, per cast | `tools/set_suna_time_offset.py` | step02 |
| `suna_response.csv` | casts with a slow SUNA response | `tools/set_suna_slow_casts.py` | step02b, step06, step08 |
| `suna_calibration.json` | the Stage 3 calibration: `fit` from this cruise's bottles, or `adopted` with its source | `tools/set_suna_calibration.py` | step03 (required) |
| `stations.csv` (or `<CRUISE>_CTD_cast_summary.csv`) | each cast's start time and position | `tools/make_stations_table.py` | step06, step07, step08 (solar time) |
| `bottle_nitrate/`, `chla/`, `winkler_oxygen/` | discrete samples (templates show the columns) | you, the bottle builders | step02b, step03, step07 |

**Before the first run**

- [ ] The PSAs in `psa/` (local, rebuilt in the Sea-Bird GUI): `01_datcnv`, `05_filter`, `02_alignctd`, `04_celltm`,
      `06a_loopedit_percentmean` (and `06b`), `07_derive`, `08b_binavg_1m_down`, `08d_binavg_1s_full`; each with
      blank input and output directories and name append (`fix_psa_output_paths.py`). The bottle PSA is generated.
- [ ] `19-8460.xmlcon` and `SNA2503C.cal` in `metadata/calibration/` (or the cruise's `Instrument_Calibration_Files/`).
- [ ] For bottle mode, each `.afm` beside its `.hex` with the same name, and `APPLY_BOTTLES` set for the cruise.
- [ ] `suna_cast_map.csv` filled, and the bottle and chlorophyll tables in the template's columns.
- [ ] step01 dry run first; the same `CRUISE_ID` in every step.

---

## 5. The steps

### step01: Sea-Bird processing (L0 to L1, and the 1 dbar bins at L3)

```text
DatCnv -> Filter -> Align CTD -> Cell Thermal Mass -> Section (Python) -> Loop Edit -> Derive -> Bin Average
```

- **Filter** (before Align): pressure 1.0 s; temperature and conductivity 0.5 s.
- **Align CTD**: temperature +0.5 s, conductivity 0 s, SBE 43 oxygen +5 s (provisional until Winkler data).
- **Cell Thermal Mass**: alpha 0.03, 1/beta 7.0 (Sea-Bird defaults; an open decision).
- **Section**, a Python step: removes the surface soak and the equilibration dip before Loop Edit, whose
  running-maximum-pressure rule the dip would otherwise poison (`suncat.section.descent_start`).
- **Loop Edit**: percent of mean speed (3 s window, 20%), minimum velocity 0.10 m/s.
- **Bin Average**: 1 dbar pressure bins, down-cast only (science), and 1 s bins over the full cast (SUNA merge).

Run it with `RUN_SBE_COMMANDS = False` first, check `_audit/01_psa_diagnostic.csv` (every PSA's `warnings` empty),
then set it `True`. Each module has a timeout, and the work folder copes with read-only and locked files.
**Bottle mode** (`APPLY_BOTTLES = True`, for casts with an `.afm`): Data Conversion writes `.cnv` and `.ros`,
Bottle Summary writes `.btl`; a cast without an `.afm` falls back to the profile PSA. P45_06 has no `.afm` (its
Niskins were fired on the move), so its bottle depths are the sample sheet's targets.

Lessons from the build, so they are not reintroduced: Section must come before Loop Edit (the order is not
interchangeable); the Section writer keeps the `.cnv` CRLF line endings byte for byte (LF-only makes SBE report
"Header line length exceeds buffer length"); "Create both data and bottle file" is `CreateFile=2` in the PSA
(`1` writes only the `.ros`); PSAs re-saved from the SBE GUI regain hardcoded output paths, so blank them again.

### The SUNA path

- **step00 / step01 (SUNA QA).** Frame inventory and instrument health; then per-frame spectral QC: warm-up,
  saturation, low signal, fit RMSE, humidity, an out-of-water gate built from each cast's own in-water frames, and
  surface breaks (the package back above 0.5 dbar for 30 s or more inside the CTD's wet interval). Writes the
  per-frame `suna_L1_clean` flag.
- **step02 (merge) and step02_suna_colocation.** SUNA frames onto the CTD 1 s time base by UTC, with the clock
  offset from `suna_time_offsets.csv`, refusing to interpolate across gaps; the co-located product gives every
  frame a pressure, a cast and a leg.
- **step04_suna_time_lag.** The SUNA-to-CTD lag by three methods, with confidence intervals, a drift test, and an
  independent check from water entry and exit.
- **step02b.** The SUNA's response time per cast from the disagreement between down- and up-cast (no bottles),
  then whether lag-corrected profiles match the bottles well enough to calibrate from them.
- **step03_suna_bottle_calibration.** The homogeneous-layer fit: only bottles where the down- and up-cast agree
  within 2 uM (at least 3 clean frames on each leg), with leave-one-out, so the lag cannot enter the gain.
- **step03_nitrate_qc** (production). Zheng et al. (2024): Stage 1 the onboard nitrate (the TSP re-derivation
  is work in progress), Stage 3 the calibration stated in `suna_calibration.json`, Stage 4 the assessment against
  bottles. Writes `L2/nitrate/<cast>_nitrate_QC.csv`, which step08 grids.

### step07: chlorophyll (and, when Winkler data exist, oxygen)

Fits `.cnv value = offset + gain x extracted chl-a` through bottles matched by depth (median of 1 m cells within
2 m on the down-cast), and applies it to every scan. Near-surface daytime bottles (10 m or shallower, 07 to 17
local solar time) are left out for quenching; the fit under a PAR rule (PAR at least 100 umol photons m-2 s-1 at
the bottle) is printed beside it. Each scan carries `chl_above_cal` where the value exceeds the highest bottle in
the fit. A fluorometer stuck at one value is treated as dead for that cast. `step07_noremoval` keeps the
quenched points; `step07_npq_comparison` overlays both. The oxygen path is coded and dormant.

### step08: L3 gridding

One table per cast on the SBE 1 dbar down-cast bins: the CTD's own bin averages unchanged, and the production
nitrate and calibrated chlorophyll averaged onto the same bins from flagged-clean samples. Per bin it adds:
- `nitrate_updown_uM` and `nitrate_hyst_flag`: down- and up-cast nitrate read at the bin centre from their own frames
  (within 2.5 dbar, 3 frames a leg), flagged where they differ by more than 2 uM;
- `chl_npq_risk` (daylight rule) and `chl_npq_risk_par` (PAR rule);
- `chl_above_cal`;
- the cast's SUNA response, with a per-cast summary in `_audit/gridding/`.

### step04 to step06: figures and tables

step04 profile figures and readable tables, with the nitrate and calibrated chlorophyll from step08's grid; step05 cline, DCM and mixed-layer diagnostics; step06 cross-shelf
sections, with step08's nitrate and red ticks where the down/up check flags a bin.

### The CTD QC layer (step00a to step00e)

Diagnostics beside the chain that never change production data: Filter effect (step00a), Filter and Loop Edit
flag maps (step00b), Align advance sweeps (step00c), a stage-by-stage change audit that reproduces step01 to 0.0
difference (step00d), and soak against no-soak near-surface quality (step00e, `soak_quality_analysis`).

---

## 6. Results and findings

| | P45_06 | P45_07 |
|---|---|---|
| SUNA Stage 3 (`N_QC = (N - a0) / a1`) | fit: a1 0.886, a0 -0.59 (27 bottles) | adopted: a1 0.89, a0 -0.39 |
| Chlorophyll (`.cnv value = offset + gain x chl-a`) | gain 1.123, offset -0.002, R2 0.67, n 15 | gain 0.515, offset +0.189, R2 0.95, n 18 |
| Gridded nitrate bins with down and up more than 2 uM apart | 30% | 41% |

What the two cruises showed, each with its evidence in [`docs/VALIDATION.md`](docs/VALIDATION.md):

- **The SUNA lags, and the down-cast is what is gridded.** Down- and up-casts differ by more than 2 uM in 30 to 41%
  of assessed bins, almost all in the nitracline; on P45_07 the up-cast reads up to 16 uM higher, probably water held in the vertically
  mounted SUNA's optical path. Bottles stay the nitracline of record, and P45_07's calibration is adopted from
  the homogeneous layer, supported by lag-aware fits to 52 bottles.
- **Three casts surfaced mid-cast** (151 to 202 s above 0.5 dbar), two of them among the three casts with held
  water: an observation, not yet a cause.
- **The SUNA's field gain (0.86 to 0.89) disagrees with its factory saltwater check (1.00).** Open.
- **P45_06's bottle depths are targets**; the homogeneous-layer check shows that does not move its calibration.
- **The ECO's chlorophyll gain differs 2.2-fold between cruises**, and P45_06's upper-layer bottles read low
  against its deeper ones; the `.xmlcon` lacks the ECO's Legacy coefficients. P45_07 lost the ECO entirely from
  cast 13 (both channels at 0 counts for 19 casts).
- **The daylight and PAR rules for quenching** agree on 87 to 95% of bins, but move the bottle fits little; the
  daylight rule is kept and both flags are carried.

The field practices these findings call for (clock synchronisation, holds at each bottle, the auto-fire log,
SUNA mounting and flushing, ECO coefficients, monitoring every channel in real time) will be collected in a CTD
operations guide, in preparation.

---

## 7. Adapting it to another CTD or sensor suite

The pipeline is general where it reads data and specific where it was configured. To use it with other
instruments:

1. **Another Sea-Bird CTD** (SBE 911plus, 25plus, another 19plus): step01 works unchanged; supply its `.xmlcon`
   and rebuild the PSAs for its sensors (Align advances, Cell Thermal Mass and Filter differ by instrument; the
   step00a to step00c notebooks help choose them). Record every value in `docs/PROCESSING_PARAMETERS.md`.
2. **Column names.** Readers resolve pressure from `prdM`/`prDM`/`prSM` and a few variants, and PAR from any column
   starting with `par`. Instrument-specific names are listed in each notebook's settings: step07's `TARGETS`
   (`sensor_cols`, for example `flECO-AFL`, `wetCDOM`, `flC`; oxygen `sbeox0Mm/L`), step04's `CTD_VARIABLES` and
   `LAB_RENAME`. Add your sensors' names there.
3. **The 1 dbar grid.** step08 requires SBE Bin Average output in pressure bins (`binavg_bintype = decibars`).
4. **A non-Sea-Bird CTD** cannot use step01. Its data can enter at L1 if converted to `.cnv`-style files with
   the same column names and a Loop Edit `flag` column, and binned to 1 dbar with the bin-type header; the readers
   in `suncat.io.cnv` define what is expected.
5. **Another nitrate sensor or SUNA.** Put its `.cal` in `metadata/calibration/` and name it in the SUNA
   notebooks' settings. The spectral QC thresholds (`QCConfig` in `suncat.nitrate.spectral_qc`) assume a SUNA V2's
   256-channel spectra; check them on a few casts with step00 and step01 before trusting the flags.
6. **No SUNA.** step01, step07, step08 (CTD and chlorophyll only) and step04 to step06 run without it.
7. **Fluorometers.** If the ECO reports Legacy counts, enter its `GetCC` coefficients in the `.xmlcon`; either
   way, step07's bottle calibration is what makes the chlorophyll quantitative.
8. **Location and time.** Nothing assumes the Gulf of Guinea: daylight is local solar time from each cast's
   longitude (`stations.csv`), and the synthetic cruises test the southern hemisphere, east longitudes and the
   equator.

---

## 8. Helper tools

| Script | What it does |
|---|---|
| `tools/run_notebook.py` | Runs a notebook from the command line with settings overridden (`--set NAME=value`), optionally at another git revision; never writes outputs into the repo |
| `tools/set_suna_time_offset.py`, `set_suna_slow_casts.py`, `set_suna_calibration.py` | Write (or `--show`) the per-cruise SUNA settings in section 4 |
| `tools/make_stations_table.py` | `stations.csv` from the `.cnv` NMEA header fixes, else the nutrient file through the pairing audit |
| `tools/estimate_response_time.py` | A sensor's response time per cast from its own profiles (SUNA or chlorophyll) |
| `tools/diagnose_suna_lag.py` | The SUNA's clock, kind of response and bottle match in one run |
| `tools/compare_match_estimators.py`, `compare_bottle_pairs.py` | Which sensor value best represents a bottle's water; step02b's pairs against the lag diagnosis' |
| `tools/inspect_eco_channels.py` | What the ECO's third channel carries, and scans at the Legacy full range |
| `tools/make_synthetic_cruise.py` | A synthetic cruise built to differ from P45 wherever the code could assume it |
| `tools/make_baseline.py` | Freezes the current products as the regression baseline (local only) |
| `tools/migrate_layout.py` | Moves a cruise between the data-level layouts, checked file by file |
| `build_bottle_nitrate_table.py` | Bottle mode: `.btl` firings to `<CRUISE>_bottle_ctd.csv`, joined to lab nitrate by cast and Niskin |
| `build_bottle_table.py` | Without `.btl`: pairs bottles to casts by time, carries hand flags, reports the CTD's stops against the sheet depths |
| `combine_readable_profiles.py` | Stacks step04's readable tables into one cruise table |
| `cnv_to_csv.py`, `verify_reproducibility.py`, `fix_psa_output_paths.py`, `step03_run_report.py` | `.cnv` to CSV; compare two step01 runs; blank PSA output paths; a readable step03 report |

---

## 9. Tests, baseline and CI

`pytest -q` runs the suite (245 tests on a machine with the cruise data, Nov 2026; tests that need cruise data or
the Sea-Bird software skip elsewhere). GitHub Actions runs it on Python 3.10 to 3.12 on every push. Notebook tests
run the real notebooks on synthetic cruises through `tools/run_notebook.py`. `tools/make_baseline.py` freezes the
current L2 products under `tests/baseline/local/` (git-ignored); `tests/test_baseline_products.py` then reports any
product that changes, so a change that moves a number is reviewed before the baseline is frozen again.

Workflow: each piece of work on a fresh branch from `main`, one squash-merged pull request each. Versions are
calendar versions, `YYYY.M.D` (non-zero-padded), tagged `vYYYY.M.D` on `main` and set in `pyproject.toml`.

---

## 10. Known gaps and planned work

- **Oxygen.** The SBE 43 path in step07 waits for Winkler titrations.
- **SUNA Stage 1 re-derivation** (`STAGE1_MODE = "tsp"`, Sakamoto et al. 2009) is work in progress and needs
  `SNA2503E.cal`, the file the onboard nitrate was computed with; keep `"onboard"`.
- **Open questions** carried in VALIDATION: the SUNA's factory against field gain, the ECO's Legacy coefficients
  and third channel, P45_06's upper-layer chlorophyll, the Cell Thermal Mass coefficients.
- **Package migration.** Parts of the notebooks still carry their own code that belongs in `suncat`.
- **Two root variables.** `project_paths.py` reads `CTD_PIPELINE_ROOT`, `suncat.config` reads `SUNCAT_CTD_ROOT`; one
  should go.

---

## References

**CTD and SUNA nitrate**

- Johnson, K. S., & Coletti, L. J. (2002). In situ ultraviolet spectrophotometry for high resolution and long-term
  monitoring of nitrate, bromide and bisulfide in the ocean. *Deep-Sea Research I, 49,* 1291 to 1305.
- Sakamoto, C. M., Johnson, K. S., & Coletti, L. J. (2009). Improved algorithm for the computation of nitrate
  concentrations in seawater using an in situ ultraviolet spectrophotometer. *Limnol. Oceanogr.: Methods, 7,*
  132 to 143.
- Sakamoto, C. M., et al. (2017). Pressure correction for the computation of nitrate concentrations in seawater
  using an in situ ultraviolet spectrophotometer. *Limnol. Oceanogr.: Methods, 15,* 897 to 902.
- Plant, J. N., et al. (2023). Updated temperature correction for computing seawater nitrate with in situ
  ultraviolet spectrophotometer and submersible ultraviolet nitrate analyzer nitrate sensors. *Limnol. Oceanogr.:
  Methods, 21,* 581 to 593.
- Zheng, B., et al. (2024). Bias-corrected high-resolution vertical nitrate profiles from the CTD rosette-mounted
  SUNA. *Limnol. Oceanogr.: Methods, 22,* 889 to 902.
- Sea-Bird Scientific. SBE Data Processing manual and application notes; ECO V2 user manual (Legacy output,
  appendix A); SUNA V2 manual.
- GO-SHIP / WOCE Hydrographic Programme. CTD operations guidance, including the surface soak.

**Chlorophyll fluorescence and oxygen**

- Roesler, C., et al. (2017). Recommendations for obtaining unbiased chlorophyll estimates from in situ
  chlorophyll fluorometers. *Limnol. Oceanogr.: Methods, 15,* 572 to 585.
- Lorenzen, C. J. (1966). A method for the continuous measurement of in vivo chlorophyll concentration.
  *Deep-Sea Research, 13,* 223 to 227.
- Holm-Hansen, O., et al. (1965). Fluorometric determination of chlorophyll. *ICES J. Mar. Sci., 30,* 3 to 15.
- Welschmeyer, N. A. (1994). Fluorometric analysis of chlorophyll a in the presence of chlorophyll b and
  pheopigments. *Limnol. Oceanogr., 39,* 1985 to 1992.
- Kiefer, D. A. (1973). Fluorescence properties of natural phytoplankton populations. *Marine Biology, 22,*
  263 to 269.
- Falkowski, P., & Kiefer, D. A. (1985). Chlorophyll a fluorescence in phytoplankton: relationship to
  photosynthesis and biomass. *J. Plankton Research, 7,* 715 to 731.
- Cullen, J. J. (1982). The deep chlorophyll maximum: comparing vertical profiles of chlorophyll a. *Can. J.
  Fish. Aquat. Sci., 39,* 791 to 803.
- Schmechtig, C., et al. (2018). BGC-Argo quality control manual for the chlorophyll-a concentration.
- Carpenter, J. H. (1965). The Chesapeake Bay Institute technique for the Winkler dissolved oxygen method.
  *Limnol. Oceanogr., 10,* 141 to 143.

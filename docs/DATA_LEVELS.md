# Data levels (L0 to L3): definitions, layout and migration

Agreed November 2026 following Prof. Amala's and Prof. Mahu's review: **calibrate and flag before binning**.
This document is the reference for what each level means and where every product lives. The layout itself
is defined once in code, in `project_paths.py` (`PRODUCTS`, `product_dir()`), and is pinned by
`tests/test_layout.py`; the table below is generated from that registry.

## Definitions

| Level | Meaning | Resolution | Examples |
|---|---|---|---|
| **L0** | Raw, as logged | instrument | CTD `.hex` (+ `.afm`), SUNA frames, bottle log sheets |
| **L1** | Cleaned per cast | native (CTD ~4 Hz, SUNA ~1 Hz; 1 s time base) | DatCnv `.cnv`, SBE stages Filter to Derive, CTD 1 s time base, SUNA spectral QC tables, out-of-water removed |
| **L2** | Merged, calibrated and flagged | native (1 s) | CTD and SUNA merged and time-aligned, nitrate on the bottle scale, chl-a calibrated against discrete samples, quality flags carried, down/up legs identified |
| **L3** | Gridded | 1 m down-cast bins | science-ready gridded CTD, nitrate and chl-a profiles |

**The principle:** averaging happens once, at L3, and only on calibrated, flagged data. Bad frames are
excluded *before* a bin mean is formed, and per-sample corrections (for example the temperature and salinity
dependence of the SUNA nitrate, which varies within a 1 m bin) are applied at native resolution. This also
matches the common convention in ocean data management that a gridded product carries the highest level.

## Decisions

1. **SBE processing stages** (Filter, Align CTD, Cell Thermal Mass, Section, Loop Edit, Derive) are L1:
   cleaned per cast, at native resolution. `L1/CTD/<stage>/`.
2. **The CTD 1 s product** is L1. It resamples the CTD to the SUNA's ~1 Hz rate as the time base for the
   merge; it is not a depth grid. `L1/CTD/CTD_1s/`.
3. **L3 gridding:** the CTD 1 m down-cast grid stays the SBE Bin Average output (community standard;
   flagged Loop Edit scans are excluded by Bin Average). Calibrated nitrate and chl-a are gridded in Python
   onto the same 1 m bins, excluding flagged samples, into one combined L3 table per cast.
4. **L3 legs:** down-cast only, as the science product. Up-cast data stay available at L2 with their leg flag.
5. **Figures** are not data levels: they move to `_figures/` (git-ignored, like every cruise folder).

## Level map (generated from `project_paths.PRODUCTS`)

Bold marks a folder that moves. `v1` is the layout the pipeline was built on; `v2` is the agreed layout.

| Product (`key`) | v1 folder (built on) | v2 folder (agreed) | Written by | Read by |
|---|---|---|---|---|
| `ctd_raw` | L0/CTD | L0/CTD | field data | step01 |
| `suna_raw` | L0/SUNA | L0/SUNA | field data | step01_suna_spectral_qc, step02 (merge, co-location), step03 |
| `ctd_datcnv` | L1/CTD | L1/CTD | step01 (DatCnv) | step00c, step00d, step00e, soak_quality, step02_suna_colocation, step04_suna_time_lag |
| `ctd_filter` | L2/CTD/Filter | **L1/CTD/Filter** | step01 | step00d |
| `ctd_align` | L2/CTD/Align_CTD | **L1/CTD/Align_CTD** | step01 | step00a, step00b, step00d |
| `ctd_celltm` | L2/CTD/Cell_tm | **L1/CTD/Cell_tm** | step01 | step00d |
| `ctd_section` | L2/CTD/Section | **L1/CTD/Section** | step01 | step00d |
| `ctd_loopedit` | L2/CTD/Loop_Edit | **L1/CTD/Loop_Edit** | step01 | step00d |
| `ctd_derived` | L2/CTD/Derived | **L1/CTD/Derived** | step01 | step00d |
| `ctd_1s` | L2/CTD/CTD_1s | **L1/CTD/CTD_1s** | step01 (BinAvg 1 s, full cast) | step02_merge_ctd_suna |
| `suna_l1` | L1/SUNA | L1/SUNA | step01_suna_spectral_qc | step02_merge_ctd_suna, step03, step04_suna_time_lag |
| `suna_l2` | L2/SUNA | L2/SUNA | step02_merge_ctd_suna (merged 1 s), step02_suna_colocation (co-located) | step02b, step02c, step03 (both), step04 (both), step05_profile_diagnostics, step06 |
| `ctd_bottle` | L2/CTD/Bottle | L2/CTD/Bottle | step01 (bottle mode) | build_bottle_nitrate_table |
| `nitrate_qc` (`*_nitrate_QC.csv`) | L3 | **L2/nitrate** | step03_nitrate_qc | baseline |
| `suna_calibrated` | L3/SUNA | **L2/SUNA_calibrated** | step03_suna_bottle_calibration, step02b (only if its correction works) | step04_plot_profiles, step06 |
| `chl_calibrated` | L3/productivity | **L2/productivity** | step07_productivity_calibration | baseline |
| `chl_calibrated_noNPQ` | L3/productivity_noNPQ | **L2/productivity_noNPQ** | step07_noremoval | baseline |
| `ctd_grid_1m_down` | L2/CTD/CTD_1m_down | **L3/CTD/CTD_1m_down** | step01 (BinAvg 1 m, down-cast) | step00d, step04_plot_profiles, step05_profile_diagnostics, step06, step07 (all three) |
| `grid_1m_down` (`<cast>_L3_1m_down.csv`) | L3/gridded | L3/gridded | step08_grid_l3 (`suncat.gridding`) | step04 to step06 after phase 4 |
| `qc_plots` | L2/CTD/qc_plots | **_figures/qc** | step00a, step00b, step00c, step00d | people |
| `profile_plots` | L2/CTD/profile_plots | **_figures/profiles** | step04_plot_profiles, step05_profile_diagnostics | people |
| `profiles_table` (`combined_all_casts_*.csv`) | L2/CTD | **L3/CTD** | combine_readable_profiles | people |

## Findings from the inventory (to resolve in phase 4)

- **Which nitrate feeds the figures.** The README names `step03_nitrate_qc` as the production nitrate and
  `step03_suna_bottle_calibration` as a cross-check, but `step04_plot_profiles` and `step06` read the
  cross-check product (`suna_calibrated`), and nothing downstream reads `nitrate_qc`. The L3 grid should take
  its nitrate from the production product, and the plots should then read the L3 grid.
- **chl-a is calibrated on 1 m bins today.** `step07` fits and applies the bottle calibration on the 1 m
  down-cast grid, so `chl_calibrated` is binned, not native resolution. Until phase 4 applies the fit at
  native resolution and grids afterwards, `L2/productivity` will hold a binned product. The fit itself does
  not need to change.
- **`step05_productivity_calibration.ipynb`** writes to the same folder as `step07_productivity_calibration`
  and appears to be an earlier copy. Confirm it is retired so two notebooks cannot overwrite one product.

## Migration, in phases

Each phase is its own pull request and keeps the regression baseline green unless the change is intended.

1. **Define** (this document, the registry in `project_paths.py`, `tests/test_layout.py`, the `_figures`
   rule in `.gitignore`). No step changes; `ACTIVE_LAYOUT = "v1"`.
2. **Centralize.** Each notebook and script asks `product_dir(CRUISE_ID, "<key>")` for its folders instead
   of building `"L2"`/`"L3"` paths. Under v1 every name resolves to today's folder (pinned by
   `test_v1_is_exactly_todays_folders`), so outputs are unchanged and the baseline stays green.
3. **Flip**, in two pull requests. (a) `tools/migrate_layout.py` and the cruise-layout guard in
   `product_dir`, with `ACTIVE_LAYOUT` still `"v1"`; then migrate the data (below). (b) `ACTIVE_LAYOUT = "v2"`.
   Re-running every step then reproduces byte-identical products in their new places, which the baseline
   confirms.
4. **Science.** The flag-aware Python gridder (calibrated nitrate and chl-a onto the SBE 1 m bins, one L3
   table per cast), step07 applying its calibration at native resolution, and the nitrate source above. These
   change products on purpose; the baseline reports exactly what moved, for review before re-freezing.
5. **Documentation.** README level table and folder tree, `PROCESSING_PARAMETERS.md`, and the slide.

## Migrating a cruise

```
python tools/migrate_layout.py                    # dry run: what moves, what no product claims
python tools/migrate_layout.py --apply            # every cruise to v2 (one cruise: --cruise P45_07)
python tools/migrate_layout.py --to v1 --apply    # roll back
```

Files are assigned to products by the registry above (folder plus migrate glob, longest folder wins). A
cruise is planned in full before anything moves and is refused if any destination already holds different
content. Each file is copied, its SHA-256 checked against the source, renamed into place, and only then is
the source deleted; every move is logged to `_audit/layout_migration_<UTC>.csv` as it happens, so an
interrupted run resumes. Files no product claims are listed and left alone. A cruise that finishes cleanly
gets `_audit/layout.json`, and the frozen baseline (`tests/baseline/local/`) is remapped with its cruise.

`product_dir()` checks that marker (or, without it, which layout's folders hold files) and raises if the code's
layout differs from the cruise's, so a step can never write new outputs beside data in the other layout.
`python project_paths.py` prints each cruise's stored layout.

## Using it in code

```python
from project_paths import product_dir

merged_dir = product_dir(CRUISE_ID, "suna_l2")                 # read
out_dir    = product_dir(CRUISE_ID, "nitrate_qc", create=True)  # write
```

Never build `"L1"`, `"L2"` or `"L3"` paths by hand in a step. To change where a product lives, change the
registry, the test, and this document together.

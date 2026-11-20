# _EXAMPLE: synthetic cruise for testing and review

**This is fabricated data, not a real cruise.** It exists so anyone can run the pipeline end to end without the
confidential P45 measurements. Three casts (E1, E2, E3) on 20 Apr 2026 along a made-up cross-shelf section,
generated with a fixed random seed.

Planted structure, so you know what a correct run should recover:
- a thermocline (surface 29 C to deep 19 C), a deep chlorophyll maximum near 25 m, and a nitracline below 30 m
  rising to about 20 uM at depth;
- a SUNA that reads low on purpose, `SUNA = -1.10 + 0.879 x true nitrate`;
- daytime non-photochemical quenching: surface (< 10 m) fluorescence suppressed to 80%;
- an ECO fluorometer that reads 1.14 x true chlorophyll.

## What to run

step01 (Sea-Bird processing) is skipped: it needs `.hex` files and the SBE software. The example ships step01's
outputs instead (`L1/CTD/*.cnv`, `L1/CTD/Derived/*_der.cnv`, `L1/CTD/CTD_1s/*_1s.cnv`,
`L3/CTD/CTD_1m_down/*_1m_down.cnv`), and its settings (`metadata/suna_cast_map.csv`, `metadata/suna_calibration.json`
in `fit` mode, `metadata/stations.csv`). Run, with `CRUISE_ID = "_EXAMPLE"`:

```bash
python tools/run_notebook.py processing_scripts/notebook/suna/step01_suna_spectral_qc.ipynb --set CRUISE_ID=_EXAMPLE
```

1. `suna/step00_suna_inventory`, `suna/step01_suna_spectral_qc`
2. `step02_merge_ctd_suna`, `suna/step02_suna_colocation`
3. `suna/step03_suna_bottle_calibration`, then `step03_nitrate_qc`
4. `productivity/step07_productivity_calibration`
5. `step08_grid_l3`
6. `step04_plot_profiles`, `step05_profile_diagnostics`, `step06_hydrographic_sections`

`tests/test_example_cruise.py` runs this path on a copy and checks the numbers below. The SUNA `.cal` comes from
the shared `metadata/calibration/` (run from the repository, or copy that folder next to `cruises/` under your
data root). The spectra are not physically calibrated, so keep `STAGE1_MODE = "onboard"`.

## Expected results (a correct run, Nov 2026)

| Step | Result | Planted |
|---|---|---|
| SUNA spectral QC | 543 of 633 frames L1-clean (85.8%) | |
| step02 merge | 3 of 3 casts paired, SUNA coverage of the CTD 1 s scans 100% | |
| step03_suna_bottle_calibration (homogeneous layer) | `SUNA = -1.131 + 0.8788 x bottle`, n 7, R2 0.9999 | -1.10, 0.879 |
| step03_nitrate_qc (Stage 3, fit mode) | a1 0.882, a0 -1.15 (10 bottles); RMSE 1.76 to 0.04 uM | -1.10, 0.879 |
| step07 chlorophyll | gain 1.127 (95% CI 1.124 to 1.131), offset +0.001, n 6; 3 surface bottles left out for quenching | 1.14 |
| step08 | 3 casts gridded; no bin flagged by the down/up check; 7 of 187 chlorophyll bins above the highest bottle | |

The SUNA here has no lag, so the down- and up-casts agree; the real cruises do not (docs/VALIDATION.md).

## Layout notes

- Raw SUNA frames are in `L0/SUNA`; `L1/SUNA` also holds copies from the earlier layout, beside the QC tables
  step01 writes there.
- Discrete tables are filled: `metadata/bottle_nitrate/`, `metadata/chla/`. `metadata/winkler_oxygen/` is empty, so
  the oxygen path stays dormant. There is no PAR channel, so the PAR rule for quenching is not assessed.
- The CTD 1 s files carry no Loop Edit flag column, so step08 falls back to its 1 s record-gap rule and says so.

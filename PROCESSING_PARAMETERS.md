# Processing parameters — record of what the pipeline actually does

**Purpose.** The executable configuration lives in the Sea-Bird PSA files, and those are
**not** tracked in git (they are environment-specific; see `README.md`). This file is the
tracked, human- and machine-readable record of every processing parameter, so the repo
still fully specifies how the data were processed even when the `psa/` folder is absent
from a clone. If you change a PSA, update the matching value here in the same commit.

Cruise: P45_06 / P45_07 - Instrument: SBE 19plus V2, config `19-8460.xmlcon` -
Software: SBE Data Processing v7.26.7.

---

## Chain

```
DatCnv → Filter → Align CTD → Cell Thermal Mass → Section (python) → Loop Edit → Derive → Bin Average
```

`APPLY_WILDEDIT = False`, `APPLY_FILTER = True`, `APPLY_CELLTM = True`, `APPLY_SECTION = True`.
`APPLY_BOTTLES` is per cruise (False for P45_06; True for P45_07 and later cruises with `.afm`).

## Per-module parameters

| Stage | PSA | Parameters |
|-------|-----|------------|
| DatCnv | `01_datcnv.psa` | Output vars: timeS, prDM, t090C, c0S/m, 2x fluorescence, beam attenuation, PAR, SBE 43 O2 [V], SBE 43 O2 [µmol/l]. **O2: Tau correction OFF, hysteresis ON, window 2 s.** Descent-rate window 2 s. Nominal lat/lon 5 / 0. Full cast (down+up), no binning. |
| Filter | `05_filter.psa` | Low-pass A = **1.0 s** on **pressure**; Low-pass B = **0.5 s** on **temperature** and **conductivity**. All other channels: no filter. |
| Align CTD | `02_alignctd.psa` | Temperature **+0.5 s**; conductivity **0 s**; SBE 43 oxygen **+5 s**. (19plus V2 row.) |
| Cell Thermal Mass | `04_celltm.psa` | alpha (thermal anomaly amplitude) **0.03**; 1/beta (time constant) **7.0**. *(OPEN — see below.)* |
| Section | *(none — Python)* | Motion-based descent onset; remove leading surface soak + equilibration dip, keep downcast + upcast. Onset threshold 0.10 db over ~1 s; surface-skip 8 db; degeneracy fallback to surface crossing. CRLF preserved; `# nvalues` updated. |
| Loop Edit | `06a_loopedit_percentmean.psa` (default) | **Percent of mean speed**, window **3 s**, **20 %**, minimum velocity **0.10 m/s**, surface soak removal **off**, exclude scans marked bad **on**. |
| Loop Edit (alt) | `06b_loopedit_fixedmin.psa` | **Fixed minimum velocity 0.10 m/s** (rough-weather alternative), otherwise as above. |
| Derive | `07_derive.psa` | Practical salinity, density (sigma-theta); depth / sound velocity / potential temperature as configured. |
| Bin Average (science) | `08b_binavg_1m_down.psa` | 1 m depth bins, **downcast only** (exclude scans marked bad → upcast dropped via Loop Edit flags). |
| Bin Average (SUNA) | `08d_binavg_1s_full.psa` | 1 s time bins, **full cast** (for SUNA time-matching). |

## Bottle mode (APPLY_BOTTLES = True)

- DatCnv uses `01_datcnv_bottles.psa` (auto-derived from `01_datcnv.psa`; **CreateFile = 2** =
  "Create both data and bottle file"; **ScanRangeSource = 0** = Auto-Fire Module `.afm`). It
  produces `<cast>.cnv` and `<cast>.ros`. The file self-heals if a stale `CreateFile != 2` is found.
- `.afm` staged beside the `.hex` by matching basename. Casts with no `.afm` fall back to the
  profile PSA automatically (`.cnv` only).
- Bottle Summary (`10_bottlesum.psa`, if present) → `L2/CTD/Bottle/<cast>.btl`.
- `build_bottle_nitrate_table.py` parses the `.btl` files (+ optional lab nitrate joined on
  `cast_id` + `bottle_no`) into `<CRUISE>_bottle_ctd.csv` and the step03 `<CRUISE>_bottle_nitrate.csv`.

## Runtime / safety

- `RUN_SBE_COMMANDS` — dry run when False (writes command plan only).
- `STOP_ON_ERROR = True` — halt at first failing module (set False to process all casts and
  collect failures in `_audit/02_failed_sbe_processing_steps.csv`).
- `MODULE_TIMEOUT_SECONDS = 600` — kill a hung SBE module rather than orphan it.
- Work-folder cleanup clears read-only, retries, and falls back to a timestamped folder if locked.

## Measured effect of the Filter (P45_06, 11 casts)

Reduction in the std of the scan-to-scan change (dX/dt) after filtering:
conductivity 33 to 54 % (mean 46 %); pressure 12 to 34 % (mean 24 %); temperature ~0 %.
Source: `step00a_filter_qa.py` (`_audit/qc_filter_dpdt_summary.csv`).

---

## Open decisions (tracked, not silently defaulted)

1. **Cell Thermal Mass coefficients.** Using SBE defaults alpha 0.03 / (1/beta) 7.0. Some groups
   use 0.04 / 8.0 for the 19-series. **Confirm with Drew** before treating CTM-sensitive salinity as final.
2. **Oxygen alignment.** Held at +5 s (mid-range); the swept optimum was inconclusive. Final O2
   alignment defers to bottle-O2 calibration, so the O2 product is provisional.
3. **Nitrate calibration coverage.** Zheng (2024) targets ~50 bottles spanning low/mid/high nitrate.
   A single cruise with few bottles is statistically under-constrained; step03 warns. **Label the
   `no3_qc` product provisional** until a cruise has adequate bottle coverage.
4. **Wild Edit / despiking.** Left out of the chain; the T + C low-pass Filter handles the spiking
   that matters for salinity. Revisit only if a specific cast needs it.

## Tests

`tests/test_ctd_lib.py` (run `pytest`) covers the core parsers and seawater functions in
`ctd_lib.py`: `.cnv` parsing, the Section step (row removal + CRLF preservation), PSS-78 salinity,
potential temperature, and pressure→depth, against published reference values.

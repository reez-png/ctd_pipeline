# Processing parameters: a record of what the pipeline actually does

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
`APPLY_BOTTLES` is per cruise (False for P45_06, which has no `.afm`: its Niskins were fired on the move without the
auto-fire log; True for P45_07 and later cruises with `.afm`).

## Instrument configuration (`19-8460.xmlcon`, both cruises)

The SBE 19plus V2 (S/N 8460) runs in profiling mode at 4 Hz with three external voltage channels (indices 3 to 5) and a
WET Labs RS-232 serial sensor that supplies three channels (indices 6 to 8). Each sensor's calibration is the one in the `.xmlcon`; the table says which `.cnv` column it becomes.

| `.xmlcon` sensor index | Sensor | Serial | Calibration in the `.xmlcon` | `.cnv` column |
|---|---|---|---|---|
| 0, 1 | temperature, conductivity | 8460 | 03-Aug-25 | `tv290C`, `c0S/m` |
| 2 | strain-gauge pressure | 8460 | 04-Sep-25 | `prdM` |
| 3 | SBE 43 dissolved oxygen | 4704 | 08-Jul-25 (Sea-Bird 2007 equation) | `sbeox0V`, `sbeox0Mm/L`, `sbox0Mm/Kg` |
| 4 | Satlantic log PAR | 2419 | 28 Jul 2025 | `par/sat/log` |
| 5 | WET Labs C-Star, 25 cm path | 2655 | 10 Sep 2025 | `CStarAt0` |
| 6 (serial) | ECO V2 chlorophyll, 470 nm excitation | 10043 | scale 7.345e-3, dark 54 (sheet of 8 Sep 2025) | `flECO-AFL` |
| 7 (serial) | ECO V2 backscatter, 700 nm | 10043 | scale 1.676e-6, dark 55 (beta sheet) | not exported |
| 8 (serial) | ECO V2, most likely chlorophyll at 435 nm excitation | 10043 | scale 10, dark 55: a placeholder, no sheet has it | `flECO-AFL1` (not used) |

**The ECO V2 sent Legacy (ECO Classic) counts.** On P45_07 channel 1 reaches exactly 4130 counts on casts 07 and 11
and channel 3 on casts 07 and 10: the full range of the ECO V2's Legacy output (ECO V2 manual, appendix A.3). In Legacy
mode the dark count is 50 and the scale factors are the ones the sensor's `GetCC` command reports for the Legacy
range set; the `.xmlcon` holds the characterization sheet's instead. The bottle calibration (step07) absorbs any
linear difference, so the calibrated chlorophyll stands, but the `.cnv` fluorescence cannot be taken as the factory scale.
Record `GetCC` whenever Legacy mode is entered or its range changed, and enter those coefficients in the `.xmlcon`.

The SBE 55 ECO water sampler that fires the Niskins is not in the `.xmlcon`; on P45_07 it wrote an `.afm` per cast.

## Per-module parameters

| Stage | PSA | Parameters |
|-------|-----|------------|
| DatCnv | `01_datcnv.psa` | Output vars (as named in the `.cnv`): `timeS`, `prdM` (strain-gauge pressure), `tv290C` (ITS-90), `c0S/m`, `flECO-AFL` and `flECO-AFL1` (ECO channels 1 and 3, see below), `CStarAt0`, `par/sat/log`, `sbeox0V`, `sbeox0Mm/L`. The ECO's channel 2 (700 nm backscatter, `TurbidityMeter` in the `.xmlcon`) is not selected, so it is in the `.hex` but in no `.cnv`. **O2: Tau correction OFF, hysteresis ON, window 2 s.** Descent-rate window 2 s. Nominal lat/lon 5 / 0. Full cast (down+up), no binning. |
| Filter | `05_filter.psa` | Low-pass A = **1.0 s** on **pressure**; Low-pass B = **0.5 s** on **temperature** and **conductivity**. All other channels: no filter. |
| Align CTD | `02_alignctd.psa` | Temperature **+0.5 s**; conductivity **0 s**; SBE 43 oxygen **+5 s**. (19plus V2 row.) |
| Cell Thermal Mass | `04_celltm.psa` | alpha (thermal anomaly amplitude) **0.03**; 1/beta (time constant) **7.0**. *(OPEN, see below.)* |
| Section | *(none: Python)* | Motion-based descent onset; remove leading surface soak + equilibration dip, keep downcast + upcast. Onset threshold 0.10 db over ~1 s; surface-skip 8 db; degeneracy fallback to surface crossing. CRLF preserved; `# nvalues` updated. |
| Loop Edit | `06a_loopedit_percentmean.psa` (default) | **Percent of mean speed**, window **3 s**, **20 %**, minimum velocity **0.10 m/s**, surface soak removal **off**, exclude scans marked bad **on**. |
| Loop Edit (alt) | `06b_loopedit_fixedmin.psa` | **Fixed minimum velocity 0.10 m/s** (rough-weather alternative), otherwise as above. |
| Derive | `07_derive.psa` | `sal00` (practical salinity), `sigma-é00` (sigma-theta), `depSM` (depth, salt water), `potemp090C`, `sbox0Mm/Kg` (oxygen per kg); `svCM` (sound velocity) on P45_06 only. |
| Bin Average (science) | `08b_binavg_1m_down.psa` | 1 dbar **pressure** bins (`binavg_bintype = decibars`; step08's gridder requires it), **downcast only** (exclude scans marked bad → upcast dropped via Loop Edit flags). |
| Bin Average (SUNA) | `08d_binavg_1s_full.psa` | 1 s time bins, **full cast** (for SUNA time-matching). |

## Bottle mode (APPLY_BOTTLES = True)

- DatCnv uses `01_datcnv_bottles.psa` (auto-derived from `01_datcnv.psa`; **CreateFile = 2** =
  "Create both data and bottle file"; **ScanRangeSource = 0** = Auto-Fire Module `.afm`). It
  produces `<cast>.cnv` and `<cast>.ros`. The file self-heals if a stale `CreateFile != 2` is found.
- `.afm` staged beside the `.hex` by matching basename. Casts with no `.afm` fall back to the
  profile PSA automatically (`.cnv` only).
- Bottle Summary (`10_bottlesum.psa`) → `L2/CTD/Bottle/<cast>.btl`. Averaged variables: pressure,
  temperature, conductivity, SBE 43 oxygen (mean + std per bottle firing). Salinity is deliberately
  **not** in the PSA; it is computed downstream (next line) from the averaged C, T, P.
- `build_bottle_nitrate_table.py` parses the `.btl` files (columns read from the `Bottle ... Date`
  title row) into `<CRUISE>_bottle_ctd.csv`. With no `LAB_FILE` it also writes
  `<CRUISE>_lab_nitrate_TEMPLATE.csv`, pre-keyed with every bottle so the lab fills only `nitrate_uM`.
  With a filled `LAB_FILE` (joined on `cast_id` + `bottle_no`) it writes the step03
  `<CRUISE>_bottle_nitrate.csv` plus an `_UNMATCHED.csv` for bottles/lab rows that did not pair.
  If the `.btl` has no salinity column it computes practical salinity (PSS-78) from the averaged
  conductivity/temperature/pressure via `ctd_lib`, tagging each row `salinity_source` = `btl` or `computed`.

## Runtime / safety

- `RUN_SBE_COMMANDS`: dry run when False (writes command plan only).
- `STOP_ON_ERROR = True`: halt at first failing module (set False to process all casts and
  collect failures in `_audit/02_failed_sbe_processing_steps.csv`).
- `MODULE_TIMEOUT_SECONDS = 600`: kill a hung SBE module rather than orphan it.
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

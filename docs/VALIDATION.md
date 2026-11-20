# Validation record

What the pipeline recovers, so a run can be checked against known values. Two layers: a synthetic
example anyone can reproduce against a planted truth, and the development cruise (P45_06) with real
results. All calibration numbers here are from few bottles on a single cruise and are provisional; the
limitations are stated at the end.

## 1. Synthetic example (reproducible by anyone)

Run `cruises/_EXAMPLE` (set `CRUISE_ID = "_EXAMPLE"`, run the SUNA QA, step02, step03, step07). It is
generated with a fixed seed and a planted structure, so a correct run recovers these:

| Quantity | Planted truth | A correct run recovers |
|----------|---------------|------------------------|
| SUNA spectral QC, frames L1-clean | clean, in-water spectra | 80 to 89% per cast |
| SUNA to CTD merge coverage | 2 s SUNA vs 1 s CTD | about 83% of CTD scans |
| SUNA bottle-bias fit (step03) | SUNA = -1.10 + 0.879 x bottle | offset near -1.14, gain near 0.873, R2 about 0.999 (n = 10) |
| Chlorophyll fit (step07) | sensor = 1.14 x true chl (factory over-read) | gain near 1.15, corrected factor near 0.87 |

The synthetic spectra are not physically calibrated, so keep `STAGE1_MODE = "onboard"` (the "tsp"
re-derivation is not meaningful on fabricated spectra). See `cruises/_EXAMPLE/README_EXAMPLE.md`.

## 2. Development cruise P45_06 (real data)

Results observed during the build. Bottle counts are small (7 to 10), so treat offset as reasonably
constrained and gain as provisional.

**SUNA nitrate (step03, production QC):** bottles matched to merged frames by depth.
- Stage 3 bias fit: SUNA = -0.73 + 0.849 x bottle, n = 10 (8 low, 1 mid, 1 high nitrate).
- RMSE vs bottle: 1.62 uM (Stage 1, onboard) to 1.18 uM (Stage 3, bias-corrected).
- Reference: Zheng et al. (2024) report final-QC RMSE 0.34 to 0.78 uM.

**SUNA nitrate (standalone cross-check, step03_suna_bottle_calibration):** depth-matched on the
L1-clean colocated frames.
- Fit: SUNA = -1.12 + 0.879 x bottle, R2 0.992, n = 7 (2 casts).
- RMSE vs bottle: 1.68 uM (raw) to 0.50 uM (corrected). Mean bias 0.00 uM after correction.

The two fits differ (offset -0.73 vs -1.12) because they select frames differently (all merged frames
vs L1-clean colocated frames) with only 7 to 10 bottles. **The L3 product from `step03_nitrate_qc` is
the number to report;** the standalone is an explorer and cross-check.

**SUNA spectral QC (step01_suna):** about 89% of frames L1-clean on P45_06, consistent with the
synthetic example.

**Chlorophyll (step07):** the in-situ ECO fluorometer over-reads relative to extracted chlorophyll-a by
a factory factor near 1.14 (corrected factor near 0.88). R2 is moderate, as expected for a fluorescence
proxy against discrete extractions, and the fit is provisional at this bottle count.

## 3. Processing reproducibility (unit tests, run in CI)

- `tests/test_ctd_lib.py` checks the `.cnv` parser, the Section (soak + dip removal) step, and the
  seawater functions (PSS-78 salinity, potential temperature, pressure to depth) against published
  reference values, so a regression is caught before it reaches data.
- `tests/test_bottle_table.py` checks the bottle-table builder.
- The stage-change audit (`step00d`) rebuilds the SBE chain and confirms every physical stage from
  DatCnv through Derive matches the stored products to a maximum absolute difference of 0.0, with the
  Section step byte-for-byte identical to the pipeline.
- The `test_loopedit_*` scripts are local Sea-Bird integration diagnostics, not unit tests; they run
  only on a machine with the SBE software and are skipped by CI.

## 4. Limitations (stated plainly)

- **Few bottles, one cruise.** The nitrate and chlorophyll calibrations rest on 7 to 10 bottles from
  P45_06. The gain, in particular, leans on the deep, high-nitrate bottles of a single cast. These are
  cruise-specific and provisional until more cruises add bottles.
- **P45_06 bottle depths are targets.** Its Niskins were fired on the up-cast while the rosette was moving, with no
  auto-fire log, so each depth is the sample sheet's target, not where the bottle closed; the CTD record shows no
  hold at any of them. A depth error matters only where nitrate changes with depth: the homogeneous-layer fit,
  which uses only bottles where the down- and up-cast SUNA agree within 2 uM (at least 3 clean frames on each leg),
  gives gain 0.862 and offset -0.52 (13 bottles, 4 above 2 uM), against 0.886 and -0.59 from all 27, so the error
  does not move the calibration measurably; it adds scatter at the nitracline. A few near-surface samples may come from a separate hand-deployed Niskin; the
  mixed-layer bottles all read 0 to 0.07 uM, so that changes no comparison there.
- **P45_07's gain rests on few homogeneous bottles.** Its up-casts read up to 16 uM above its down-casts in the
  nitracline (probably water held in the optical path: the SUNA was mounted vertically), so only 6 of its 61 bottles pass the
  homogeneous-layer test, and only 2 of those are above 2 uM: cast 08 and cast 15, both near 16 uM, which alone
  imply gains of 0.89 and 1.03. The stated 0.89 / -0.39 is that fit without cast 15 (5 bottles); with it the fit is
  0.957 / -0.53. The 0.89 is supported instead by step02b's fits to 52 bottles (lag-free 0.878, 95% CI 0.814 to
  0.941; lag-corrected 0.908, 0.843 to 0.974; down-cast only 0.837, 0.764 to 0.909) and by P45_06 (same SUNA,
  0.862 to 0.886). Cast 30's 79 m bottle is flagged: 12.8 uM where the SUNA reads about 17.4 uM in uniform water
  and station P3 at 80 m has 18.1, close to the same cast's 30 m bottle (12.4), so probably water from about 30 m.
- **The gridded nitrate is the down-cast, and the down-cast lags in the nitracline.** step08 compares the two
  legs at every 1 m bin (`nitrate_updown_uM`, flagged in `nitrate_hyst_flag` where they differ by more than
  2 uM). On P45_07 1029 of 2480 assessed bins are flagged (41%), on P45_06 176 of 582 (30%), almost all in the
  nitracline (about 8 to 40 dbar on P45_07's shelf casts). The slow casts are among the most flagged (P45_06
  casts 07 and 05: 91% and 71%; P45_07 casts 01, 10, 14 and 09: 59 to 77%), but casts with an ordinary response
  are flagged too wherever the nitracline is sharp (P45_07 cast 12: 74%): with a 1 uM per metre gradient and the
  package at 0.5 to 1 m/s, a 6 s response already puts the legs several uM apart. Where nitrate rises with depth
  and the lag is first-order, a flagged bin's down-cast value reads low by about half the difference; where water
  was held in the optical path (P45_07 cast 14) the split is unknown. The bottles are the nitracline of record.
- **Surface breaks and held water (an observation, three casts).** step01 finds three casts where the package
  came back to the surface for 151 to 202 s inside the CTD's wet interval: P45_07 casts 14 and 20, P45_06 cast 07.
  Two of them are among the three casts whose SUNA behaved as if water was held in its optical path (P45_07 cast
  14; P45_06 cast 07, with cast 05 the third, which had no break), and cast 20 was not. Taking the SUNA out of the
  water and back in mid-cast could leave air or surface water in the open path, but three casts cannot show that.
- **Onboard nitrate, not re-derived.** The onboard nitrate was computed with `SNA2503E.cal`, which is
  not currently in hand. The bottle correction absorbs the mean temperature and salinity bias, so the
  corrected product is usable, but a physics-based Sakamoto TCSS re-derivation is deferred until that
  cal is retrieved.
- **Fluorescence is a proxy.** The ECO fluorometer measures in-vivo chlorophyll fluorescence, a proxy
  for chlorophyll-a, not the pigment directly; step07 calibrates it against extracted samples, and the
  moderate R2 reflects that proxy relationship.
- **Daytime NPQ.** Surface daytime fluorescence is suppressed by non-photochemical quenching; step07
  flags and (by default) excludes those points. The sensitivity to that choice is small on P45_06
  (see step07_npq_comparison).

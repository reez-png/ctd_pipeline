# Validation record

What the pipeline recovers, so a run can be checked against known values. Three layers: synthetic data with a
planted truth that anyone can run, the two real cruises (P45_06, April 2026; P45_07) with the numbers the
pipeline gives on them, and the tests. Every calibration here rests on a modest number of bottles from one
shelf; the limitations are stated in section 4. Numbers are as of 20 Nov 2026.

## 1. Synthetic data (reproducible by anyone)

**`cruises/_EXAMPLE`** (three casts; the SUNA is planted to read `SUNA = -1.10 + 0.879 x true nitrate`, the
fluorometer `1.14 x true chlorophyll`). The whole path after step01 runs on it (`tests/test_example_cruise.py`):

| Notebook | A correct run gives |
|---|---|
| `suna/step00_suna_inventory` | 633 frames in 3 files; 90 warm-up frames flagged |
| `suna/step01_suna_spectral_qc` | 543 of 633 frames L1-clean (85.8%) |
| `step02_merge_ctd_suna` | 3 of 3 casts paired |
| `suna/step02_suna_colocation` | 633 of 633 frames matched to casts; clock check 1.0 s |
| `suna/step03_suna_bottle_calibration` | 7 homogeneous bottles of 10 matched; `SUNA = -1.131 + 0.8788 x bottle`, R2 0.9999 |
| `step03_nitrate_qc` (fit mode) | a1 0.882, a0 -1.15 (10 bottles); RMSE 1.76 to 0.04 uM |
| `productivity/step07_productivity_calibration` | gain 1.127 (95% CI 1.124 to 1.131), n 6 |
| `step08_grid_l3` | 3 casts gridded; no bin flagged by the down/up check (the example's SUNA has no lag) |

The fits recover the planted SUNA gain to within 0.3% and offset to within 0.06 uM, and the chlorophyll gain to
within 1.1%. Until Nov 2026 the example ran only the SUNA QA path (step02 found no SUNA files on a case-sensitive
file system; no `suna_calibration.json`, `Derived` files or bin-type header).

**`tools/make_synthetic_cruise.py`** writes cruises built to differ from P45 wherever the code could assume it
(another cruise id, three-digit cast files, `prDM`, a cast without Loop Edit flags, mean-pressure bins, the
southern hemisphere east of Greenwich, the equator). `tests/test_synthetic_cruise.py` runs step07, step08 and
step06 on them, and checks that the down/up nitrate check raises no false flags.

The synthetic spectra are not physically calibrated, so keep `STAGE1_MODE = "onboard"`.

## 2. The two cruises

Both cruises used the same instruments under the same `19-8460.xmlcon` (`docs/PROCESSING_PARAMETERS.md`).

### SUNA nitrate (step03, Stage 3 of Zheng et al. 2024: `N_QC = (N - a0) / a1`)

| | P45_06 | P45_07 |
|---|---|---|
| Calibration mode (`metadata/suna_calibration.json`) | fit, from this cruise's bottles | adopted |
| a1 (gain), a0 (offset) | 0.886, -0.59 uM (27 bottles; S2 20 m flagged) | 0.89, -0.39 uM |
| Homogeneous-layer fit (bottles where down- and up-cast agree within 2 uM) | 0.862, -0.52 (13 bottles, 4 above 2 uM) | 0.891, -0.40 without cast 15 (5 bottles); 0.957, -0.53 with it |
| step02b (lag-aware) | down-cast fit 0.839 | 52 bottles: lag-free 0.878 (95% CI 0.814 to 0.941), lag-corrected 0.908 (0.843 to 0.974) |
| Bottle depths | sample-sheet targets (fired on the move, no `.afm`) | fired depths from the `.afm` firing record |
| Flagged bottles | S2 20 m | cast 30, 79 m (probably water from about 30 m) |

**Why the modes differ.** P45_06's down- and up-casts agree closely enough that a fit through all its bottles is
sound, and the homogeneous-layer check moves it by less than its own uncertainty. On P45_07 the up-casts read up to
16 uM above the down-casts in the nitracline, so a depth-matched fit through all bottles mixes lag into the gain;
its calibration is stated instead, and supported by step02b and by P45_06 (same SUNA).

### Chlorophyll fluorescence (step07, `.cnv value = offset + gain x extracted chl-a`)

| | P45_06 | P45_07 |
|---|---|---|
| Fit (daylight NPQ rule, the adopted fit) | gain 1.123 (95% CI 0.695 to 1.550), offset -0.002, R2 0.671, n 15 | gain 0.515 (0.454 to 0.576), offset +0.189, R2 0.945, n 18 |
| RMSE against bottles, `.cnv` value to calibrated | 0.217 to 0.186 mg/m3 | 0.620 to 0.258 mg/m3 |
| NPQ-flagged points kept (step07_noremoval) | gain 1.168, R2 0.713, n 17 | |
| PAR rule instead (PAR at least 100 umol photons m-2 s-1 at the bottle) | 5 bottles change; gain 1.031 (0.501 to 1.561), R2 0.652, n 12 | 2 bottles change; gain 0.509 (0.454 to 0.563), R2 0.960, n 16 |
| Highest bottle in the fit (`chl_cal_max`) | 1.2 mg/m3 | 3.63 mg/m3 |
| Kept down-cast scans above it | 2 of 7060 (cast 08) | 215 of 11249 |
| Gridded bins above it (step08) | 0 of 931 | 38 of 1431 |
| Fluorometer data | all 11 casts | casts 01 to 12; none from cast 13 (section 4) |

### Gridding (step08)

| | P45_06 | P45_07 |
|---|---|---|
| Nitrate bins where down- and up-cast differ by more than 2 uM | 176 of 582 assessed (30%) | 1029 of 2480 (41%) |
| NPQ risk bins: daylight rule, PAR rule; bins agreeing | 42, 136; 811 of 931 (87%) | 184, 248; 2785 of 2945 (95%) |

## 3. Processing reproducibility (tests, run in CI)

- The pytest suite (245 passed and 1 skipped on a machine with the cruise data, Nov 2026) covers the `.cnv` reader,
  sectioning, the seawater functions against published reference values, bottle pairing and stop detection, the
  SUNA spectral QC (surface breaks, per-cast windows), water entry, the homogeneous-layer calibration, the gridder
  (down/up check, chlorophyll above range, NPQ by PAR) and the data-level layout. GitHub Actions runs it on Python
  3.10 to 3.12 on every push.
- `tests/test_baseline_products.py` compares the L2 products with a frozen copy (`tools/make_baseline.py`,
  git-ignored), so a change that moves a number is seen before it reaches data.
- The stage-change audit (`step00d`) rebuilds the SBE chain and confirms every physical stage from DatCnv
  through Derive matches the stored products to a maximum absolute difference of 0.0, with the Section step
  byte-for-byte identical to step01.
- The `test_loopedit_*` scripts are local Sea-Bird integration diagnostics; they need the SBE software and are
  skipped in CI.

## 4. Limitations (stated plainly)

**Nitrate**

- **P45_06 bottle depths are targets.** Its Niskins were fired on the up-cast while the rosette was moving, with no
  auto-fire log, so each depth is the sample sheet's target, not where the bottle closed; the CTD record shows no
  hold at any of them. A depth error matters only where nitrate changes with depth: the homogeneous-layer fit,
  which uses only bottles where the down- and up-cast SUNA agree within 2 uM (at least 3 clean frames on each leg),
  gives gain 0.862 and offset -0.52 (13 bottles, 4 above 2 uM), against 0.886 and -0.59 from all 27, so the error
  does not move the calibration measurably; it adds scatter at the nitracline. A few near-surface samples may come
  from a separate hand-deployed Niskin; the mixed-layer bottles all read 0 to 0.07 uM, so that changes no
  comparison there.
- **P45_07's gain rests on few homogeneous bottles.** Its up-casts read up to 16 uM above its down-casts in the
  nitracline (probably water held in the optical path: the SUNA was mounted vertically), so only 6 of its 61
  bottles pass the homogeneous-layer test, and only 2 of those are above 2 uM: cast 08 and cast 15, both near
  16 uM, which alone imply gains of 0.89 and 1.03. The stated 0.89 / -0.39 is that fit without cast 15 (5 bottles);
  with it the fit is 0.957 / -0.53. The 0.89 is supported instead by step02b's fits to 52 bottles (lag-free 0.878,
  95% CI 0.814 to 0.941; lag-corrected 0.908, 0.843 to 0.974; down-cast only 0.837, 0.764 to 0.909) and by P45_06
  (same SUNA, 0.862 to 0.886). Cast 30's 79 m bottle is flagged: 12.8 uM where the SUNA reads about 17.4 uM in
  uniform water and station P3 at 80 m has 18.1, close to the same cast's 30 m bottle (12.4), so probably water
  from about 30 m.
- **The field gain disagrees with the factory check (open).** Sea-Bird's saltwater verification of the SUNA
  (20 Aug 2025) read -0.03 uM in pure water, 1.79 uM in nitrate-free seawater and 41.83 uM in a 40 uM standard: a
  gain of (41.83 - 1.79) / 40 = 1.00. Against the bottles the SUNA reads 0.86 to 0.89 on both cruises, about 11%
  low. Possible causes, not yet separated: the bottle analysis reading high, the onboard algorithm's fixed
  temperature and salinity, or the onboard nitrate having been computed with `SNA2503E.cal` rather than the
  verified `SNA2503C.cal`. Whether the lab analysis was run against a certified reference material is being
  checked. The bottle-corrected product is on the lab's scale either way.
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
- **Onboard nitrate, not re-derived.** The onboard nitrate was computed with `SNA2503E.cal`, which is not in hand.
  The bottle correction absorbs the mean temperature and salinity bias, so the corrected product is usable, but a
  Sakamoto TCSS re-derivation is deferred until that file is retrieved.

**Chlorophyll**

- **The `.cnv` fluorescence is not on the factory scale, and its gain differs 2.2-fold between cruises.** The ECO
  V2 sent Legacy (ECO Classic) counts (it reaches exactly 4130 counts, the Legacy full range, on P45_07), whose
  coefficients come from the sensor's `GetCC`, while the `.xmlcon` holds the characterization sheet's
  (`docs/PROCESSING_PARAMETERS.md`). The same sensor under the same `.xmlcon` gives gain 1.12 on P45_06 and 0.515
  on P45_07. A Legacy range changed between the cruises would do this; so could a different phytoplankton
  community. Without the `GetCC` records the two cannot be separated. Each cruise's bottle calibration stands on
  its own, but the gains are comparable neither with each other nor with the factory calibration (Roesler et al.
  2017 found factory values about 2 times extracted chlorophyll globally).
- **On P45_06 the upper-layer bottles read low against the deeper ones (unexplained).** Among the bottles in the
  fit, the sensor reads 0.21 to 0.65 times the bottle value at 20 to 25 m and 0.87 to 1.96 times at 30 to 55 m
  (cast 08 at 20 m: 0.044 against 0.16 mg/m3). Quenching under bright light explains part at most: the PAR rule
  drops four of these bottles (PAR 104 to 144), but casts 01, 02 and 10 read low at 20 to 25 m with less light, and
  cast 01's 5 m bottle reads 0.2 times the bottle with PAR at 0. Lower fluorescence per unit chlorophyll in the
  upper-layer community, or bottle-depth error in a steep gradient (the depths are targets), would also do it. The
  single line mostly joins the two groups, so it probably under-corrects the upper layer.
- **NPQ: the daylight rule is kept.** In the profiles the two rules disagree (P45_07: PAR at least 100 reaches 15
  to 30 dbar on bright casts, below the daylight rule's 10 dbar, and never on dim casts 02, 04, 05, 07, 12, 13 and
  26). For the bottles it barely matters: the PAR rule moves P45_07's gain from 0.515 to 0.509, and P45_06's from
  1.123 to 1.031 with confidence intervals that overlap almost entirely and no gain in R2. step08 carries both
  flags per bin (`chl_npq_risk`, `chl_npq_risk_par`).
- **P45_07 lost the ECO from cast 13.** From cast 13 to 31 both ECO channels read exactly 0 counts: the CTD
  received nothing from the ECO, so 19 casts have no chlorophyll and the bottles of 12 of them cannot enter the fit.
- **The third ECO channel is not used.** `flECO-AFL1` is most likely the 435 nm chlorophyll channel, entered in the
  `.xmlcon` with a placeholder scale (10) and dark count (55). With particles (beam attenuation) held fixed, it
  follows the 470 nm chlorophyll on all 11 P45_06 casts (median partial rank correlation 0.60, against 0.14 for beam
  attenuation with chlorophyll held fixed) and on 7 of the 11 P45_07 casts with a usable record
  (`tools/inspect_eco_channels.py`). It sometimes reads below the Legacy dark level of 50 counts (down to 20),
  which is unexplained.
- **Above the calibration's range.** Where the corrected value exceeds the highest bottle in the fit
  (`chl_above_cal`), the line is extrapolated; on P45_07 that is 215 kept down-cast scans and 38 gridded bins.
- **Fluorescence is a proxy.** The ECO measures in-vivo chlorophyll fluorescence, not the pigment; the reference
  is extracted fluorometric chlorophyll, not HPLC, and can read high with chlorophyll b or phaeopigments present
  (Welschmeyer 1994).

**Oxygen**

- The SBE 43 has no Winkler calibration yet (the oxygen path in step07 is dormant), and its alignment (+5 s) is
  provisional until it has one.

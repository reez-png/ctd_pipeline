# _EXAMPLE: synthetic cruise for testing and review

**This is fabricated data, not a real cruise.** It exists so anyone can run the pipeline end to end
without the confidential P45 measurements. Three casts (E1, E2, E3) along a made-up cross-shelf
section, generated with a fixed random seed.

Planted structure, so you know what a correct run should recover:
- thermocline (surface 29 C to deep 19 C), a deep chlorophyll maximum near 25 m, and a nitracline
  below 30 m rising to about 20 uM at depth;
- a SUNA that reads low on purpose, `SUNA = -1.10 + 0.879 x true_nitrate`, so step03 should recover
  roughly that offset and gain;
- daytime non-photochemical quenching: surface (< 10 m) fluorescence suppressed to 80%;
- an ECO fluorometer that over-reads by a factory factor of 1.14, so step07 should recover a gain
  near 1.15 (corrected factor about 0.87).

## What to run

step01 (Sea-Bird processing) is **skipped**: it needs the binary `.hex` and the SBE software. The
example ships step01's outputs instead (`L1/CTD`, `L2/CTD/CTD_1s`, `L2/CTD/CTD_1m_down`), so you start
at the SUNA QA and the merge. In each notebook set `CRUISE_ID = "_EXAMPLE"`, then run, in order:

1. `suna/step00_suna_inventory` and `suna/step01_suna_spectral_qc`  (reads `L0/SUNA`, writes the
   `suna_L1_clean` flags into `L1/SUNA`)
2. `step02_merge_ctd_suna`      -> `L2/SUNA/*_SUNA_1s.csv`
3. `step03_nitrate_qc`          -> `L3/*_nitrate_QC.csv`
4. `step04_plot_profiles`, `step05_profile_diagnostics`, `step06_hydrographic_sections`
5. `productivity/step07_productivity_calibration`

The SUNA `.cal` is resolved from the shared `metadata/calibration/SNA2503C.cal` (used only for the
wavelength grid). The fabricated spectra are not physically calibrated, so keep `STAGE1_MODE = "onboard"`
in step03 (the default); the "tsp" re-derivation is not meaningful on synthetic spectra.

## Expected results (a correct run)

- SUNA spectral QC: roughly 80 to 89% of frames L1-clean per cast.
- step03 bottle-bias fit: offset near -1.1 uM, gain near 0.88, R2 about 0.99.
- step07 chlorophyll fit: gain near 1.15, corrected (factory-to-true) factor near 0.87.
- SUNA coverage of the CTD scans: about 80%.

## Layout notes

- Raw SUNA frames are placed in both `L0/SUNA` and `L1/SUNA`: the SUNA QA notebooks read `L0/SUNA`,
  while step02/step03 read `L1/SUNA` (on a real cruise, `L1/SUNA` is the UCI export of `L0`).
- Discrete tables are filled: `metadata/bottle_nitrate/`, `metadata/chla/`, and `metadata/station_coords.csv`
  (for step06). `metadata/winkler_oxygen/` is empty, the oxygen path stays dormant.

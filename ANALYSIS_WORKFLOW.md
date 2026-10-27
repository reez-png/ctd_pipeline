# Analysis Workflow and FAQ

A task-oriented companion to `README.md`. The README is the reference manual (every parameter, every
folder); this file answers "I have raw files, now what do I actually run, and what comes out?" and
"I only care about one product, what is the shortest path?".

Everything here assumes the same `CRUISE_ID` across the steps you run, and that you have followed the
install and "before the first run" checklist in `README.md`.

---

## 1. The big picture

```text
RAW FILES                     PROCESSING                         PRODUCTS
---------                     ----------                         --------
CTD .hex + .xmlcon  ->  step01 (SBE chain)            ->  L1/CTD/*.cnv  (readable full cast)
  (+ .afm bottles)                                         L2/CTD/.../*.cnv  (processed, 1 m + 1 s)
                                                           L2/CTD/Bottle/*.btl  (bottle mode)

SUNA raw .csv + .cal -> step00_suna, step01_suna (QA)  ->  suna_L1_clean flags (good/bad frames)
                        step02 (merge on UTC time)     ->  L2/SUNA/*_SUNA_1s.csv
                        step03 (Zheng nitrate QC)      ->  L3/<cast>_nitrate_QC.csv  (corrected nitrate)

discrete bottles:
  nitrate table       -> used by step03 (bottle bias)  ->  calibrated nitrate
  extracted chl-a     -> step07 (chl calibration)      ->  L3/productivity/*  (calibrated chlorophyll)
  Winkler oxygen      -> step07 (oxygen, when present)  ->  L3/productivity/*  (calibrated oxygen)

from the products:    step04 (profiles), step05 (diagnostics), step06 (sections)
```

Two rules worth repeating: **Data Conversion (step01) writes L1; everything after writes L2 or L3**, and
**step04 onward reads the products, not the raw files**, so the order matters.

---

## 2. Stage by stage: what it needs, what it makes

| Step | You run | It needs | It produces |
|------|---------|----------|-------------|
| **step01** | `step01_run_sbe_processing_v2.py` (PyCharm, press Run) | `L0/CTD/*.hex`, `*.xmlcon`, the `psa/` chain; `.afm` for bottle mode | `L1/CTD/*.cnv`, `L2/CTD/<stage>/*.cnv` (incl. `CTD_1m_down`, `CTD_1s`), `.btl` in bottle mode, audit CSVs |
| **step00_suna_inventory** | notebook | raw `L1/SUNA/*.csv`, the SUNA `.cal` | frame inventory + health summary in `_audit/suna/` |
| **step01_suna_spectral_qc** | notebook | raw `L1/SUNA/*.csv` | `suna_L1_clean` per-frame flags (good vs bad spectra), QC summary in `_audit/suna/` |
| **step02** | `step02_merge_ctd_suna.ipynb` | `L2/CTD/CTD_1s/*_1s.cnv`, `L1/SUNA/*.csv`, `suna_cast_map.csv`, the SUNA QC flags | `L2/SUNA/<cast>_SUNA_1s.csv` (CTD + SUNA on one 1 s UTC time base) |
| **step03** | `step03_nitrate_qc.ipynb` | `L2/SUNA/<cast>_SUNA_1s.csv`, raw SUNA, the SUNA `.cal`, **bottle nitrate table** | `L3/<cast>_nitrate_QC.csv` (corrected nitrate), assessment tables in `_audit/` |
| **step03_suna_bottle_calibration** | notebook (optional) | same bottle nitrate inputs | interactive SUNA bottle-bias calibration inspector (checks the step03 fit) |
| **step04** | `step04_plot_profiles.ipynb` | `L2/CTD/CTD_1m_down/*.cnv` (+ `L3` for corrected nitrate) | per-cast and multi-cast profile `.png`, `<cast>_readable.csv` |
| **step05** | `step05_profile_diagnostics.ipynb` | the 1 m downcast profiles (+ step02 for nitracline) | annotated multi-panel `.png`, `<cast>_diagnostics.txt`, `profile_diagnostics.csv` |
| **step06** | `step06_hydrographic_sections.ipynb` | the 1 m downcast profiles, a **station coordinate table** (cast to lat/lon) + transect order | contoured section `.png` per variable |
| **step07** | `step07_productivity_calibration.ipynb` | the 1 m downcast CTD (sensor channels), a **discrete chl-a table** (and/or **Winkler oxygen**) | `L3/productivity/*` calibrated product, fit tables + scatter in `_audit/productivity/` |

Discrete sample tables live per cruise under `cruises/<CRUISE_ID>/metadata/`:
- `bottle_nitrate/` for nitrate (drives step03)
- `chla/` for extracted chlorophyll-a (drives step07 chlorophyll)
- `winkler_oxygen/` for Winkler oxygen (drives step07 oxygen, empty until P45_08)
- `station_coords.csv` for the step06 section coordinates

To start a cruise, copy `cruises/_TEMPLATE/` to `cruises/<CRUISE_ID>/` and create its `L0/CTD` and `L0/SUNA` folders. Every step resolves these paths through `project_paths.py`, so the only per-cruise switch is `CRUISE_ID`.

---

## 3. "I just want X": the shortest path

### I just want clean, calibrated nitrate profiles
1. `step01` (produces the 1 s CTD needed for the SUNA merge, and the 1 m downcast for plotting)
2. `step01_suna_spectral_qc` (flags bad SUNA frames; `step00_suna_inventory` first if you want a health check)
3. `step02` (merge)
4. Fill `metadata/bottle_nitrate/` with your lab nitrate, then `step03` (calibration + QC)
5. `step04` for the figures

Without a bottle nitrate table, step03 still runs but `no3_qc` carries **no bottle correction**, so fill
the table before you trust the numbers.

### I just want the productivity (chlorophyll) calibration
This path does **not** need SUNA at all.
1. `step01` (you only need `L2/CTD/CTD_1m_down/*.cnv`, which carries the ECO fluorometer channel)
2. Put your extracted chl-a samples in `metadata/chla/<CRUISE>_chla.csv`
3. `step07_productivity_calibration.ipynb`

You get the fit (`sensor = a + b * bottle`), the corrected chlorophyll product in `L3/productivity/`, and
the fit statistics (offset, gain with 95 % CI, R2, RMSE raw vs corrected) in `_audit/productivity/`.
Use `step07_noremoval` and `step07_npq_comparison` only if you want to see the near-surface
non-photochemical-quenching sensitivity; the default `step07` (NPQ-flagged points excluded) is the one
to report.

### I just want physical profiles or sections (T, S, density, oxygen)
1. `step01`
2. `step04` (profiles) and/or `step05` (diagnostics)
3. `step06` for cross-shelf sections (needs the station coordinate table)

No SUNA, no bottles required.

### I finally have Winkler oxygen (P45_08 onward)
1. Drop the Winkler table into `metadata/winkler_oxygen/`
2. In `step07`, enable the oxygen target (it is coded but dormant until the data exist)
3. Re-run step07; you get the SBE 43 vs Winkler calibration with the same machinery as chlorophyll

---

## 4. FAQ

**Do I need to run every step in order?**
No, but each step reads the previous stage's output, so you cannot skip a dependency. step04 to step07
read products, not raw files, so step01 must run first. step07 and the physical-profile steps do not need
SUNA; the nitrate steps (step02, step03) do.

**Do I need SUNA for the productivity calibration?**
No. step07 calibrates the CTD-mounted sensors (ECO fluorometer, SBE 43) against discrete bottle samples.
It only needs step01's 1 m downcast CTD and a discrete sample table.

**Where do my discrete (bottle) samples go?**
Per cruise, under `cruises/<CRUISE_ID>/metadata/`: nitrate in `bottle_nitrate/`, extracted chlorophyll-a
in `chla/`, Winkler oxygen in `winkler_oxygen/`, and step06 coordinates in `station_coords.csv`. The
column templates (`*_TEMPLATE.csv`) show the expected columns. These data files are git-ignored; only the
templates are tracked.

**What is the minimum I need for a usable nitrate product?**
A bottle nitrate table. Without it, Stages 2 and 3 of the Zheng QC are skipped and `no3_qc` is
uncorrected. The bottles are matched to SUNA frames by depth (within 2 m), not by a single cast timestamp.

**Why does step04 show "Nitrate (corrected)" and not the raw SUNA nitrate?**
step04 reads the corrected nitrate from the L3 QC product (`nitrate_corr_uM`), so the plotted profile is
the bottle-calibrated one. Run step03 before step04 if you want that panel populated.

**My ECO puck has flECO-AFL and flECO-AFL1. Which one is chlorophyll?**
`flECO-AFL` is the live chlorophyll fluorescence channel. `flECO-AFL1` is an unused/dead slot on this
instrument (its values are noise, correlation with the bottles ~0.2), so step07 uses `flECO-AFL`.

**Is the fluorometer measuring chlorophyll-a?**
It measures in-vivo chlorophyll fluorescence, a proxy for chlorophyll-a, not the pigment concentration
directly. That is exactly why step07 calibrates it against extracted (discrete) chlorophyll-a. Report it
as calibrated fluorescence-based chlorophyll, and keep the accuracy hierarchy in mind
(HPLC > extracted fluorometric > in-situ fluorescence).

**What is NPQ and which step07 result do I trust?**
Non-photochemical quenching suppresses daytime near-surface fluorescence, so shallow daytime points read
low. The default `step07` flags those points (shallower than 10 m, between 07:00 and 17:00) and excludes
them from the fit. That is the result to report. `step07_noremoval` keeps them and
`step07_npq_comparison` overlays both, for sensitivity only.

**The SUNA file id does not match my cast id.**
Map them in `cruises/<CRUISE>/metadata/suna_cast_map.csv` (`cast_id,suna_file`). step02 uses it; a
filename fallback and an auto-proposed map exist if you do not fill it.

**Should I run step01 live the first time?**
No. Start with `RUN_SBE_COMMANDS = False` (dry run, writes only the command plan), check
`_audit/01_psa_diagnostic.csv` (every PSA should have an empty `warnings` column), then set it `True`.

**Nothing changes when I re-run a step.**
Check `CRUISE_ID` is the same everywhere, and that the previous stage actually wrote its output (the
tables in Section 2 say what each step expects as input).

**Will any of this commit my cruise data to git?**
No. All raw and processed data, the PSAs, and the figures are git-ignored. Only code, docs, the small
config files (`suna_cast_map.csv`, the calibration references) and the discrete-sample templates are
tracked. See `.gitignore`.

# Analysis workflow and FAQ

A task-oriented companion to `README.md`. The README is the reference (instruments, levels, every step, the
run order); this file answers "I only want one product: what is the shortest path?" and the questions that come
up when running a cruise. Use the same `CRUISE_ID` in every step.

---

## 1. The big picture

```text
RAW FILES                    PROCESSING                                  PRODUCTS
CTD .hex + .xmlcon      ->   step01 (SBE chain)                     ->   L1/CTD stages, L3/CTD/CTD_1m_down
  (+ .afm bottles)                                                       L2/CTD/Bottle/*.btl (bottle mode)
SUNA frames + .cal      ->   SUNA step00, step01 (QA)               ->   L1/SUNA (suna_L1_clean per frame)
                             step02 merge, SUNA step02 co-location  ->   L2/SUNA
                             SUNA step04 (clock), step02b (response)->   metadata: time offsets, slow casts
bottle nitrate          ->   SUNA step03 (homogeneous-layer fit)    ->   metadata/suna_calibration.json
                             step03_nitrate_qc                      ->   L2/nitrate (production nitrate)
extracted chl-a         ->   step07                                 ->   L2/productivity (calibrated chl)
everything above        ->   step08                                 ->   L3/gridded (one table per cast)
from the products:           step04 (profiles), step05 (diagnostics), step06 (sections)
```

The README's section 4 gives the full order with the tool to run after each step.

---

## 2. "I just want X": the shortest path

### Calibrated nitrate profiles

1. `step01` (the 1 s CTD time base and the 1 dbar bins).
2. `suna/step01_suna_spectral_qc` (`suna/step00_suna_inventory` first for a health check).
3. `step02_merge_ctd_suna` and `suna/step02_suna_colocation`.
4. `suna/step04_suna_time_lag`: if the SUNA clock is off, `tools/set_suna_time_offset.py` and rerun step 3.
5. `suna/step02b_suna_response_correction`, then `tools/set_suna_slow_casts.py`.
6. The bottle table (`build_bottle_nitrate_table.py` with `.btl`, `build_bottle_table.py` without).
7. `suna/step03_suna_bottle_calibration`, then state the calibration: `tools/set_suna_calibration.py <CRUISE> fit`,
   or `adopted --a1 ... --a0 ... --source "..."` when the lag makes a fit through all bottles unsafe (P45_07).
8. `step03_nitrate_qc`, then `step08_grid_l3` for the gridded profiles.

step03 stops if `suna_calibration.json` is missing: every cruise must state its calibration.

### The chlorophyll calibration (no SUNA needed)

1. `step01`.
2. Extracted chl-a in `metadata/chla/<CRUISE>_chla.csv` (template columns).
3. `tools/make_stations_table.py <CRUISE> --write` (positions for the solar-time daylight test).
4. `productivity/step07_productivity_calibration`.

You get the fit, every scan calibrated in `L2/productivity/`, and the fit table and scatter in
`_audit/productivity/`. The default step07 (near-surface daytime bottles left out) is the one to report;
`step07_noremoval` and `step07_npq_comparison` show the sensitivity to that choice, and step07 prints the fit under
a PAR rule beside it.

### Physical profiles or sections (T, S, density, oxygen)

1. `step01`.
2. `step04` (profiles) and `step05` (diagnostics).
3. `step06` for sections (needs `stations.csv` or the cruise's cast summary). step08 first if you want nitrate.

### Winkler oxygen arrives

Put the table in `metadata/winkler_oxygen/` and rerun step07: the oxygen target activates when its file exists,
with the same machinery as chlorophyll.

---

## 3. FAQ

**Which nitrate is the one to report?**
`L2/nitrate/<cast>_nitrate_QC.csv` from `step03_nitrate_qc`, gridded by step08 into `L3/gridded/`. The
`L2/SUNA_calibrated` product from `suna/step03_suna_bottle_calibration` is the cross-check (step04 still plots it;
README, known gaps).

**Why do some gridded nitrate bins carry `nitrate_hyst_flag`?**
The SUNA lags, so in the nitracline the down-cast reads low and the up-cast high. step08 compares the two legs at
every bin and flags where they differ by more than 2 uM. The value is the down-cast; treat flagged bins as
uncertain by about half the difference, and use the bottles there.

**Fit or adopted calibration?**
`fit` when the down- and up-casts agree well enough that bottles anywhere in the profile can be used (P45_06).
`adopted` when the lag is large (P45_07): the gain then comes from the homogeneous layer, supported by step02b's
lag-aware fits, and `source` records that.

**Where do my discrete samples go?**
`cruises/<CRUISE_ID>/metadata/`: nitrate in `bottle_nitrate/`, chlorophyll in `chla/`, Winkler oxygen in
`winkler_oxygen/`. The `*_TEMPLATE.csv` files show the columns. A bottle you know is bad gets a reason in the `flag`
column, which leaves it out of the fits; `build_bottle_table.py` carries those flags when it rebuilds the table.

**My bottles have no firing record.**
Without `.afm`/`.btl` the depths are the sample sheet's targets. `build_bottle_table.py` reports where the CTD record
shows the rosette held still against those depths (`STOP_DEPTHS = "report"`), and can use them (`"apply"`) when the
holds are real bottle stops.

**My ECO has `flECO-AFL` and `flECO-AFL1`. Which is chlorophyll?**
`flECO-AFL`, the 470 nm channel, which step07 calibrates. On these cruises `flECO-AFL1` is most likely the 435 nm
chlorophyll channel with a placeholder scale in the `.xmlcon`; `tools/inspect_eco_channels.py` shows how it
behaves. It is not used.

**Is the fluorometer measuring chlorophyll-a?**
It measures in-vivo fluorescence, a proxy. step07 calibrates it against extracted chlorophyll-a, and only the
calibrated value is chlorophyll. Values above the highest bottle in the fit carry `chl_above_cal`.

**What is NPQ, and which flag do I use?**
Non-photochemical quenching lowers daytime fluorescence near the surface. step07 leaves out bottles at 10 m or
shallower between 07 and 17 local solar time. step08 carries that rule per bin (`chl_npq_risk`) and a light-based
one (`chl_npq_risk_par`, PAR at least 100 umol photons m-2 s-1), which is the better guide in a given profile.

**The SUNA file id does not match my cast id.**
Map them in `metadata/suna_cast_map.csv` (`cast_id,suna_file`); step02 proposes a map from time overlap.

**Should I run step01 live the first time?**
No: `RUN_SBE_COMMANDS = False` first, check `_audit/01_psa_diagnostic.csv`, then set it `True`.

**Nothing changes when I rerun a step.**
Check `CRUISE_ID`, that the previous step wrote its output, and that you are not reading a per-cruise setting you
meant to change (the files in `metadata/`; the `tools/set_suna_*.py` scripts show theirs with `--show`).

**Will any of this commit my cruise data?**
No. Data, PSAs, figures, `_audit` and the per-cruise settings are git-ignored; only code, docs, the calibration
references, `suna_cast_map.csv` and the templates are tracked.

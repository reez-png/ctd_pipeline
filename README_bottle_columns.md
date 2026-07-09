# Bottle Nitrate Table — Required Columns

`step03_nitrate_qc.ipynb` calibrates the SUNA against discrete bottle nitrate samples
(Zheng 2024, Stages 2–3). To do that it needs a bottle table with **explicit numeric depth** and a
**cast key**. The qualitative layer labels in the raw WHOI report (`SURF`, `ML`, `DCM`,
`DEEP LAYER`) cannot place a bottle at a metre, so they must be converted to the schema below.

---

## 1. File location and name

```text
ctd_pipeline_v2/ctd/cruises/<CRUISE_ID>/metadata/bottle_nitrate/<CRUISE_ID>_bottle_nitrate.csv
```

Example: `cruises/P45_05/metadata/bottle_nitrate/P45_05_bottle_nitrate.csv`

Running Section 2 of step03 writes a ready-to-fill `<CRUISE_ID>_bottle_nitrate_TEMPLATE.csv` in that
folder. Fill it, save **without** the `_TEMPLATE` suffix, and rerun step03 from Section 6.

**One row per replicate** (do not pre-average replicates — the notebook does the replicate QC).

---

## 2. Columns

| Column | Required? | Type | Description |
|--------|-----------|------|-------------|
| `cruise_id` | **Yes** | text | Cruise of this sample, e.g. `P45_05`. Must equal `CRUISE_ID` in the notebook. Rows from other cruises are ignored. |
| `cast_id` | **Yes** | text | The CTD cast this bottle was fired on, e.g. `P45_05_CTD_03`. **Must match a cast that has a `*_SUNA_1s.csv` in `L2/SUNA/`.** This is how a bottle is paired to SUNA scans. |
| `depth_m` | **Yes** | number | Niskin firing depth in **metres**. This replaces SURF/ML/DCM labels and is the field used to find co-located SUNA scans. |
| `nitrate_uM` | **Yes** | number | Nitrate (NO₃ + NO₂) in **µmol/L (µM)**. Leave **blank** if below detection (put `BDL` in `flag`). |
| `time_utc` | Optional | datetime | UTC time the Niskin fired, ISO format `YYYY-MM-DDTHH:MM:SS`. If present, matching uses time + depth (more accurate). If absent, matching falls back to depth only. |
| `station` | Optional | text | Station label, e.g. `J2`, `P3`, `S1`. Used only as a cross-check / for your own bookkeeping. |
| `replicate_id` | Optional | text/int | Same value for replicates of one sample. Enables the replicate QC (Zheng): replicate groups whose nitrate spread exceeds 0.5 µM are excluded. |
| `flag` | Optional | text | Quality flag. Use `BDL` (below detection), `NES` (not enough sample), etc. Keeps non-numeric notes **out of** `nitrate_uM`. |

---

## 3. Critical rules

1. **`depth_m` must be numeric and real.** No `SURF`/`ML`/`DCM` text. If you only have layer labels,
   read the actual firing depth off the CTD profile / bottle-fire log for that cast and enter the
   metre value.
2. **`cast_id` must match an existing SUNA cast.** If `L2/SUNA/` has `P45_05_CTD_03_SUNA_1s.csv`, the
   `cast_id` must be exactly `P45_05_CTD_03`. A bottle whose `cast_id` has no matching SUNA file is
   silently skipped.
3. **Below-detection values go in `flag`, not `nitrate_uM`.** Writing `<0.015` into the numeric column
   breaks the regression. Leave `nitrate_uM` blank and set `flag = BDL`.
4. **One cruise per table.** Do not mix `P45_05` and `P45_06` rows expecting both to calibrate — the
   bias correction is fit per cruise. The notebook filters to `CRUISE_ID` only.
5. **Replicates as separate rows**, sharing a `replicate_id`. The notebook averages and QCs them.

---

## 4. Example

```csv
cruise_id,cast_id,station,depth_m,time_utc,nitrate_uM,replicate_id,flag
P45_05,P45_05_CTD_01,J1,2,2025-11-12T09:14:00,0.42,J1-2m,
P45_05,P45_05_CTD_01,J1,2,2025-11-12T09:14:00,0.39,J1-2m,
P45_05,P45_05_CTD_01,J1,25,2025-11-12T09:18:00,5.10,J1-25m,
P45_05,P45_05_CTD_01,J1,60,2025-11-12T09:22:00,15.30,J1-60m,
P45_05,P45_05_CTD_02,J2,5,,,J2-5m,BDL
P45_05,P45_05_CTD_02,J2,80,,22.10,J2-80m,
```

Row 5 is below detection: `nitrate_uM` blank, `flag = BDL`.

---

## 5. Mapping the WHOI report to this table

Your raw report (`Marrec__client__2026.xlsx`) encodes station + layer in the sample ID, e.g.
`P45-05-J2-NUT-US DCM#3`. To convert each row:

| Source field | → | Table column |
|--------------|---|--------------|
| `P45-05` in the ID | → | `cruise_id` = `P45_05` |
| `J2` in the ID | → | `station`; and decode to a `cast_id` (your station→cast mapping) |
| `DCM` / `SURF` / `ML` / `DEEP` label | → | look up the **metre depth** for that layer on that cast → `depth_m` |
| `Nitrate+NO2` column | → | `nitrate_uM` (blank + `flag=BDL` if `<0.015`) |

The station→cast mapping (`J1/J2/J3`, `P1/P2/P3`, `S1/S2/S3` → `P45_05_CTD_XX`) is something only you
have from the cruise log — the notebook cannot infer it.

---

## 6. Sourcing depth and time from the SBE 55 auto-fire log (.afm)

If you used the SBE 55 ECO auto-fire sampler, the `.afm` file it produces (one per cast) records, for
each bottle: **bottle number, firing time, and firing pressure in decibars** (the 19plus V2 transmits
pressure in dbar to the auto-fire module). This directly supplies the `depth_m` and `time_utc` columns.

**Recommended path — use SBE Bottle Summary, not a hand parser:**
1. Run the `.afm` files through **SBE Data Processing → Bottle Summary** (the same GUI used for CTD
   processing). It converts the firing records into a calibrated table with pressure, temperature, and
   salinity at each firing, using the instrument config file.
2. Read the firing pressure/depth and time from that output into `depth_m` / `time_utc`.
3. Join the lab nitrate values by **bottle number** to fill `nitrate_uM`.

This avoids reverse-engineering the raw `.afm` hex format and uses Sea-Bird's validated conversion.

Note: the `.afm` does **not** contain nitrate — it is firing geometry only. Nitrate always comes from
the lab, joined by bottle number.

## 7. How the table is used

- **Stage 0 (match-up):** for each bottle, the notebook finds SUNA scans within ±1 m (`MATCH_DEPTH_TOL_M`)
  or ±30 s (`MATCH_TIME_TOL_S`) and averages them, after dropping replicate groups that disagree by
  more than 0.5 µM (`REPLICATE_MAX_DIFF`).
- **Stage 2:** match-ups with `nitrate_uM < 2.5` are used to fit the temperature-residual coefficients.
- **Stage 3:** all match-ups are used to fit the cruise bias (slope a1, intercept a0). The notebook
  **warns** if there are fewer than ~10 bottles or if they do not span the low / mid / high nitrate
  subranges — the single-cruise, few-bottle case is statistically under-constrained (Zheng targets ~50
  bottles across three subranges).

These tolerances are settings near the top of step03 and can be adjusted for your cast geometry.

# Cruise Data Collection Guide — CTD + SUNA Nitrate Pipeline

**Purpose:** a field-ready checklist of exactly what to collect on each cast, where each file goes in
the project, and how to name it, so the data drops straight into the processing pipeline afterward.
This is the collection reference; for processing see `README.md`.

Replace `<CRUISE>` below with the cruise id, e.g. `P45_07`. All paths are under
`ctd_pipeline\cruises\<CRUISE>\`.

---

## 1. Per-cast checklist (do this for EVERY cast)

| # | Collect | From | Save to folder | Name as |
|---|---------|------|----------------|---------|
| 1 | CTD raw data | CTD offload (`.hex`) | `L0\CTD\` | `<CRUISE>_CTD_NN.hex` |
| 2 | Bottle auto-fire log | SBE 55 offload (`.afm`) | `L0\CTD\` (beside the `.hex`) | `<CRUISE>_CTD_NN.afm` |
| 3 | SUNA raw log | SUNA offload (`.csv`, SATSLF frames) | `L0\SUNA\` (working) → later `L1\SUNA\` | keep instrument name, e.g. `A00000NN.CSV` |
| 4 | Water samples for nitrate | Niskin bottles → lab | (physical) | label bottle #, cast, target depth |

`NN` = zero-padded cast number (`01`, `02`, …). **Do the SUNA and `.afm` offload after every cast** —
the SBE 55 and SUNA overwrite memory between deployments, so an un-offloaded cast is lost.

> **Where the `.afm` goes changed.** Put the `.afm` **directly in `L0\CTD\`, next to its `.hex`, with
> the same base name** (`<CRUISE>_CTD_NN.hex` ↔ `<CRUISE>_CTD_NN.afm`). step01's bottle mode
> (`APPLY_BOTTLES = True`) finds the `.afm` by that matching name and runs Data Conversion in
> "Create both data and bottle file" mode automatically. If you also want a separate archive copy, you
> may keep one in `metadata\bottle_fire\`, but the working copy the pipeline reads must be in `L0\CTD\`.
> (On P45_06 there were no bottle firings, so no `.afm`; casts run profile-only.)

---

## 2. Once-per-cruise / once-per-instrument

| Collect | From | Save to | Name as |
|---------|------|---------|---------|
| CTD config file | instrument / manufacturer | `Instrument_Calibration_Files\` | `19-8460.xmlcon` (this instrument) |
| SUNA calibration | instrument / SUNACom | `metadata\calibration\` (shared, not per-cruise) | `SNA2503C.cal` (SN 2503) — must match deployed SUNA serial |
| Lab nitrate results | WHOI Nutrient Analytical Facility | `metadata\bottle_nitrate\` | `<CRUISE>_lab_nitrate.xlsx` (raw lab sheet) |
| Deployment log (see §4) | you, during cruise | `metadata\` | `<CRUISE>_deployment_log.csv` |

---

## 3. Target folder layout after a cruise

```text
cruises\<CRUISE>\
├── Instrument_Calibration_Files\   19-8460.xmlcon
├── L0\
│   ├── CTD\      <CRUISE>_CTD_01.hex ... _NN.hex
│   │             <CRUISE>_CTD_01.afm ... _NN.afm   (beside each .hex, bottle casts only)
│   └── SUNA\     A00000NN.CSV (raw SUNA logs, instrument-named)
├── metadata\
│   ├── <CRUISE>_deployment_log.csv
│   ├── suna_cast_map.csv           (built/confirmed during processing)
│   ├── bottle_fire\               (optional archive of raw .afm; the working copies live in L0\CTD\)
│   └── bottle_nitrate\
│       ├── <CRUISE>_lab_nitrate.xlsx
│       └── <CRUISE>_bottle_nitrate.csv   (filled table for step03)
└── (L1, L2, L3 are created by the pipeline)
```

(`metadata\calibration\SNA2503C.cal` sits at the project level, shared across cruises.)

> **These collection folders are created manually.** Unlike `L1\`, `L2\`, and `L3\` (which the
> pipeline creates automatically), the collection folders — `L0\CTD\`, `L0\SUNA\`,
> `Instrument_Calibration_Files\`, `metadata\`, and `metadata\bottle_nitrate\` — must be created by
> hand when you set up a new cruise. `metadata\bottle_fire\` is optional (an archive only).

---

## 4. The deployment log — the thing that was missing on P45_06

Keep a simple running log during the cruise. This one file removes almost all downstream guesswork
(cast↔SUNA pairing, bottle↔depth). One row per bottle firing is ideal; at minimum one row per cast.

Suggested columns:

```csv
cast_id,station,suna_file,cast_start_utc,bottle_no,fire_time_utc,target_depth_m
P45_07_CTD_01,J1,A0000021.CSV,2026-11-12T09:14:00,1,2026-11-12T09:16:00,2
P45_07_CTD_01,J1,A0000021.CSV,2026-11-12T09:14:00,2,2026-11-12T09:20:00,25
```

- `suna_file` — which raw SUNA log belongs to this cast. This alone makes `suna_cast_map.csv` exact
  instead of inferred from time overlap.
- `bottle_no` + `target_depth_m` — so lab nitrate (keyed by bottle number) can be joined to a depth.

---

## 5. Two mistakes from P45_06 to avoid

1. **SUNA not offloaded per cast.** Only 6 of 11 casts had SUNA data. Offload the SUNA `.csv` right
   after each cast, same as the `.hex`.
2. **Bottle configuration overwritten between casts.** Keep a per-cast record (the deployment log
   above, plus the `.afm` files in `L0\CTD\`) so which bottle fired at which depth is never lost.

---

## 6. After the cruise — handing off to the pipeline

1. Place files per §3 (remember: `.afm` beside each `.hex` in `L0\CTD\`).
2. Process CTD: run `step01` (see `README.md`). For a bottle cruise, set `APPLY_BOTTLES = True` — step01
   then runs Data Conversion in create-both mode and, if `psa\10_bottlesum.psa` is present, runs
   **SBE Bottle Summary** on each `.ros` to write `L2\CTD\Bottle\<cast>.btl` automatically. Casts with
   no `.afm` are processed profile-only.
3. Merge SUNA: run `step02`. It builds/confirms `suna_cast_map.csv` (use the deployment log's
   `suna_file` column to confirm it).
4. Bottles: read firing depth/time from the `.btl` (or run the `.afm` through SBE Bottle Summary by
   hand), join the lab nitrate by bottle number, and fill `<CRUISE>_bottle_nitrate.csv` (schema in
   `README_bottle_columns.md`).
5. QC nitrate: run `step03`. With the bottle table present, Stages 2–3 apply and produce the
   calibrated `no3_qc` product in `L3\`.

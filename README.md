# OMI CTD + SUNA Processing Pipeline (v2)

Cruise-centric processing of Sea-Bird CTD casts and rosette-mounted SUNA nitrate, following the
June 2026 OMI specification (Marrec / Lucas) and the SUNA quality-control method of
Zheng et al. (2024).

This pipeline turns raw instrument files into a quality-controlled nitrate product, organised by
cruise under an `L0 → L1 → L2 → L3` data-level scheme.

---

## Installation

```bash
git clone <your-repo-url>
cd ctd_pipeline_v2/ctd
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
# source .venv/bin/activate
pip install -r requirements.txt
```

For CTD processing (step01) you also need Sea-Bird **SBE Data Processing** installed (Windows) — the
pipeline calls its executables. steps 02–05 are pure Python.

> **Data is not included in this repository.** The `cruises/*/L0..L3/` folders hold instrument data and
> are excluded via `.gitignore`. Only code, PSAs, calibration reference, docs, and small config files are
> tracked. Place your cruise data locally following `CRUISE_DATA_COLLECTION.md`.

---

## 1. What each level means

| Level | Contents | Binning |
|-------|----------|---------|
| **L0** | Raw `.hex` (CTD) and `.bin` / raw `.csv` (SUNA), renamed | none |
| **L1** | Readable **full-cast** `.cnv` (Data Conversion only) and SUNA `.csv` | **none** |
| **L2** | All CTD post-processing + 1 m / 1 s binning, plus the merged CTD+SUNA product | 1 m down, 1 s full |
| **L3** | SUNA quality-controlled nitrate product | — |

The single most important rule: **Data Conversion writes to L1; everything after it writes to L2.**
L1 is readable but unbinned. (In v1 the binned products were wrongly placed under L1 — this is fixed.)

---

## 2. Folder layout (per cruise)

```text
ctd_pipeline_v2/ctd/
├── psa/                              # Sea-Bird PSA setup files (see note below)
│   ├── 01_datcnv.psa                 # REQUIRED — must be created (see Section 6)
│   ├── 02_alignctd.psa
│   ├── 04_celltm.psa                 # optional, gated off by default
│   ├── 06_loopedit.psa
│   ├── 07_derive.psa
│   ├── 08b_binavg_1m_down.psa
│   ├── 08d_binavg_1s_full.psa
│   └── 09_asciiout_semicolon.psa     # optional
├── metadata/
│   ├── calibration/                  # 19-8460.xmlcon (CTD), SNA2503C.cal (SUNA)
│   └── bottle_nitrate/               # bottle template + filled table (see README_bottle_columns.md)
├── processing_scripts/notebook/
│   ├── step01_run_sbe_processing_v2.ipynb
│   ├── step02_merge_ctd_suna.ipynb
│   └── step03_nitrate_qc.ipynb
├── _sbe_work/                        # short, space-free staging area for Sea-Bird CLI
└── cruises/
    └── P45_06/                       # one folder per cruise (05 = Dec 2025, 06 = Apr 2026)
        ├── P45_06_CTD_cast_summary.csv
        ├── Instrument_Calibration_Files/
        ├── L0/
        │   ├── CTD/                  # *.hex (renamed)
        │   └── SUNA/                 # raw SUNA *.csv (SATSLF frames) / *.bin
        ├── L1/
        │   ├── CTD/                  # *.cnv  full cast, readable, NO binning
        │   └── SUNA/                 # *.csv  from SeaBird UCI
        ├── L2/
        │   ├── CTD/
        │   │   ├── Align_CTD/        # *_al.cnv
        │   │   ├── Loop_Edit/        # *_loop.cnv
        │   │   ├── Derived_Parameter/# *_der.cnv   (unbinned hub)
        │   │   ├── CTD_1m_down/      # *_1m_down.cnv  (science profiles)
        │   │   ├── CTD_1s/           # *_1s.cnv       (for SUNA merge)
        │   │   └── profile_plots/    # step04/05 figures, lab CSVs, diagnostics .txt
        │   └── SUNA/                 # *_SUNA_1s.csv  (merged CTD+SUNA @1s)
        ├── L3/                       # *_nitrate_QC.csv  (final product)
        └── _audit/                   # command logs, inventories, QC summaries
```

---

## 3. The five notebooks, in order

The pipeline is five steps. steps 01–03 must run in sequence (each reads the previous stage's output);
steps 04–05 are visualization/diagnostics and run after step02 (they read CTD profiles + merged nitrate,
and do not require step03). `CRUISE_ID` must be the same in all of them.

### step01 — Sea-Bird processing  →  L1 + L2
Runs Data Conversion (to L1) then Align → [Wild Edit] → [Cell Thermal Mass] → [Filter] → Loop Edit →
Derive (to L2), then Bin Average into `CTD_1m_down` (downcast science profiles) and `CTD_1s`
(full cast, for SUNA matching). Bracketed steps are **gated off by default** — see the decisions memo.

> **Two formats are provided for step01:**
> - `step01_run_sbe_processing_v2.py` — **use this in PyCharm Community.** A flat top-to-bottom
>   script; edit the settings near the top and press Run (green ▶). Recommended.
> - `step01_run_sbe_processing_v2.ipynb` — the notebook version (PyCharm Professional / Jupyter).
>
> Both contain identical logic. step02 and step03 remain notebooks (run them via a browser Jupyter
> server in Community: `jupyter notebook` in the PyCharm terminal).

- **Inputs:** `L0/CTD/*.hex`, `Instrument_Calibration_Files/*.xmlcon`, PSA files in `psa/`
- **Outputs:** `L1/CTD/*.cnv`, `L2/CTD/<process>/*.cnv`, audit logs in `_audit/`
- **Safety:** start with `RUN_SBE_COMMANDS = False` (dry run, writes the command plan only), inspect
  the command log, then set `True`.

### step02 — CTD ↔ SUNA merge  →  L2/SUNA
Interpolates the SUNA `.csv` onto the CTD `*_1s.cnv` time base (UTC), producing a 1 s product with
CTD + SUNA columns. Handles SUNA clock offset and refuses to interpolate across data gaps.

- **Inputs:** `L1/SUNA/*.csv`, `L2/CTD/CTD_1s/*_1s.cnv`
- **Output:** `L2/SUNA/<cast>_SUNA_1s.csv`

> **Pairing CTD casts to SUNA files.** SUNA files are usually named with the instrument's own
> sequential id (e.g. `A0000010.CSV`), which does not match the CTD cast id. step02 resolves the
> pairing in three steps:
>
> 1. **Explicit map** (preferred): a CSV at `cruises/<CRUISE_ID>/metadata/suna_cast_map.csv` with two
>    columns, `cast_id,suna_file` (e.g. `P45_06_CTD_01,A0000010.CSV`). This is the auditable record of
>    which SUNA log belongs to which cast.
> 2. **Name match** (fallback): if no map, it tries to match a SUNA filename to the cast id.
> 3. **Auto-propose**: if anything is still unpaired, step02 reads each SUNA file's frame times, overlaps
>    them with each CTD cast's time window, and writes a suggested `suna_cast_map_PROPOSED.csv` with a
>    confidence flag. Review it, fix any low-confidence rows, save it as `suna_cast_map.csv`, and rerun.
>
> The actual data alignment within a paired cast is done by **UTC time**, not filename — the SUNA
> scans are interpolated onto the CTD 1 s timestamps (`SUNA_CLOCK_OFFSET_S` corrects clock drift).

> **The same `suna_cast_map.csv` is used by BOTH step02 and step03.** It lives at
> `cruises/<CRUISE_ID>/metadata/suna_cast_map.csv`. step03 reads it to find each cast's correct raw
> SUNA log (the `A00000NN.CSV` files don't name-match the casts). Create it once (step02 can propose
> it automatically); both notebooks then use it.

### step03 — SUNA nitrate QC  →  L3
Applies the Zheng (2024) post-processing: (Stage 1) starting nitrate, (Stage 2) low-nitrate
temperature residual, (Stage 3) cruise-specific bottle bias, (Stage 4) assessment vs bottles.

- **Inputs:** `L2/SUNA/<cast>_SUNA_1s.csv`, raw `L1/SUNA/*.csv`, `metadata/calibration/SNA2503C.cal`,
  and the **bottle nitrate table** (see `README_bottle_columns.md`)
- **Output:** `L3/<cast>_nitrate_QC.csv` + assessment tables in `_audit/`

### step04 — Profile plots + lab-friendly CSVs  →  L2/CTD/profile_plots
Plots vertical profiles from the 1 m downcast CTD data and, where a cast has SUNA, the merged nitrate.
Produces per-cast single-variable figures (T, S, O₂, fluorescence, beam attenuation, density, nitrate
vs depth) and per-variable multi-cast overlays. Also writes a **lab-friendly CSV** per cast with
readable column names (`Temperature_C`, `Salinity_PSU`, `Nitrate_uM`, …) for sharing with colleagues.

- **Inputs:** `L2/CTD/CTD_1m_down/*_1m_down.cnv`, `L2/SUNA/*_SUNA_1s.csv` (for nitrate, where present)
- **Outputs:** `L2/CTD/profile_plots/<cast>__<variable>.png`, `_overlay__<variable>.png`,
  `<cast>_readable.csv`

### step05 — Profile diagnostics (thermocline / DCM / MLD / nitracline)  →  L2/CTD/profile_plots + _audit
Computes vertical-structure diagnostics and produces annotated multi-panel figures (3×2 grid when
nitrate is present) with shaded interpretation bands and cited feature depths. Diagnostics: thermocline,
halocline, pycnocline, oxycline, DCM (fluorescence max), MLD by temperature (0.2 °C) and density
(0.03 kg/m³) criteria, and the nitracline reported two ways — **onset** (first depth exceeding the
surface minimum by 1 µM, the ecological top of the nitracline; primary, shaded band) and **max gradient**
(steepest dN/dz; reference line). The nitrate panel overlays raw points, a smoothed trend curve, and the
nitracline. A per-cast `.txt` sidecar records every feature depth with its citation and a method
disclaimer.

- **Inputs:** same as step04
- **Outputs:** `L2/CTD/profile_plots/<cast>__annotated.png`, `<cast>_diagnostics.txt`, and a combined
  `_audit/profile_diagnostics.csv`
- **Note:** cline depths are the strongest smoothed gradient (method-dependent); MLD is threshold-based.
  The on-figure disclaimer and per-feature citations reflect that these are automated, single-profile
  diagnostics — first-pass guidance, not definitive layer boundaries.

---

## 4. Does it produce a CSV? Yes — here is exactly what and when

| Notebook | Output file | Folder | Produced when |
|----------|-------------|--------|---------------|
| step01 (.py or .ipynb) | `<cast>_1m_down.cnv`, `<cast>_1s.cnv` | `L2/CTD/...` | `RUN_SBE_COMMANDS = True` and PSAs present |
| step02 | `<cast>_SUNA_1s.csv` | `L2/SUNA/` | step01 produced the 1 s CNV |
| step03 | `<cast>_nitrate_QC.csv` | `L3/` | step02 produced the merged file |
| step04 | `<cast>__<var>.png`, `<cast>_readable.csv` | `L2/CTD/profile_plots/` | step01 (+ step02 for nitrate) done |
| step05 | `<cast>__annotated.png`, `<cast>_diagnostics.txt`, `profile_diagnostics.csv` | `profile_plots/` + `_audit/` | step01 (+ step02 for nitrate) done |

**step03's `_nitrate_QC.csv` has these columns:**

| Column | Meaning |
|--------|---------|
| `utc_time` | UTC timestamp of the SUNA scan |
| `temp_c` | CTD temperature (°C), nearest-time matched |
| `salinity` | CTD practical salinity |
| `no3_onboard` | SUNA firmware nitrate (µM) |
| `no3_stage1` | Stage-1 nitrate (onboard, or TSP if enabled) |
| `no3_stage2` | after low-nitrate temperature residual correction |
| `no3_qc` | **final quality-controlled nitrate** after bottle bias correction |

> Without a bottle table, step03 still writes `no3_onboard` / `no3_stage1`, but `no3_stage2` and
> `no3_qc` will not carry the corrections (Stages 2–3 are skipped). The real accuracy gain comes from
> the bottle calibration, so fill the bottle table before relying on `no3_qc`.

---

## 5. Two decisions still open (see the decisions memo)

1. **Cell Thermal Mass** — gated OFF; recommended ON for the stratified Gulf of Guinea shelf
   (set `APPLY_CELLTM = True`, add `04_celltm.psa`).
2. **Despiking** — Wild Edit and Filter are gated OFF, leaving no despiking in the chain.
   Recommended: re-enable Wild Edit (`APPLY_WILDEDIT = True`) or apply the salinity despike in step02
   (`APPLY_SALINITY_DESPIKE = True`). See `CTD_pipeline_v2_decisions_memo.docx`.

---

## 5b. Running in PyCharm Community Edition

Community has weak notebook support, so step01 is provided as a plain `.py` script.

1. **Open the project** at the top level: `ctd_pipeline_v2`.
2. **Create the virtual environment:** Settings → Project → Python Interpreter → Add Interpreter →
   Add Local Interpreter → Virtualenv → New, location `ctd_pipeline_v2/.venv`, base Python 3.10/3.11.
3. **Install packages** in the PyCharm terminal (prompt must show `(.venv)`):
   ```
   pip install pandas numpy matplotlib openpyxl
   ```
   (Add `jupyter nbformat` only if you also want to run step02/03 notebooks via the browser.)
4. **Run step01:** open `step01_run_sbe_processing_v2.py`, edit the settings block at the top
   (`CRUISE_ID`, `CTD_ROOT`, `SBE_BIN_DIR`, `RUN_SBE_COMMANDS = False` for the dry run), then press Run.
5. **Run step02 / step03 / step04 / step05:** in the terminal run `jupyter notebook`, open the `.ipynb`
   in the browser that launches, and run the cells there.

---

## 6. Before the first run — checklist

- [ ] **Create `01_datcnv.psa`** in `psa/` (export from the SBE GUI: variable selection, full cast
      down+up, no binning). The other four required PSAs are present.
- [ ] Put the CTD config (`19-8460.xmlcon`) in the cruise `Instrument_Calibration_Files/`.
- [ ] Put the SUNA cal (`SNA2503C.cal`) in `metadata/calibration/`.
- [ ] Set the same `CRUISE_ID` in all three notebooks.
- [ ] Confirm the cruise has **both** SUNA data and bottle nitrate for that cruise — do not calibrate
      one cruise's SUNA against another's bottles (bias is fit per cruise).
- [ ] Fill the bottle table per `README_bottle_columns.md`.
- [ ] step01: dry-run first (`RUN_SBE_COMMANDS = False`), inspect `_audit/02_sbe_processing_command_log.csv`,
      then set `True`.

---

## 7. Fixes and lessons from the P45_06 real-data run

Running on real data surfaced several issues that are now fixed. Recorded here so they are not
re-introduced:

- **PSA output paths were hardcoded.** The six PSAs had frozen `OutputDir`/`OutputFile`/`NameAppend`
  values (two even pointed at the old v1 project). Fixed with `fix_psa_output_paths.py`, which blanks
  those fields so the runner's command-line flags control output. Re-run that script if you ever
  re-save a PSA from the SBE GUI.
- **SUNA files are raw SATSLF frames, not tabular exports.** step02's reader now auto-detects the
  286-field `SATSLF` format (parsing year-doy + decimal hour + nitrate) and falls back to tabular
  parsing for clean UCI exports.
- **SUNA↔cast pairing is by time overlap, via a map.** The `A00000NN.CSV` names don't match casts, so
  both step02 and step03 use `suna_cast_map.csv`. A single long SUNA log can serve several casts — the
  time-based merge splits it correctly per cast window.
- **step03 raw-file resolution.** step03 now uses the map to find each cast's raw SUNA file (it
  previously fell back to the alphabetically-first file, which would have processed every cast against
  the wrong SUNA data).
- **19plus V2 specifics.** Pressure is strain-gauge (correct for this frame, not Digiquartz). CTD
  columns are `tv290C` (temperature) and `sal00` (salinity); step03's column matching handles these.

## 8. Known limitation — Stage 1 TSP

step03's `STAGE1_MODE` defaults to `"onboard"`, which uses the SUNA firmware's
temperature/salinity/pressure-corrected nitrate as the starting point. This is a validated approach
(e.g. OOI's NES nitrate dataset works this way) and the bottle-based Stages 2–3 do the main accuracy
correction.

A `"tsp"` mode that re-derives nitrate from the raw 256-channel spectra is included but **marked
work-in-progress**: the least-squares fit over the narrow 217–240 nm window needs the
`OPTICAL_WAVELENGTH_OFFSET` baseline parameterization (Sakamoto 2009; Bio-Argo DAC manual) to be
numerically stable. Until that is implemented and validated, keep `STAGE1_MODE = "onboard"`.

---

## 9. Next-cruise data-collection checklist

The pipeline is proven. What separated P45_06 from a fully calibrated product was **data collection**,
not code. For the next cruise, ensure:

**SUNA**
- [ ] Offload the SUNA `.csv` (raw SATSLF frames) **after every cast**, alongside the `.hex`. On
      P45_06 only 6 of 11 casts had SUNA data because this wasn't done consistently.
- [ ] Keep a simple deployment log: which cast, which SUNA file, start/stop times. This makes the
      `suna_cast_map.csv` trivial to confirm (rather than inferring from time overlap).
- [ ] Confirm the SUNA `.cal` file matches the deployed instrument serial (currently `SNA2503C.cal`
      for SN 2503). It lives in `ctd/metadata/calibration/`.

**Bottles (needed for the calibrated nitrate product — Stages 2–3)**
- [ ] Offload the SBE 55 auto-fire log (`.afm`) after every cast. It records bottle number, firing
      time, and firing pressure (decibars) per bottle.
- [ ] Process the `.afm` files through **SBE Data Processing → Bottle Summary** to get a calibrated
      bottle table (pressure, T, S at each firing). This is the recommended path rather than parsing
      the raw `.afm`.
- [ ] Do **not** let the bottle configuration overwrite between casts — keep a per-cast record of which
      Niskin fired at which target, so lab nitrate can be joined back to firing depth.
- [ ] Get the discrete nitrate values from the lab (WHOI Nutrient Analytical Facility), keyed so each
      value maps to a bottle number / cast / depth.
- [ ] Aim for ~50 bottles across the cruise spanning low / mid / high nitrate (Zheng's guidance); the
      pipeline warns if there are too few or they don't span the range.

**Then**
- [ ] Fill `<CRUISE_ID>_bottle_nitrate.csv` (see `README_bottle_columns.md`) and run step03 — Stages
      2–3 will apply, producing the calibrated `no3_qc` product in L3.

With those in hand, the next cruise is the first fully quality-controlled nitrate dataset.

---

## 10. References

- Sakamoto, C. M., Johnson, K. S., & Coletti, L. J. (2009). Improved algorithm for nitrate from a UV
  spectrophotometer. *Limnol. Oceanogr.: Methods, 7,* 132–143.
- Sakamoto, C. M., et al. (2017). Pressure correction for nitrate computation. *Limnol. Oceanogr.:
  Methods, 15,* 897–902.
- Plant, J. N., et al. (2023). Updated temperature correction for seawater nitrate. *Limnol.
  Oceanogr.: Methods, 21,* 581–593.
- Zheng, B., et al. (2024). Bias-corrected high-resolution vertical nitrate profiles from the CTD
  rosette-mounted SUNA. *Limnol. Oceanogr.: Methods, 22,* 889–902.
- Sea-Bird Scientific. SBE Data Processing / Seasoft module documentation.

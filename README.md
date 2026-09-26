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
cd ctd_pipeline
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
| **L0** | Raw `.hex` (CTD), `.afm` (bottle auto-fire, beside each `.hex`), and `.bin` / raw `.csv` (SUNA) | none |
| **L1** | Readable **full-cast** `.cnv` (Data Conversion only) and SUNA `.csv` | **none** |
| **L2** | All CTD post-processing + 1 m / 1 s binning, plus the merged CTD+SUNA product and bottle products | 1 m down, 1 s full |
| **L3** | SUNA quality-controlled nitrate product | — |

The single most important rule: **Data Conversion writes to L1; everything after it writes to L2.**
L1 is readable but unbinned.

---

## 2. Folder layout (per cruise)

```text
ctd_pipeline/
├── psa/                                # Sea-Bird PSA setup files (see Section 6)
│   ├── 01_datcnv.psa                   # REQUIRED — Data Conversion (profile: .cnv only)
│   ├── 01_datcnv_bottles.psa           # AUTO-GENERATED when bottle mode is on (Create both + AFM)
│   ├── 05_filter.psa                   # Filter: P 1.0 s (A); T + C 0.5 s (B)
│   ├── 02_alignctd.psa                 # Align CTD: T +0.5 s, C 0 s, O2 +5 s
│   ├── 04_celltm.psa                   # Cell Thermal Mass (alpha 0.03, 1/beta 7.0 — confirm w/ Drew)
│   ├── 06a_loopedit_percentmean.psa    # Loop Edit — percent of mean speed, minV 0.10 (default)
│   ├── 06b_loopedit_fixedmin.psa       # Loop Edit — fixed minimum velocity 0.10 (rough-weather alt)
│   ├── 07_derive.psa                   # Derive (salinity, density, ...)
│   ├── 08b_binavg_1m_down.psa          # Bin Average 1 m, downcast only (science)
│   ├── 08d_binavg_1s_full.psa          # Bin Average 1 s, full cast (for SUNA merge)
│   ├── 09_asciiout_semicolon.psa       # optional ASCII export
│   └── 10_bottlesum.psa                # optional — Bottle Summary (.ros -> .btl), bottle mode only
├── metadata/
│   ├── calibration/                    # 19-8460.xmlcon (CTD), SNA2503C.cal (SUNA)
│   └── bottle_nitrate/                 # bottle template + filled table (see README_bottle_columns.md)
├── processing_scripts/notebook/
│   ├── step01_run_sbe_processing_v2.py     # CTD processing (run this in PyCharm Community)
│   ├── step02_merge_ctd_suna.ipynb
│   └── step03_nitrate_qc.ipynb
├── cnv_to_csv.py                       # standalone: convert a folder of .cnv to .csv
├── verify_reproducibility.py           # standalone: compare two runs' .cnv data columns
├── _sbe_work/                          # short, space-free staging area for the Sea-Bird CLI
└── cruises/
    └── P45_06/                         # one folder per cruise (06 = Apr 2026, 07 = ...)
        ├── Instrument_Calibration_Files/
        ├── L0/
        │   ├── CTD/                    # *.hex (renamed) + matching *.afm beside each .hex (bottle casts)
        │   └── SUNA/                   # raw SUNA *.csv (SATSLF frames) / *.bin
        ├── L1/
        │   ├── CTD/                    # *.cnv  full cast, readable, NO binning
        │   └── SUNA/                   # *.csv  from SeaBird UCI
        ├── L2/
        │   ├── CTD/
        │   │   ├── Filter/             # *_filt.cnv
        │   │   ├── Align_CTD/          # *_al.cnv
        │   │   ├── Cell_tm/            # *_ctm.cnv
        │   │   ├── Section/            # *_sec.cnv   (soak + equilibration dip removed)
        │   │   ├── Loop_Edit/          # *_loop.cnv
        │   │   ├── Derived/            # *_der.cnv   (unbinned hub)
        │   │   ├── CTD_1m_down/        # *_1m_down.cnv  (science profiles)
        │   │   ├── CTD_1s/             # *_1s.cnv       (for SUNA merge)
        │   │   ├── Bottle/             # *.btl          (bottle summary, bottle mode only)
        │   │   └── profile_plots/      # step04/05 figures, lab CSVs, diagnostics .txt
        │   └── SUNA/                   # *_SUNA_1s.csv  (merged CTD+SUNA @1s)
        ├── L3/                         # *_nitrate_QC.csv  (final product)
        └── _audit/                     # command logs, inventories, QC summaries
```

> **Note:** `.ros` bottle files (raw bottle-fire scans) are written by Data Conversion beside the L1
> `.cnv`; `Bottle Summary` reads them into `L2/CTD/Bottle/*.btl`. Both are produced only in bottle mode.

---

## 3. The processing chain (step01)

The Amala-approved linear chain, run per cast:

```text
DatCnv → Filter → Align CTD → Cell Thermal Mass → Section → Loop Edit → Derive → Bin Average
```

- **DatCnv** → readable full-cast `.cnv` in `L1/CTD/`.
- **Filter** (before Align): pressure 1.0 s (Filter A); temperature and conductivity 0.5 s (Filter B).
  The T and C low-pass reduces conductivity/salinity spiking ahead of Cell Thermal Mass and Loop Edit.
- **Align CTD**: temperature **+0.5 s**, conductivity **0 s**, SBE 43 oxygen **+5 s** (19plus V2 row).
- **Cell Thermal Mass**: conductivity cell thermal-mass correction (after Align).
- **Section** (a **non-SBE Python step** inside step01): removes the leading surface soak +
  equilibration dip and keeps downcast + upcast, so Loop Edit's running-maximum-pressure rule is not
  poisoned by the dip. Downcast selection is deferred to Bin Average.
- **Loop Edit**: percent of mean speed (window 3 s, 20 %), minimum velocity **0.10 m/s**, surface soak
  removal off. A fixed-minimum-velocity variant (`06b`) is kept for rough weather.
- **Derive**: salinity, density, and derived variables.
- **Bin Average**: 1 m depth bins, **downcast only** (science profiles) and 1 s time bins, full cast
  (for SUNA matching).

Module gates (top of `step01_run_sbe_processing_v2.py`):

```python
APPLY_WILDEDIT = False   # not in the agreed chain
APPLY_CELLTM   = True     # on (after Align)
APPLY_FILTER   = True     # on (before Align)
APPLY_SECTION  = True     # Python soak+dip removal before Loop Edit
APPLY_BOTTLES  = False    # per cruise: True for cruises with .afm bottle files (see Section 3b)
```

> **step01 format:** use `step01_run_sbe_processing_v2.py` in PyCharm Community — a flat top-to-bottom
> script; edit the settings block at the top and press Run (green ▶). step02 and step03 remain
> notebooks (run them via `jupyter notebook` in the PyCharm terminal).

- **Inputs:** `L0/CTD/*.hex`, `Instrument_Calibration_Files/*.xmlcon`, PSA files in `psa/`
- **Outputs:** `L1/CTD/*.cnv`, `L2/CTD/<process>/*.cnv`, audit logs in `_audit/`
- **Safety:** start with `RUN_SBE_COMMANDS = False` (dry run, writes the command plan only), inspect
  `_audit/01_psa_diagnostic.csv` (every PSA should have an empty `warnings` column) and
  `_audit/02_sbe_processing_command_log.csv`, then set `True`.
- **Robustness built in:** the work-folder cleanup clears read-only files and retries (and falls back
  to a fresh timestamped work folder if a leftover is locked); each SBE module has a hard timeout
  (`MODULE_TIMEOUT_SECONDS`, default 600 s) so a stuck module cannot hang the run. `STOP_ON_ERROR`
  halts the batch at the first failure — set it `False` to process every cast and collect all failures
  in `_audit/02_failed_sbe_processing_steps.csv`.

### 3b. Bottle mode (cruises with an SBE auto-fire module)

For cruises where bottles were fired (each cast has an `.afm`), set `APPLY_BOTTLES = True`. Then step01:

1. Runs Data Conversion in **"Create both data and bottle file"** mode with scan-range source =
   **Auto-Fire Module (.afm)**, using an auto-generated `01_datcnv_bottles.psa` (derived from
   `01_datcnv.psa`; it self-heals if stale). This produces `<cast>.cnv` **and** `<cast>.ros`.
2. Stages each cast's `.afm` beside its `.hex` in the work folder (matched by identical basename).
3. **Per cast:** a cast with an `.afm` uses the bottle PSA; a cast **without** one automatically falls
   back to the profile PSA (`.cnv` only) — so a mix of bottle and non-bottle casts runs cleanly.
4. If `psa/10_bottlesum.psa` exists, runs **Bottle Summary** on each `.ros` → `L2/CTD/Bottle/<cast>.btl`.
   If it is absent, step01 still produces `.cnv` + `.ros` and skips only the `.btl` (with a note).

`APPLY_BOTTLES = False` (the default, e.g. P45_06) leaves the profile chain completely unchanged.

---

## 4. The five notebooks / scripts, in order

steps 01–03 run in sequence (each reads the previous stage's output); steps 04–05 are
visualization/diagnostics after step02. `CRUISE_ID` must be the same in all of them.

### step01 — Sea-Bird processing → L1 + L2
The chain in Section 3. Also handles bottle mode (Section 3b).

### step02 — CTD ↔ SUNA merge → L2/SUNA
Interpolates the SUNA `.csv` onto the CTD `*_1s.cnv` time base (UTC), producing a 1 s product with
CTD + SUNA columns. Handles SUNA clock offset and refuses to interpolate across data gaps.

- **Inputs:** `L1/SUNA/*.csv`, `L2/CTD/CTD_1s/*_1s.cnv`
- **Output:** `L2/SUNA/<cast>_SUNA_1s.csv`

> **Pairing CTD casts to SUNA files.** SUNA files carry the instrument's own sequential id (e.g.
> `A0000010.CSV`), which does not match the CTD cast id. step02 resolves the pairing via
> `cruises/<CRUISE_ID>/metadata/suna_cast_map.csv` (`cast_id,suna_file`), a filename fallback, or an
> auto-proposed `suna_cast_map_PROPOSED.csv` from time overlap. Alignment within a paired cast is by
> **UTC time** (`SUNA_CLOCK_OFFSET_S` corrects drift). The same map is used by step02 and step03.

### step03 — SUNA nitrate QC → L3
Applies the Zheng (2024) post-processing: (Stage 1) starting nitrate, (Stage 2) low-nitrate
temperature residual, (Stage 3) cruise-specific bottle bias, (Stage 4) assessment vs bottles.

- **Inputs:** `L2/SUNA/<cast>_SUNA_1s.csv`, raw `L1/SUNA/*.csv`, `metadata/calibration/SNA2503C.cal`,
  and the **bottle nitrate table** (see `README_bottle_columns.md`)
- **Output:** `L3/<cast>_nitrate_QC.csv` + assessment tables in `_audit/`

### step04 — Profile plots + lab-friendly CSVs → L2/CTD/profile_plots
Vertical profiles from the 1 m downcast CTD data and, where a cast has SUNA, the merged nitrate.
Per-cast single-variable figures and per-variable multi-cast overlays, plus a lab-friendly `<cast>_readable.csv`.

### step05 — Profile diagnostics → L2/CTD/profile_plots + _audit
Thermocline / halocline / pycnocline / oxycline, DCM, MLD (temperature and density criteria), and the
nitracline (onset + max gradient), as annotated multi-panel figures with a per-cast `.txt` sidecar.

---

## 4b. Does it produce a CSV? Yes

| Notebook/script | Output file | Folder | Produced when |
|-----------------|-------------|--------|---------------|
| step01 | `<cast>_1m_down.cnv`, `<cast>_1s.cnv` (+`.ros`/`.btl` in bottle mode) | `L2/CTD/...` | `RUN_SBE_COMMANDS = True`, PSAs present |
| step02 | `<cast>_SUNA_1s.csv` | `L2/SUNA/` | step01 produced the 1 s CNV |
| step03 | `<cast>_nitrate_QC.csv` | `L3/` | step02 produced the merged file |
| step04 | `<cast>__<var>.png`, `<cast>_readable.csv` | `L2/CTD/profile_plots/` | step01 (+ step02 for nitrate) done |
| step05 | `<cast>__annotated.png`, `<cast>_diagnostics.txt`, `profile_diagnostics.csv` | `profile_plots/` + `_audit/` | step01 (+ step02) done |
| `cnv_to_csv.py` (helper) | `<name>.csv` per `.cnv` | chosen output folder | run on demand |

`step03`'s `_nitrate_QC.csv` columns: `utc_time`, `temp_c`, `salinity`, `no3_onboard`, `no3_stage1`,
`no3_stage2`, `no3_qc` (final). Without a bottle table, Stages 2–3 are skipped and `no3_qc` carries no
bottle correction — fill the bottle table before relying on `no3_qc`.

---

## 5. Decisions (resolved this cycle)

These were open in v1 and are now settled for the Gulf of Guinea shelf work:

1. **Filter** — ON, before Align (P 1.0 s; T + C 0.5 s). Measured effect on the P45_06 casts:
   conductivity dX/dt std reduced 33–54 % (mean 46 %), pressure 12–34 % (mean 24 %), temperature ~0 %.
2. **Cell Thermal Mass** — ON, after Align (`APPLY_CELLTM = True`, `04_celltm.psa`). Coefficients
   `alpha 0.03 / (1/beta) 7.0` are the SBE defaults; **confirm with Drew** (some groups use 0.04 / 8.0).
3. **Section before Loop Edit** — added, because our staged equilibration dip made Loop Edit over-flag
   (34–81 % on the full cast). Removing the leading soak + dip first drops that to a few percent.
4. **Align advances** — temperature +0.5 s (fixed); conductivity **0 s** (the swept optimum was
   marginal and inconsistent across casts, so no residual advance); oxygen +5 s (mid-range, pending
   bottle-O2 calibration). SBE 43 **Tau correction OFF**, hysteresis ON.
5. **Loop Edit minimum velocity** — lowered 0.25 → **0.10 m/s** for the slow shelf descents.
6. **Wild Edit / despiking** — left OUT of the chain; the T + C low-pass Filter handles the spiking that
   matters for salinity. Revisit only if a specific cast needs it.

---

## 5b. Running in PyCharm Community Edition

1. **Open the project** at the top level: `ctd_pipeline` (e.g. `C:\Projects\ctd_pipeline`).
2. **Create the virtual environment:** Settings → Project → Python Interpreter → Add → Virtualenv →
   New, location `ctd_pipeline/.venv`, base Python 3.10/3.11 (3.12+ also fine).
3. **Install packages** in the PyCharm terminal (prompt shows `(.venv)`): `pip install -r requirements.txt`.
4. **Run step01:** open `step01_run_sbe_processing_v2.py`, edit the settings block
   (`CRUISE_ID`, `SBE_BIN_DIR`, `APPLY_BOTTLES`, `RUN_SBE_COMMANDS = False` for the dry run), press Run.
5. **Run step02 / step03 / step04 / step05:** in the terminal run `jupyter notebook`, open the `.ipynb`
   in the browser, and run the cells there.

---

## 6. Before the first run — checklist

- [ ] **`01_datcnv.psa`** present in `psa/` (variable selection, full cast down+up, no binning,
      SBE 43 O2 Tau OFF / hysteresis ON). The bottle variant is auto-generated when needed.
- [ ] The other chain PSAs present: `05_filter`, `02_alignctd`, `04_celltm`,
      `06a_loopedit_percentmean` (+`06b`), `07_derive`, `08b_binavg_1m_down`, `08d_binavg_1s_full`.
- [ ] Every PSA has **blank** Input/Output directory and Name append (so the runner's CLI flags control
      output). Confirm via `_audit/01_psa_diagnostic.csv` — `warnings` column empty.
- [ ] CTD config (`19-8460.xmlcon`) in the cruise `Instrument_Calibration_Files/`.
- [ ] SUNA cal (`SNA2503C.cal`) in `metadata/calibration/`.
- [ ] For bottle mode: each cast's `.afm` sits **beside its `.hex` in `L0/CTD/`** (same basename), and
      `10_bottlesum.psa` present if you want `.btl` summaries.
- [ ] Same `CRUISE_ID` across step01/02/03; `APPLY_BOTTLES` matches the cruise.
- [ ] Fill the bottle table per `README_bottle_columns.md` (for the calibrated nitrate product).
- [ ] step01: dry-run first, inspect the audit CSVs, then `RUN_SBE_COMMANDS = True`.

---

## 7. Fixes and lessons from the real-data runs

Recorded so they are not re-introduced:

- **Section before Loop Edit** is essential given the staged equilibration dip; sectioning must happen
  *before* Loop Edit, not after (order is not interchangeable). The Section step preserves the SBE
  `.cnv` **CRLF** line endings byte-for-byte (writing LF-only makes SBE report
  "Header line length exceeds buffer length").
- **Per-stage L2 folders.** Each module writes to its own folder (`Filter/`, `Align_CTD/`, `Cell_tm/`,
  `Section/`, `Loop_Edit/`, `Derived/`). The derived hub was renamed `Derived_Parameter/` → `Derived/`.
- **Bottle mode `CreateFile`.** "Create both data and bottle file" is `CreateFile=2` in the PSA;
  `CreateFile=1` is bottle-only (writes `.ros`, no `.cnv`). The auto-generated bottle PSA uses `2` and
  self-heals a stale one.
- **No-afm casts in bottle mode** fall back to the profile PSA automatically (bottle-only DatCnv would
  otherwise halt with no `.cnv`).
- **Work-folder cleanup** clears read-only (raw `.hex` are often archived read-only), retries transient
  locks, and falls back to a fresh timestamped work folder. A **per-module timeout** prevents a stuck
  SBE dialog from hanging the run and orphaning a process.
- **PSA output paths were hardcoded** in the original set; blanked so the runner's CLI flags control
  output. Re-blank (or re-save cleanly) if you ever re-export a PSA from the SBE GUI.
- **SUNA files are raw SATSLF frames**, paired to casts by time overlap via `suna_cast_map.csv`.
- **19plus V2 specifics.** Strain-gauge pressure (not Digiquartz); temperature `tv290C`, conductivity
  `c0S/m`, salinity `sal00`/`sal78` downstream.

## 8. Known limitation — Stage 1 TSP

`step03`'s `STAGE1_MODE` defaults to `"onboard"` (SUNA firmware T/S/P-corrected nitrate as the starting
point; bottle-based Stages 2–3 do the main accuracy correction). A `"tsp"` mode that re-derives from
raw spectra is included but **work-in-progress** (needs the `OPTICAL_WAVELENGTH_OFFSET` baseline,
Sakamoto 2009); keep `"onboard"` until validated.

---

## 9. Helper tools (standalone utilities)

Run directly (`python <script>.py`); most have a settings block at the top.

| Script | What it does | When to run |
|--------|--------------|-------------|
| `cnv_to_csv.py` | Pick an input folder of `.cnv` and an output folder; converts every `.cnv` to `.csv` (headers from the `# name` lines; values copied through exactly). Pure standard library. | Any time someone needs the `.cnv` as plain CSV. |
| `verify_reproducibility.py` | Compares two step01 output trees by `.cnv` **data columns** (ignoring the run-timestamp header lines); reports per-variable max abs diff and flag counts. | To confirm a rerun reproduces the same numbers. |
| `build_bottle_table.py` | Pairs discrete nutrient bottle samples to CTD/SUNA casts by time, with confidence/reason; writes the audit trail and the step03-ready `<CRUISE>_bottle_nitrate.csv`. | After nutrient lab data arrives, before step03. |
| `step03_run_report.py` | One readable report from step03's CSVs: bottle coverage, RMSE progression, per-cast medians, fit coefficients, verdict. | After step03. |
| `combine_readable_profiles.py` | Stacks step04's per-cast readable CSVs into one cruise-wide table. | After step04. |
| `fix_psa_output_paths.py` | Blanks hardcoded OutputDir/OutputFile/NameAppend in the PSAs. | Once, or after re-saving any PSA from the SBE GUI. |

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

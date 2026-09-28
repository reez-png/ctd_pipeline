# Walkthrough: rerun the .hex files on the new workflow, then check reproducibility

**Cruise:** P45_06 - **Instrument:** SBE 19plus V2 (`19-8460.xmlcon`) - **Software:** SBE Data Processing v7.26.7

This covers the two files:

- `step01_run_sbe_processing_v2.py` - updated to the agreed chain.
- `verify_reproducibility.py` - checks a rerun reproduces the same numbers.

One honest note first: I could not run any of this myself. This session runs in a Linux
container with no Sea-Bird binaries and no access to your `.hex` files, so the actual
processing happens on your Windows machine. What I did do is rewrite the chain, and validate
all the new Python logic here (the Section step and the reproducibility comparator) on
synthetic `.cnv` files. So the Python parts are tested; the SBE binary calls are what you run.

---

## 1. What changed in step01

Old canonical chain (what the script had):

```
DatCnv -> AlignCTD -> [WildEdit off] -> [CellTM off] -> [Filter off] -> LoopEdit -> Derive
```

New chain (Amala-approved, what the script does now):

```
DatCnv -> Filter -> AlignCTD -> CellTM -> Section (python) -> LoopEdit -> Derive
                                                                              |
                                              Bin Average: 1 m downcast  +  1 s full (SUNA)
```

Concretely:

- **Filter is ON and moved BEFORE Align** (pressure 1.0 s Filter A; temperature and
  conductivity 0.5 s Filter B). PSA: `05_filter.psa`. Output suffix `_filt.cnv`.
- **Cell Thermal Mass is ON and runs AFTER Align.** PSA: `04_celltm.psa`. Suffix `_ctm.cnv`.
- **Section is a new NON-SBE Python step** between Cell Thermal Mass and Loop Edit. It removes
  the leading surface soak and equilibration dip and keeps downcast + upcast, so Loop Edit's
  running-max pressure rule is not poisoned by the dip. Suffix `_sec.cnv`. This is the
  `section_cnv()` function inside step01; it reproduces the motion-based descent onset we
  validated in step00b/step00c and preserves the data rows byte-for-byte (only a contiguous
  tail is kept; the header `# nvalues` is updated).
- **Loop Edit now points at `06a_loopedit_percentmean.psa`** (percent of mean speed, minV
  0.10 m/s), with `06_loopedit.psa` kept as a fallback name.
- **WildEdit is left out** of the canonical chain.
- Bin Average unchanged: `08b_binavg_1m_down` (science, downcast only) and `08d_binavg_1s_full`
  (for the SUNA merge).

Gating flags at the top of the script (section 5 of the settings):

```python
APPLY_WILDEDIT = False
APPLY_CELLTM   = True     # on
APPLY_FILTER   = True     # on, before Align
APPLY_SECTION  = True     # python section before Loop Edit
```

---

## 2. Before you run: PSA files

The script copies each module's PSA from your `psa/` folder into the no-space work folder.
Make sure these exist in `C:\Projects\ctd_pipeline\psa\`:

| Module | PSA file | Notes |
|--------|----------|-------|
| DatCnv | `01_datcnv.psa` | Tau correction OFF, hysteresis ON (the version I sent) |
| Filter | `05_filter.psa` | P 1.0 s, T + C 0.5 s (the version I sent) |
| Align CTD | `02_alignctd.psa` | T +0.5 s, C 0 s, O2 +5 s |
| Cell Thermal Mass | `04_celltm.psa` | your existing CellTM alpha / 1-beta |
| Loop Edit | `06a_loopedit_percentmean.psa` | percent of mean, minV 0.10 m/s |
| Derive | `07_derive.psa` | salinity, density, etc. |
| Bin Average (down) | `08b_binavg_1m_down.psa` | 1 m, downcast only |
| Bin Average (1 s) | `08d_binavg_1s_full.psa` | 1 s, full cast |

Drop in the three I edited (`01_datcnv.psa`, `05_filter.psa`, and `06a_loopedit_percentmean.psa`)
so they match this run. The Section step needs no PSA.

---

## 3. Dry run first (no Sea-Bird binary is launched)

In the settings block near the top:

```python
TEST_SINGLE_CAST_ONLY = True
TEST_CAST_ID          = "P45_06_CTD_01"
RUN_SBE_COMMANDS      = False     # dry run
```

Run the script (PyCharm: right-click -> Run). It writes the command plan and audit tables to
`cruises\P45_06\_audit\` without touching any binary. Check:

- `00_preflight_module_check.csv` - every enabled module shows `status = OK`, and
  `05_section` shows `OK (python step)` with exe `(python:python_section)`.
- `01_psa_diagnostic.csv` - **most important.** Any PSA with a warning ("OutputDir is
  hardcoded", "AppendOutputFile=1") will silently ignore the command-line output path and is
  the usual cause of "return code 0 but no output". If a PSA is flagged, open it once in the
  SBE Data Processing GUI, set the output to use the input folder + base name with no append,
  and save.
- `02_sbe_processing_command_log.csv` - the exact command lines, in order. Confirm the chain
  reads DatCnv -> 02_filter -> 03_alignctd -> 04_celltm -> 05_section -> 06_loopedit ->
  07_derive, then the two bin branches.

---

## 4. Real run on CTD_01

Flip one switch and run again:

```python
RUN_SBE_COMMANDS = True
```

Then check in `_audit\`:

- `03_l1_l2_deliverable_index.csv` - `exists = True` for the L1 `.cnv`, the `_1m_down.cnv`
  and the `_1s.cnv`.
- `per_cast_sbe_outputs\P45_06_CTD_01\logs\05_section.log` - the Section step's report:
  rows in, rows kept, rows removed (soak + dip), onset scan and pressure. Sanity check that
  the onset pressure is a couple of db (just below the surface after the dip), not 11 db.
- `per_cast_sbe_outputs\P45_06_CTD_01\logs\06_loopedit.log` - Loop Edit should now flag a
  small fraction, not the 34 to 81% we saw before sectioning.

The stage files land in `cruises\P45_06\L2\CTD\...`:
`_filt.cnv` and `_al.cnv` and `_ctm.cnv` and `_sec.cnv` in `Align_CTD\`, `_loop.cnv` in
`Loop_Edit\`, `_der.cnv` in `Derived_Parameter\`, and the binned products in `CTD_1m_down\`
and `CTD_1s\`.

---

## 5. Check it reproduces (the "same things" test)

Sea-Bird writes a run timestamp into every `.cnv` header, so two runs are never byte-identical.
The test that matters is whether the DATA columns come out the same. `verify_reproducibility.py`
does exactly that: it matches files by name, and for each pair reports row count, per-variable
maximum absolute difference, and bad-flag counts.

1. After your first real run, copy the L2 tree to a baseline:
   `cruises\P45_06\L2`  ->  `cruises\P45_06\L2_baseline` (Explorer copy is fine).
2. Run step01 again (same `.hex`, same PSAs). It overwrites `L2`.
3. Edit the two paths at the top of `verify_reproducibility.py`:

   ```python
   BASELINE_DIR = Path(r"C:\Projects\ctd_pipeline\cruises\P45_06\L2_baseline")
   RERUN_DIR    = Path(r"C:\Projects\ctd_pipeline\cruises\P45_06\L2")
   GLOB_PATTERN = "*.cnv"     # or "*_1m_down.cnv" to check only the science product
   ```

   Run it. A deterministic pipeline gives `max_abs_diff = 0.0` for every file and prints
   `REPRODUCIBLE`. Any row with `identical = False` names the file and the worst variable, so
   you can see exactly what moved (usually a PSA or setting that changed between runs).

The report is saved to `_audit\reproducibility_report.csv`.

---

## 6. All casts

Once CTD_01 looks right:

```python
TEST_SINGLE_CAST_ONLY = False
RUN_SBE_COMMANDS      = True
```

Run once more to process every `.hex` under `L0\CTD\`. `STOP_ON_ERROR = True` halts at the
first failing cast so a bad file cannot quietly corrupt the batch; the log names the cast and
module. CTD_04 is the one we flagged earlier (anomalous conductivity) - expect to inspect it.

---

## 7. Gotchas (already handled, worth knowing)

- **No-space work folder.** SBE command-line modules choke on spaces in paths, so the script
  stages inputs into `C:\Projects\ctd_pipeline\_sbe_work` and runs from there. Do not move it
  under a path with spaces.
- **OneDrive placeholders.** A zero-byte `.hex`, `.xmlcon` or `.psa` (Files-On-Demand not
  downloaded) raises a clear error. Right-click the folder -> "Always keep on this device".
- **The Section step is deterministic** but depends on the pressure column and sample rate in
  the header; it reads both automatically and falls back to the surface crossing if a cast has
  a near-bottom stall (the CTD_05 failure mode from the QC notebook), so it cannot collapse to
  a 2 m sliver.

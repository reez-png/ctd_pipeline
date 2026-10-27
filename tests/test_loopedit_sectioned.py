"""
test_loopedit_sectioned.py

Confirms the sectioning fix for Loop Edit on shallow shelf casts with an equilibration dip.

For each cast it runs Loop Edit THREE ways and reports flagged-scan %:
  (a) full raw cast       - leading surface soak + equilibration dip + descent + upcast, nothing cut
  (b) sectioned cast      - leading soak and equilibration dip removed, but the whole real cast
                            (downcast AND upcast) kept. This is the actual pipeline input to
                            Loop Edit -> Bin Average.
  (c) sectioned downcast  - the descent only (soak, dip AND upcast removed). The clean-descent
                            diagnostic, and effectively what Bin Average keeps as the downcast.

Note downcast-mode Loop Edit flags the upcast by design, so (b) stays high while (c) is the low,
meaningful number that shows loop-edit settings are fine once the dip is out.

Loop edit SETTINGS are unchanged (minV 0.10 m/s, percent-mean window 3 s / fixed-min). Only the INPUT
changes: feed it the sectioned cast, not the raw cast. Sectioning FIRST matters because on the raw
cast the equilibration dip drives Loop Edit's running maximum pressure to ~11 db before the real
descent, so the top of the descent is wrongly flagged as a reversal. Removing the leading soak+dip
before Loop Edit fixes that. (The sectioned cast still contains the upcast, which downcast-mode
Loop Edit flags by design; that is resolved later at Bin Average, not here.)

Section start is found by MOTION, not a fixed depth line: a fixed soak-depth cutoff fails because
the post-dip equilibration plateau can sit above or below any chosen depth. We smooth pressure and,
walking back from the deepest scan, find the descent onset (last scan still genuinely descending),
then keep everything from there to the end of the record (descent + upcast).

Because SBE Loop Edit reads a .cnv file, we write the sectioned cast to a temporary .cnv, preserving
the original header, then run Loop Edit on it.

This is a LOCAL diagnostic script, not a unit test: it runs the Sea-Bird software on real cruise
data. The executable part is under `if __name__ == "__main__":` so pytest can import the module
without launching Sea-Bird (so CI, and any machine without SBE, are unaffected).

USAGE: python test_loopedit_sectioned.py
"""
from pathlib import Path
import re, subprocess
import numpy as np, pandas as pd

CTD_ROOT = Path(r"C:\Projects\ctd_pipeline")
CRUISE   = "P45_06"
CASTS    = ["P45_06_CTD_03", "P45_06_CTD_06"]
SBE_BIN  = Path(r"C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32")
PSA_DIR  = CTD_ROOT / "psa"
ALIGN_DIR= CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "Align_CTD"
FILTER_PSA          = PSA_DIR / "05_filter.psa"
LOOPEDIT_PCTMEAN    = PSA_DIR / "06a_loopedit_percentmean.psa"
LOOPEDIT_FIXEDMIN   = PSA_DIR / "06b_loopedit_fixedmin.psa"
WORK     = CTD_ROOT / "_sbe_work" / "qc"

# Motion-based descent-onset detection (must match step00b_qc_filter_loopedit.ipynb).
DESCENT_SMOOTH_SCANS = 24     # ~1 s median smoothing at 24 Hz, to beat pressure noise
DESCENT_MIN_DROP_DB  = 0.10   # min pressure gain over the window for a scan to count as descending


def read_cnv_full(path):
    """Return (header_lines, cols, data_df, rows_raw) so a sectioned .cnv can be re-written."""
    txt = Path(path).read_text(encoding="latin-1", errors="replace").splitlines()
    names, ds = {}, None
    nr = re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*[:=]", re.I)
    for i, l in enumerate(txt):
        s = l.strip()
        if s.startswith("#") or s.startswith("*"):
            m = nr.search(l)
            if m: names[int(m.group(1))] = m.group(2).strip()
            if s.upper().startswith("*END*"): ds = i + 1; break
    header = txt[:ds]
    cols = [names[k] for k in sorted(names)]
    rows_raw = [l for l in txt[ds:] if l.strip()]
    data = pd.DataFrame([r.split()[:len(cols)] for r in rows_raw if len(r.split()) >= len(cols)],
                        columns=cols).apply(pd.to_numeric, errors="coerce")
    return header, cols, data, rows_raw


def pcol_of(cols):
    for c in ("prdM", "prDM", "prM"):
        if c in cols: return c
    return cols[0]


def descent_start(pv, win=DESCENT_SMOOTH_SCANS, min_drop=DESCENT_MIN_DROP_DB):
    """Positional index of the true descent onset for pressure array `pv` (soak + equilibration dip
    excluded). Walk back from the deepest scan keeping every scan still genuinely descending."""
    pv = np.asarray(pv, float)
    if pv.size == 0 or np.isnan(pv).all():
        return 0
    imax = int(np.nanargmax(pv))
    if imax < win + 1:
        return 0
    seg = pv[:imax + 1]                                   # downcast side only (never peek upcast)
    s = pd.Series(seg).rolling(win, center=True, min_periods=1).median().to_numpy()
    fwd = np.full(len(s), np.nan)
    fwd[:len(s) - win] = s[win:] - s[:len(s) - win]       # pressure gain over the next `win` scans
    descending = np.where(np.isnan(fwd), True, fwd > min_drop)
    i = imax
    while i > 0 and descending[i - 1]:
        i -= 1
    return i


def cast_bounds(data, pcol):
    """(start, end) of the real cast: start = descent onset (soak + dip removed), end = last scan of
    the record (upcast retained)."""
    pv = pd.to_numeric(data[pcol], errors="coerce").to_numpy(float)
    return descent_start(pv), len(pv) - 1


def write_sectioned_cnv(src_cnv, out_cnv, start, end):
    """Write a .cnv with the same header but only rows [start:end] (the sectioned cast:
    leading soak + equilibration dip removed, downcast + upcast kept)."""
    header, cols, data, rows_raw = read_cnv_full(src_cnv)
    keep = rows_raw[start:end + 1]
    out_cnv.write_text("\n".join(header + keep) + "\n", encoding="latin-1")
    return out_cnv


def run(exe, psa, in_cnv, stem):
    subprocess.run([str(SBE_BIN / exe), f"/p{psa}", f"/i{in_cnv}", f"/o{WORK}", f"/f{stem}", "/s"],
                   capture_output=True, text=True)
    out = WORK / f"{stem}.cnv"
    return out if out.exists() else None


def flag_pct(cnv):
    _, cols, data, _ = read_cnv_full(cnv)
    if "flag" not in data.columns or len(data) == 0: return None
    f = pd.to_numeric(data["flag"], errors="coerce").fillna(0) != 0
    return 100 * f.mean(), int(f.sum()), len(data)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    print(f"{'cast':<16} {'method':<14} {'full raw cast':>16} {'sectioned':>16} {'sect. downcast':>16}")
    print("-" * 82)
    for cast in CASTS:
        aligned = ALIGN_DIR / f"{cast}_al.cnv"
        if not aligned.exists():
            print(f"{cast}: no aligned cnv"); continue
        filt = run("FilterW.exe", FILTER_PSA, aligned, f"{cast}_filt")
        src = filt if filt else aligned

        # section the FILTERED cnv two ways, both with the leading soak + equilibration dip removed:
        #   sect  = descent onset -> end of record (downcast + upcast)
        #   dcast = descent onset -> bottom          (descent only)
        _, cols, data, _ = read_cnv_full(src)
        pcol = pcol_of(cols)
        pv = pd.to_numeric(data[pcol], errors="coerce").to_numpy(float)
        start = descent_start(pv); imax = int(np.nanargmax(pv)); end = len(pv) - 1
        sect  = write_sectioned_cnv(src, WORK / f"{cast}_sect.cnv",  start, end)
        dcast = write_sectioned_cnv(src, WORK / f"{cast}_dcast.cnv", start, imax)

        for label, psa in [("percent-mean", LOOPEDIT_PCTMEAN), ("fixed-min", LOOPEDIT_FIXEDMIN)]:
            outs = {k: run("LoopEditW.exe", psa, v, f"{cast}_le_{k}")
                    for k, v in (("full", src), ("sect", sect), ("dcast", dcast))}
            def fmt(o):
                r = flag_pct(o) if o else None
                return f"{r[0]:.1f}% ({r[1]}/{r[2]})" if r else "n/a"
            print(f"{cast:<16} {label:<14} {fmt(outs['full']):>16} {fmt(outs['sect']):>16} {fmt(outs['dcast']):>16}")

    print("\n'sectioned'      = Loop Edit on the cast with only the leading soak + equilibration dip")
    print("                   removed (downcast AND upcast kept) - the real pipeline input.")
    print("'sect. downcast' = same, but the upcast is also removed (descent only) - the clean-descent")
    print("                   diagnostic, and effectively what Bin Average keeps as the downcast.")
    print("Sectioning FIRST stops the dip from poisoning Loop Edit's running max-pressure and wrongly")
    print("flagging the top of the descent. The upcast is flagged by design, so 'sectioned' stays high")
    print("while 'sect. downcast' is the low, meaningful number.")


if __name__ == "__main__":
    main()

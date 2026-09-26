"""
test_loopedit_sectioned.py

Confirms the fix for Loop Edit over-flagging on shallow shelf casts with an equilibration dip.

For each cast it runs Loop Edit TWO ways and reports flagged-scan %:
  (a) full raw cast   - includes soak + equilibration dip + upcast (over-flags ~80%)
  (b) true downcast   - sliced to the real descent only (dip trimmed), which is what should
                        be processed. Expected to drop to a normal ~2-10%.

Loop edit SETTINGS are unchanged (0.25 m/s, percent-mean, window 3 s). Only the INPUT changes:
feed it the sectioned true-downcast, not the raw cast. This proves the problem was cast geometry
(un-trimmed equilibration dip), not loop-edit configuration.

Because SBE Loop Edit reads a .cnv file, we write the sectioned downcast to a temporary .cnv,
preserving the original header, then run Loop Edit on it.

USAGE: python test_loopedit_sectioned.py
"""
from pathlib import Path
import re, subprocess
import numpy as np, pandas as pd

CTD_ROOT = Path(r"C:\Projects\ctd_pipeline")
CRUISE   = "P45_06"
CASTS    = ["P45_06_CTD_03", "P45_06_CTD_06"]
NEAR_SURFACE_DB = 2.0
SBE_BIN  = Path(r"C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32")
PSA_DIR  = CTD_ROOT / "psa"
ALIGN_DIR= CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "Align_CTD"
FILTER_PSA          = PSA_DIR / "05_filter.psa"
LOOPEDIT_PCTMEAN    = PSA_DIR / "06a_loopedit_percentmean.psa"
LOOPEDIT_FIXEDMIN   = PSA_DIR / "06b_loopedit_fixedmin.psa"
WORK     = CTD_ROOT / "_sbe_work" / "qc"; WORK.mkdir(parents=True, exist_ok=True)


def read_cnv_full(path):
    """Return (header_lines, cols, data_df) so we can re-write a sectioned .cnv."""
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


def true_descent_slice(data, pcol):
    """Indices [start:bottom] of the true downcast (dip trimmed)."""
    p = pd.to_numeric(data[pcol], errors="coerce").to_numpy(float)
    imax = int(np.nanargmax(p))
    i = imax
    while i > 0 and p[i] > NEAR_SURFACE_DB:
        i -= 1
    return i, imax


def write_sectioned_cnv(src_cnv, out_cnv, start, bottom):
    """Write a .cnv with the same header but only rows [start:bottom]."""
    header, cols, data, rows_raw = read_cnv_full(src_cnv)
    keep = rows_raw[start:bottom + 1]
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


print(f"{'cast':<16} {'method':<14} {'full cast':>14} {'true downcast':>16}")
print("-" * 64)
for cast in CASTS:
    aligned = ALIGN_DIR / f"{cast}_al.cnv"
    if not aligned.exists():
        print(f"{cast}: no aligned cnv"); continue
    filt = run("FilterW.exe", FILTER_PSA, aligned, f"{cast}_filt")
    src = filt if filt else aligned

    # section the FILTERED cnv to the true downcast
    _, cols, data, _ = read_cnv_full(src)
    pcol = pcol_of(cols)
    start, bottom = true_descent_slice(data, pcol)
    sect = write_sectioned_cnv(src, WORK / f"{cast}_down.cnv", start, bottom)

    for label, psa in [("percent-mean", LOOPEDIT_PCTMEAN), ("fixed-min", LOOPEDIT_FIXEDMIN)]:
        full_out = run("LoopEditW.exe", psa, src,  f"{cast}_le_full")
        down_out = run("LoopEditW.exe", psa, sect, f"{cast}_le_down")
        fp = flag_pct(full_out) if full_out else None
        dp = flag_pct(down_out) if down_out else None
        fs = f"{fp[0]:.1f}% ({fp[1]}/{fp[2]})" if fp else "n/a"
        dsx = f"{dp[0]:.1f}% ({dp[1]}/{dp[2]})" if dp else "n/a"
        print(f"{cast:<16} {label:<14} {fs:>14} {dsx:>16}")

print("\nThe 'true downcast' column is loop edit run on the sectioned descent (dip trimmed).")
print("If it drops to ~2-10%, it confirms: loop-edit settings were fine all along; the fix is")
print("to section the true downcast (remove the equilibration dip) BEFORE loop edit.")

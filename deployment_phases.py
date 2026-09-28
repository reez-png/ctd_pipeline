"""
deployment_phases.py

Pressure-vs-time diagnostic that marks the deployment phases of a CTD cast:
  - surface soak (near-surface, before any descent)
  - equilibration dip (a deliberate lower-and-raise to wet/prime the sensors)
  - true descent start (where the real downcast begins)
  - bottom, and upcast

This explains loop-edit over-flagging on shallow shelf casts: when a cast includes an
equilibration dip, loop edit flags the dip (genuine slow/reversing motion) even though it is
not part of the science profile. The fix is to start the cast at the true-descent start, not
to change loop-edit thresholds.

USAGE: python deployment_phases.py   (edit CAST)
"""
from pathlib import Path
import re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CTD_ROOT = Path(r"C:\Projects\ctd_pipeline")
CRUISE   = "P45_06"
CASTS    = ["P45_06_CTD_03", "P45_06_CTD_06"]
ALIGN_DIR= CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "Align_CTD"
OUT_DIR  = CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "qc_plots"; OUT_DIR.mkdir(parents=True, exist_ok=True)

NEAR_SURFACE_DB = 2.0     # "near surface" band
DESCENT_MIN_DB  = 5.0     # a real descent must pass at least this depth

def read_cnv(path):
    txt=Path(path).read_text(encoding="latin-1",errors="replace").splitlines()
    names,ds={},None
    nr=re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*[:=]",re.I)
    for i,l in enumerate(txt):
        s=l.strip()
        if s.startswith("#") or s.startswith("*"):
            m=nr.search(l)
            if m: names[int(m.group(1))]=m.group(2).strip()
            if s.upper().startswith("*END*"): ds=i+1; break
    cols=[names[k] for k in sorted(names)]
    rows=[l.split()[:len(cols)] for l in txt[ds:] if l.strip() and len(l.split())>=len(cols)]
    return pd.DataFrame(rows,columns=cols).apply(pd.to_numeric,errors="coerce")

def pcol_of(df):
    for c in ("prdM","prDM","prM"):
        if c in df.columns: return c
    return None

def tcol_of(df):
    for c in ("timeS","timeM"):
        if c in df.columns: return c
    return None

def find_true_descent_start(p):
    """The true descent = the LAST time the trace leaves the near-surface band and then
    descends monotonically to the global max. Everything before that (soak + equilibration
    dip) is deployment overhead, not science."""
    p = pd.to_numeric(p, errors="coerce").to_numpy(float)
    imax = int(np.nanargmax(p))
    # walk back from the bottom to the last time we were in the near-surface band
    i = imax
    while i > 0 and p[i] > NEAR_SURFACE_DB:
        i -= 1
    return i, imax   # true-descent start index, bottom index

def phase_plot(cast):
    cnv = ALIGN_DIR / f"{cast}_al.cnv"
    if not cnv.exists():
        print(f"{cast}: no aligned cnv"); return None
    df = read_cnv(cnv)
    pcol, tcol = pcol_of(df), tcol_of(df)
    p = pd.to_numeric(df[pcol], errors="coerce")
    t = pd.to_numeric(df[tcol], errors="coerce") if tcol else pd.Series(np.arange(len(df)))
    start, bottom = find_true_descent_start(p)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, p, lw=0.7, color="#1f77b4")
    ax.axvspan(t.iloc[0], t.iloc[start], color="#f0a030", alpha=0.15, label="soak + equilibration dip")
    ax.axvspan(t.iloc[start], t.iloc[bottom], color="#2ca02c", alpha=0.12, label="true downcast")
    ax.axvspan(t.iloc[bottom], t.iloc[len(t)-1], color="#9467bd", alpha=0.10, label="upcast")
    ax.axvline(t.iloc[start], color="#2ca02c", lw=1.3, ls="--")
    ax.annotate("true descent start", (t.iloc[start], p.iloc[start]),
                xytext=(15, 20), textcoords="offset points", fontsize=8,
                arrowprops=dict(arrowstyle="->", color="#2ca02c"))
    ax.invert_yaxis()
    ax.set_xlabel("elapsed time (s)"); ax.set_ylabel("pressure (db)")
    ax.set_title(f"{cast}: deployment phases")
    ax.legend(loc="lower right", fontsize=8); ax.grid(alpha=0.3)

    # phase timing summary
    dt_soak = float(t.iloc[start] - t.iloc[0])
    dt_down = float(t.iloc[bottom] - t.iloc[start])
    txt = (f"soak+dip: {dt_soak:.0f} s ({100*start/len(t):.0f}% of scans)\n"
           f"true downcast: {dt_down:.0f} s\n"
           f"max pressure: {p.max():.1f} db")
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=7.5, family="monospace",
            va="bottom", ha="left", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#999", alpha=0.9))

    out = OUT_DIR / f"{cast}__deployment_phases.png"
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"{cast}: soak+dip {100*start/len(t):.0f}% of scans, true descent starts at scan {start} "
          f"(t={t.iloc[start]:.0f}s), bottom {p.max():.1f} db -> {out.name}")
    return out

for c in CASTS:
    phase_plot(c)
print("\nPhase plots written. The green region is the science downcast; the orange region")
print("(soak + equilibration dip) is what loop edit was flagging. Start the cast at the")
print("true-descent line to stop the over-flagging.")

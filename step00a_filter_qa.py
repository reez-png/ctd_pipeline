"""
step00a_filter_qa.py

QA/QC for the Filter module: quantify what the low-pass Filter actually does to each cast, before
Loop Edit. Companion to step00b (which does the Filter + Loop Edit before/after panels); this one
isolates the Filter and reports the numbers you would quote in a methods section.

For each cast it:
  1. runs FilterW.exe on the aligned cnv (same PSA the pipeline uses);
  2. computes, for pressure (and temperature/conductivity if present), the scan-to-scan change
     dX/dt = gradient(X), and its standard deviation BEFORE vs AFTER filtering;
  3. reports the reduction:  reduction% = 100 * (1 - std_filtered / std_before).
     For a strain-gauge 19plus V2 the pressure reduction is expected to be small (little dithering
     to remove), which is the point: it shows Filter is doing little here, unlike a 9plus/Digiquartz.

It also reads the PSA back and prints exactly which variables the Filter is set to touch and with
which time constant, so there is no ambiguity about what was applied.

Outputs: a per-cast/-variable table to stdout, a summary CSV, and a before/after dP/dt plot per cast.

USAGE: python step00a_filter_qa.py
"""
from pathlib import Path
import re, subprocess
import xml.etree.ElementTree as ET
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---------------------------------------------------------------- paths / config
CTD_ROOT  = Path(r"C:\Projects\ctd_pipeline")
CRUISE    = "P45_06"
SBE_BIN   = Path(r"C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32")
PSA_DIR   = CTD_ROOT / "psa"
ALIGN_DIR = CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "Align_CTD"
FILTER_PSA = PSA_DIR / "05_filter.psa"
WORK      = CTD_ROOT / "_sbe_work" / "qc"; WORK.mkdir(parents=True, exist_ok=True)
QC_PLOT_DIR = CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "qc_plots"; QC_PLOT_DIR.mkdir(parents=True, exist_ok=True)

# all aligned casts in the folder; to restrict, replace with an explicit list e.g. ["P45_06_CTD_03"]
CASTS = sorted(p.name[:-len("_al.cnv")] for p in ALIGN_DIR.glob("*_al.cnv"))

# variables to test: (cnv short-name candidates, label)
VARS = [(("prdM", "prDM", "prM"), "pressure"),
        (("tv290C", "t090C", "tv290"), "temperature"),
        (("c0S/m", "cond0S/m", "c0mS/cm", "c_S/m"), "conductivity")]

rcParams.update({"figure.dpi": 110, "savefig.dpi": 300, "font.size": 9, "axes.titlesize": 10,
                 "axes.labelsize": 9, "axes.linewidth": 0.8, "xtick.direction": "in",
                 "ytick.direction": "in", "axes.grid": False, "font.family": "DejaVu Sans"})
C_BEFORE, C_FILTER = "#E8703A", "#2E6FB0"


def _despine(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def read_cnv(path):
    """Parse a .cnv keeping all scans."""
    txt = Path(path).read_text(encoding="latin-1", errors="replace").splitlines()
    names, ds = {}, None
    nr = re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*[:=]", re.I)
    for i, l in enumerate(txt):
        s = l.strip()
        if s.startswith("#") or s.startswith("*"):
            m = nr.search(l)
            if m: names[int(m.group(1))] = m.group(2).strip()
            if s.upper().startswith("*END*"): ds = i + 1; break
    if ds is None: return None
    cols = [names[k] for k in sorted(names)]
    rows = [l.split()[:len(cols)] for l in txt[ds:] if l.strip() and len(l.split()) >= len(cols)]
    return pd.DataFrame(rows, columns=cols).apply(pd.to_numeric, errors="coerce")


def pick(df, cands):
    for c in cands:
        if c in df.columns: return c
    for c in df.columns:
        if any(x.lower() in c.lower() for x in cands): return c
    return None


def describe_filter_psa(psa):
    """Read the PSA and return (list of (variable, filter-label), summary string)."""
    root = ET.parse(psa).getroot()
    def tc(tag):
        el = root.find(tag)
        return float(el.get("value")) if el is not None else None
    tcA, tcB = tc("TimeConstFilterA"), tc("TimeConstFilterB")
    names = [ci.find("Calc/FullName").get("value") for ci in root.find("CalcArray")]
    types = [int(ai.get("value")) for ai in root.find("FilterTypeArray")]
    label = {0: "none", 1: f"low-pass A ({tcA:g} s)", 2: f"low-pass B ({tcB:g} s)"}
    rows = [(n, label.get(t, f"type {t}")) for n, t in zip(names, types)]
    applied = [f"{n} -> {lab}" for n, lab in rows if not lab.startswith("none")]
    summ = "; ".join(applied) if applied else "nothing (all variables set to no filter)"
    return rows, summ


def run_filter(in_cnv, stem):
    subprocess.run([str(SBE_BIN / "FilterW.exe"), f"/p{FILTER_PSA}", f"/i{in_cnv}",
                    f"/o{WORK}", f"/f{stem}", "/s"], capture_output=True, text=True)
    out = WORK / f"{stem}.cnv"
    return out if out.exists() else None


def dxdt_std(series):
    """Std of the scan-to-scan change of a series (the high-frequency content the filter removes)."""
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    x = x[np.isfinite(x)]
    if x.size < 3: return np.nan
    return float(np.nanstd(np.gradient(x)))


def plot_dpdt(cast, before, filt, pcol):
    """Before/after dP/dt vs pressure (the visual companion to the pressure reduction number)."""
    fig, ax = plt.subplots(figsize=(4.6, 6.0))
    for df, c, lab in [(before, C_BEFORE, "before"), (filt, C_FILTER, "filtered")]:
        p = pd.to_numeric(df[pcol], errors="coerce")
        ax.plot(np.gradient(p.to_numpy(float)), p, lw=0.5, color=c, label=lab, alpha=0.9)
    ax.invert_yaxis()
    ax.set_xlabel("dP/dt (db/scan)"); ax.set_ylabel("pressure (db)")
    ax.set_title(f"{cast}: pressure dP/dt, before vs filtered")
    ax.legend(loc="upper right", fontsize=8); _despine(ax)
    out = QC_PLOT_DIR / f"{cast}__qc_filter_dpdt.png"
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    return out


# ---------------------------------------------------------------- run
rows_psa, filt_summary = describe_filter_psa(FILTER_PSA)
print(f"Filter PSA: {FILTER_PSA.name}")
print(f"  applied: {filt_summary}\n")

summary_rows = []
hdr = f"{'cast':<16} {'variable':<13} {'std before':>12} {'std filtered':>13} {'reduction':>11}"
print(hdr); print("-" * len(hdr))
for cast in CASTS:
    aligned = ALIGN_DIR / f"{cast}_al.cnv"
    if not aligned.exists():
        print(f"{cast}: no aligned cnv"); continue
    before = read_cnv(aligned)
    filt = read_cnv(run_filter(aligned, f"{cast}_filt"))
    if before is None or filt is None:
        print(f"{cast}: filter run failed"); continue
    pcol = pick(before, VARS[0][0])
    for cands, label in VARS:
        b = pick(before, cands); f = pick(filt, cands)
        if b is None or f is None:
            continue
        sb, sf = dxdt_std(before[b]), dxdt_std(filt[f])
        red = 100 * (1 - sf / sb) if (sb and np.isfinite(sb) and sb > 0) else np.nan
        print(f"{cast:<16} {label:<13} {sb:>12.4f} {sf:>13.4f} {red:>10.1f}%")
        summary_rows.append({"cast": cast, "variable": label,
                             "dxdt_std_before": round(sb, 5), "dxdt_std_filtered": round(sf, 5),
                             "reduction_pct": round(red, 1) if np.isfinite(red) else ""})
    if pcol is not None:
        png = plot_dpdt(cast, before, filt, pcol)
        print(f"{'':<16} plot -> {png.name}")

if summary_rows:
    sdf = pd.DataFrame(summary_rows)
    out_csv = QC_PLOT_DIR / "qc_filter_dpdt_summary.csv"
    sdf.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"\nsummary -> {out_csv}")

print("\nReduction% = 100 * (1 - std(dX/dt)_filtered / std(dX/dt)_before).")
print("A small pressure reduction is expected on a 19plus V2 (strain-gauge): little dithering to")
print("remove. Variables the PSA does not filter show ~0% (their before and filtered are identical).")

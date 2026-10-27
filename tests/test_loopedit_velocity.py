"""
test_loopedit_velocity.py

Diagnostic for Loop Edit over-flagging on shallow shelf casts.
Runs Loop Edit on the filtered cnv, then reports the flagged-scan % two ways:
  (a) full cast  - includes surface soak + upcast (what over-flags)
  (b) downcast only, soak trimmed at SOAK_DEPTH_DB - what actually matters

If (b) drops into the ~2-10% range, sectioning the downcast is the fix (not the velocity
threshold). Compares 0.25 and 0.10 m/s minimum velocity for completeness.

This is a LOCAL diagnostic script, not a unit test: it runs the Sea-Bird software on real
cruise data. The executable part is under `if __name__ == "__main__":` so pytest can import
the module without launching Sea-Bird (so CI, and any machine without SBE, are unaffected).

USAGE: python test_loopedit_velocity.py
"""
from pathlib import Path
import re, subprocess
import pandas as pd, numpy as np

CTD_ROOT = Path(r"C:\Projects\ctd_pipeline")
CRUISE   = "P45_06"
CAST     = "P45_06_CTD_03"
SOAK_DEPTH_DB = 2.0                    # Leah's soak depth: start downcast past this
SBE_BIN  = Path(r"C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32")
PSA_DIR  = CTD_ROOT / "psa"
ALIGN    = CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "Align_CTD" / f"{CAST}_al.cnv"
FILTER_PSA = PSA_DIR / "05_filter.psa"
WORK     = CTD_ROOT / "_sbe_work" / "qc"

TESTS = [
    ("fixedmin  0.25",    PSA_DIR / "06b_loopedit_fixedmin.psa"),
    ("fixedmin  0.10",    PSA_DIR / "06b_loopedit_fixedmin_v010.psa"),
    ("percentmean 0.25",  PSA_DIR / "06a_loopedit_percentmean.psa"),
    ("percentmean 0.10",  PSA_DIR / "06a_loopedit_percentmean_v010.psa"),
]

def run(exe, psa, in_cnv, stem):
    subprocess.run([str(SBE_BIN/exe), f"/p{psa}", f"/i{in_cnv}", f"/o{WORK}", f"/f{stem}", "/s"],
                   capture_output=True, text=True)
    out = WORK / f"{stem}.cnv"
    return out if out.exists() else None

def read_cnv(cnv):
    txt=cnv.read_text(encoding="latin-1",errors="replace").splitlines()
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

def section_downcast(df, soak=SOAK_DEPTH_DB):
    p=pcol_of(df)
    pr=pd.to_numeric(df[p],errors="coerce")
    imax=int(pr.idxmax())
    down=df.loc[:imax]
    prd=pd.to_numeric(down[p],errors="coerce")
    below=prd[prd>soak]
    return down.loc[int(below.index[0]):] if not below.empty else down

def flag_pct(df):
    if "flag" not in df.columns or len(df)==0: return None
    f=pd.to_numeric(df["flag"],errors="coerce").fillna(0)!=0
    return 100*f.mean(), int(f.sum()), len(df)


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    filt = run("FilterW.exe", FILTER_PSA, ALIGN, f"{CAST}_filt")
    src = filt if filt else ALIGN
    print(f"Cast {CAST}   soak trim: pressure > {SOAK_DEPTH_DB} db\n")
    print(f"{'method / minV':<20} {'full cast':>12} {'downcast-only':>16}")
    print("-"*52)
    for label, psa in TESTS:
        if not psa.exists():
            print(f"{label:<20}  MISSING: {psa.name}"); continue
        out = run("LoopEditW.exe", psa, src, f"{CAST}_le_test")
        if out is None:
            print(f"{label:<20}  run failed"); continue
        df = read_cnv(out)
        full = flag_pct(df)
        dwn  = flag_pct(section_downcast(df))
        fs = f"{full[0]:.1f}% ({full[1]}/{full[2]})" if full else "n/a"
        ds = f"{dwn[0]:.1f}% ({dwn[1]}/{dwn[2]})" if dwn else "n/a"
        print(f"{label:<20} {fs:>12} {ds:>16}")
    print("\nRead the downcast-only column: if it is ~2-10%, sectioning fixes the over-flagging.")
    print("The full-cast column stays high because it still counts the soak + upcast flags.")


if __name__ == "__main__":
    main()

"""
test_loopedit_velocity.py

Quick diagnostic: run Loop Edit at two minimum-velocity thresholds (0.25 and 0.10 m/s) on one
cast, print the flagged-scan percentage for each. Tests whether lowering the velocity threshold
fixes the over-flagging seen on shallow shelf casts.

USAGE: python test_loopedit_velocity.py   (edit CAST below)
"""
from pathlib import Path
import re, subprocess
import pandas as pd, numpy as np

CTD_ROOT = Path(r"C:\Projects\ctd_pipeline")
CRUISE   = "P45_06"
CAST     = "P45_06_CTD_03"          # the shallow cast that over-flagged
SBE_BIN  = Path(r"C:\Program Files (x86)\Sea-Bird\SBEDataProcessing-Win32")
PSA_DIR  = CTD_ROOT / "psa"
ALIGN    = CTD_ROOT / "cruises" / CRUISE / "L2" / "CTD" / "Align_CTD" / f"{CAST}_al.cnv"
FILTER_PSA = PSA_DIR / "05_filter.psa"
WORK     = CTD_ROOT / "_sbe_work" / "qc"; WORK.mkdir(parents=True, exist_ok=True)

# (label, psa file) pairs to compare
TESTS = [
    ("fixedmin  0.25", PSA_DIR / "06b_loopedit_fixedmin.psa"),
    ("fixedmin  0.10", PSA_DIR / "06b_loopedit_fixedmin_v010.psa"),
    ("percentmean 0.25", PSA_DIR / "06a_loopedit_percentmean.psa"),
    ("percentmean 0.10", PSA_DIR / "06a_loopedit_percentmean_v010.psa"),
]

def run(exe, psa, in_cnv, stem):
    cmd = [str(SBE_BIN / exe), f"/p{psa}", f"/i{in_cnv}", f"/o{WORK}", f"/f{stem}", "/s"]
    subprocess.run(cmd, capture_output=True, text=True)
    out = WORK / f"{stem}.cnv"
    return out if out.exists() else None

def flag_pct(cnv):
    txt = cnv.read_text(encoding="latin-1", errors="replace").splitlines()
    names, ds = {}, None
    nr = re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*[:=]", re.I)
    for i, l in enumerate(txt):
        s=l.strip()
        if s.startswith("#") or s.startswith("*"):
            m=nr.search(l)
            if m: names[int(m.group(1))]=m.group(2).strip()
            if s.upper().startswith("*END*"): ds=i+1; break
    cols=[names[k] for k in sorted(names)]
    rows=[l.split()[:len(cols)] for l in txt[ds:] if l.strip() and len(l.split())>=len(cols)]
    df=pd.DataFrame(rows,columns=cols).apply(pd.to_numeric,errors="coerce")
    if "flag" not in df.columns: return None
    f=pd.to_numeric(df["flag"],errors="coerce").fillna(0)!=0
    return 100*f.mean(), int(f.sum()), len(df)

# filter first (loop edit runs on filtered cnv, per Leah's order)
filt = run("FilterW.exe", FILTER_PSA, ALIGN, f"{CAST}_filt")
src = filt if filt else ALIGN

print(f"Cast {CAST}  (filtered input: {'yes' if filt else 'no, using aligned'})\n")
print(f"{'method / minV':<20} {'flagged %':>10} {'flagged':>9} {'total':>8}")
print("-"*50)
for label, psa in TESTS:
    if not psa.exists():
        print(f"{label:<20}  MISSING PSA: {psa.name}"); continue
    out = run("LoopEditW.exe", psa, src, f"{CAST}_le_test")
    if out is None:
        print(f"{label:<20}  run failed"); continue
    r = flag_pct(out)
    if r: print(f"{label:<20} {r[0]:>9.1f}% {r[1]:>9} {r[2]:>8}")
print("\nExpectation: a well-set threshold flags roughly 2-10%. If 0.10 drops the")
print("percentage sharply toward that range, the 0.25 default was too aggressive for")
print("shelf descent. If it stays high, the soak/upcast are the issue -> section downcast.")

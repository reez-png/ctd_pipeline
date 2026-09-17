"""
step03_run_report.py

Read the CSV audit files step03 writes and produce ONE human-readable run report:

    <CRUISE>_step03_run_report.txt

It captures what step03 currently only prints to the console and then loses:
the RMSE progression across stages, the bottle count and its nitrate-range
coverage, and a provisional/OK verdict with the reasons. Point it at a cruise's
_audit + bottle folders after running step03.

The Stage 2 / Stage 3 fit coefficient dicts are printed by step03 but not saved
to CSV, so this report reconstructs the coverage/count diagnostics from the bottle
table and the assessment CSV, and leaves a clearly-labelled slot for the fit
coefficients if you paste them in (or we add a one-line dump to step03 later).

USAGE:
    python step03_run_report.py            # uses CRUISE_ID below
"""

from __future__ import annotations
from datetime import datetime
from pathlib import Path
import pandas as pd

CTD_ROOT   = Path(r"C:\Projects\ctd_pipeline")
CRUISE_ID  = "P45_06"
CRUISE_DIR = CTD_ROOT / "cruises" / CRUISE_ID
AUDIT      = CRUISE_DIR / "_audit"
BOTTLE     = CRUISE_DIR / "metadata" / "bottle_nitrate" / f"{CRUISE_ID}_bottle_nitrate.csv"

# Zheng (2024) reference points, for context in the report.
ZHENG_RMSE_LOW, ZHENG_RMSE_HIGH = 0.34, 0.78
ZHENG_TARGET_BOTTLES = 50


def coverage(nitrate: pd.Series) -> dict:
    n = pd.to_numeric(nitrate, errors="coerce").dropna()
    return {"low(0-2)": int((n < 2).sum()),
            "mid(2-8)": int(((n >= 2) & (n < 8)).sum()),
            "high(8+)": int((n >= 8).sum()),
            "total": int(len(n))}


def main() -> None:
    lines = []
    W = lines.append
    W(f"step03 nitrate-QC run report  -  {CRUISE_ID}")
    W("=" * 52)
    W(f"generated: {datetime.now():%Y-%m-%d %H:%M}")
    W("")

    # --- bottle table: count + coverage ---
    if BOTTLE.exists():
        b = pd.read_csv(BOTTLE)
        cov = coverage(b["nitrate_uM"])
        casts = sorted(b["cast_id"].dropna().unique())
        flagged = b["flag"].notna() & (b["flag"].astype(str).str.strip() != "")
        W("BOTTLE TABLE")
        W(f"  bottles (post-QC)      : {cov['total']}")
        W(f"  casts represented      : {len(casts)}  {casts}")
        W(f"  nitrate-range coverage : low(0-2)={cov['low(0-2)']}, "
          f"mid(2-8)={cov['mid(2-8)']}, high(8+)={cov['high(8+)']}")
        if flagged.any():
            W(f"  flagged bottles        : {int(flagged.sum())} "
              f"({sorted(b.loc[flagged,'flag'].unique())})")
        W("")
    else:
        cov = {"total": 0}
        W("BOTTLE TABLE: none found -> Stages 2-3 did not run.\n")

    # --- assessment RMSE progression ---
    rmse_path = AUDIT / "step03_assessment_rmse.csv"
    if rmse_path.exists():
        r = pd.read_csv(rmse_path)
        W("ASSESSMENT (RMSE vs bottles, uM)")
        for _, row in r.iterrows():
            W(f"  {row['stage']:<28} {float(row['rmse_uM']):6.2f}")
        # is each stage improving?
        vals = pd.to_numeric(r["rmse_uM"], errors="coerce").tolist()
        improving = all(vals[i] >= vals[i+1] for i in range(len(vals)-1))
        W(f"  monotonic improvement  : {'yes' if improving else 'no'}")
        final = vals[-1] if vals else float('nan')
        W(f"  Zheng reference        : {ZHENG_RMSE_LOW}-{ZHENG_RMSE_HIGH} uM")
        W("")
    else:
        final = float('nan')
        W("ASSESSMENT: step03_assessment_rmse.csv not found.\n")

    # --- per-cast QC medians ---
    qc_path = AUDIT / "step03_qc_summary.csv"
    if qc_path.exists():
        q = pd.read_csv(qc_path)
        W("PER-CAST QC MEDIANS (no3_qc, uM)")
        for _, row in q.iterrows():
            W(f"  {row['cast_id']:<16} median={float(row['median_qc_no3']):8.2f}  n={int(row['n'])}")
        W("")

    # --- verdict, stated plainly ---
    W("VERDICT")
    reasons = []
    if cov.get("total", 0) < ZHENG_TARGET_BOTTLES:
        reasons.append(f"only {cov.get('total',0)} bottles (Zheng targets ~{ZHENG_TARGET_BOTTLES})")
    if cov.get("total", 0) and (cov.get("mid(2-8)", 0) + cov.get("high(8+)", 0)) < 4:
        reasons.append("bottles skewed to low nitrate (few mid/high points to constrain the slope)")
    if pd.notna(final) and final > ZHENG_RMSE_HIGH:
        reasons.append(f"final RMSE {final:.2f} uM is well above Zheng's {ZHENG_RMSE_LOW}-{ZHENG_RMSE_HIGH} uM")
    if reasons:
        W("  PROVISIONAL - pipeline validated end to end, calibration under-constrained:")
        for r_ in reasons:
            W(f"    - {r_}")
        W("  Treat the calibrated no3_qc as a pipeline test, not a definitive product.")
    else:
        W("  Calibration is well-constrained by the available bottles.")
    W("")
    fit_path = AUDIT / "step03_fit_coefficients.csv"
    W("FIT COEFFICIENTS")
    if fit_path.exists():
        fc = pd.read_csv(fit_path)
        for _, row in fc.iterrows():
            parts = []
            for k in row.index:
                if k == "stage" or pd.isna(row[k]) or str(row[k]).strip() == "":
                    continue
                parts.append(f"{k}={row[k]}")
            W(f"  {row['stage']}:")
            for p in parts:
                W(f"      {p}")
    else:
        W("  (step03_fit_coefficients.csv not found - re-run step03 with the")
        W("   fit-coefficient dump enabled to capture c1/c2 and a1/a0 here.)")

    out = AUDIT / f"{CRUISE_ID}_step03_run_report.txt"
    out.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nReport -> {out}")


if __name__ == "__main__":
    main()

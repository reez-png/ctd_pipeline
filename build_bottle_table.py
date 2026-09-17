"""
build_bottle_table.py

Pair discrete nutrient bottle samples to CTD/SUNA casts by time, with a stated
confidence and reasoning for every pairing, and write:

  1. <CRUISE>_bottle_pairing_audit.csv  -- every station, its matched cast, the
     time gap, whether that cast has SUNA, the confidence, and the reason.
  2. <CRUISE>_bottle_nitrate.csv        -- the step03-ready bottle table, containing
     ONLY bottles whose cast has a merged SUNA product (others cannot be calibrated
     against SUNA and are excluded, with the reason logged in the audit).

Why time-based pairing: field nutrient logs identify samples by STATION (J1, J2, ...),
not by the CTD cast NUMBER that the .hex files were later given. The reliable link is
the sample time vs the cast .cnv start_time. This tool makes that link explicit and
auditable rather than hand-entered.

Confidence rules (stated, so they can be reviewed):
  - high    : nearest cast start within HIGH_TOL_MIN of the sample time
  - medium  : within MED_TOL_MIN
  - low     : nearest match is farther than MED_TOL_MIN (flagged, not trusted)
A pairing can be high-confidence yet still be EXCLUDED from the bottle table if that
cast has no SUNA data -- the audit records both facts separately.

USAGE (edit the settings, then run):
    python build_bottle_table.py
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
CTD_ROOT   = Path(r"C:\Projects\ctd_pipeline")
CRUISE_ID  = "P45_06"
CRUISE_DIR = CTD_ROOT / "cruises" / CRUISE_ID

# Raw multi-cruise nutrient file (Pierre's OMI export).
NUTRIENT_FILE = CRUISE_DIR / "metadata" / "bottle_nitrate" / "OMI_nutrients_all_raw.csv"
# How the raw file names the cruise (it uses hyphens: "P45-06").
NUTRIENT_CRUISE_LABEL = "P45-06"

L1_CTD_ROOT  = CRUISE_DIR / "L1" / "CTD"        # cast .cnv headers (start_time)
L2_SUNA_ROOT = CRUISE_DIR / "L2" / "SUNA"       # which casts have SUNA merged products
BOTTLE_DIR   = CRUISE_DIR / "metadata" / "bottle_nitrate"

# Column names in the raw nutrient file.
COL_CRUISE, COL_STATION, COL_TIME = "cruise", "station", "date_time"
COL_DEPTH, COL_NITRATE, COL_REP   = "depth", "nitrate_nitrite", "replicate"

# Time zone: Ghana is UTC+0 (GMT), and the field logs / nutrient file read as local
# time that equals UTC here. If a future cruise logs a real offset, set it (minutes).
SAMPLE_TZ_OFFSET_MIN = 0

# Confidence tolerances (minutes between sample time and cast start_time).
HIGH_TOL_MIN = 20
MED_TOL_MIN  = 90


def read_cast_start_times() -> dict:
    """cast_id -> start_time (UTC) parsed from each L1 .cnv header."""
    out = {}
    pat = re.compile(r"#\s*start_time\s*=\s*([A-Za-z]{3}\s+\d{1,2}\s+\d{4}\s+\d{2}:\d{2}:\d{2})")
    for cnv in sorted(L1_CTD_ROOT.glob(f"{CRUISE_ID}_CTD_*.cnv")):
        cast_id = cnv.stem  # e.g. P45_06_CTD_08
        for line in cnv.read_text(encoding="latin-1", errors="replace").splitlines():
            m = pat.search(line)
            if m:
                out[cast_id] = pd.Timestamp(datetime.strptime(m.group(1), "%b %d %Y %H:%M:%S"), tz="UTC")
                break
    return out


def casts_with_suna() -> set:
    return {p.name[: -len("_SUNA_1s.csv")] for p in L2_SUNA_ROOT.glob(f"{CRUISE_ID}_CTD_*_SUNA_1s.csv")}


def confidence(gap_min: float) -> str:
    if gap_min <= HIGH_TOL_MIN:
        return "high"
    if gap_min <= MED_TOL_MIN:
        return "medium"
    return "low"


def main() -> None:
    # --- load nutrients for this cruise ---
    raw = pd.read_csv(NUTRIENT_FILE)
    nut = raw[raw[COL_CRUISE].astype(str).str.strip() == NUTRIENT_CRUISE_LABEL].copy()
    if nut.empty:
        raise SystemExit(f"No rows for {NUTRIENT_CRUISE_LABEL} in {NUTRIENT_FILE.name}")
    nut["_t"] = (pd.to_datetime(nut[COL_TIME], errors="coerce")
                 - pd.Timedelta(minutes=SAMPLE_TZ_OFFSET_MIN)).dt.tz_localize("UTC")

    cast_times = read_cast_start_times()
    suna = casts_with_suna()
    if not cast_times:
        raise SystemExit(f"No cast start_times found under {L1_CTD_ROOT}")

    # --- pair each STATION (one sampling event) to the nearest cast by time ---
    audit_rows = []
    station_to_cast = {}
    for station, grp in nut.groupby(COL_STATION):
        stime = grp["_t"].dropna().iloc[0] if grp["_t"].notna().any() else pd.NaT
        best_cast, best_gap = None, float("inf")
        for cast_id, ct in cast_times.items():
            if pd.isna(stime):
                continue
            gap = abs((stime - ct).total_seconds()) / 60.0
            if gap < best_gap:
                best_gap, best_cast = gap, cast_id
        conf = confidence(best_gap) if best_cast else "none"
        has_suna = best_cast in suna
        review_flag = ""   # written into the bottle table's flag column for non-high matches
        if has_suna and conf == "high":
            reason = f"matched {best_cast} at {best_gap:.0f} min; cast has SUNA -> usable (high confidence)"
            usable = True
        elif has_suna and conf == "medium":
            reason = f"matched {best_cast} at {best_gap:.0f} min (medium confidence, > {HIGH_TOL_MIN} min); cast has SUNA -> KEPT but FLAGGED for review"
            usable = True
            review_flag = f"REVIEW_PAIRING_{best_gap:.0f}min"
        elif not has_suna:
            reason = f"matched {best_cast} at {best_gap:.0f} min, but that cast has NO SUNA product -> excluded (cannot calibrate against SUNA)"
            usable = False
        else:  # low confidence
            reason = f"nearest cast {best_cast} is {best_gap:.0f} min away (> {MED_TOL_MIN} min) -> too weak to trust -> excluded"
            usable = False
        station_to_cast[station] = (best_cast, usable, review_flag)
        audit_rows.append({
            "cruise_id": CRUISE_ID, "station": station,
            "sample_time_utc": stime, "matched_cast": best_cast,
            "cast_start_utc": cast_times.get(best_cast),
            "time_gap_min": round(best_gap, 1) if best_cast else None,
            "confidence": conf, "cast_has_suna": has_suna,
            "usable_for_calibration": usable, "n_bottles": len(grp), "reason": reason,
        })

    audit = pd.DataFrame(audit_rows).sort_values("station")
    BOTTLE_DIR.mkdir(parents=True, exist_ok=True)
    audit_path = BOTTLE_DIR / f"{CRUISE_ID}_bottle_pairing_audit.csv"
    audit.to_csv(audit_path, index=False, encoding="utf-8-sig")

    # --- build the step03 bottle table from USABLE stations only ---
    keep_rows = []
    for _, r in nut.iterrows():
        cast_id, usable, review_flag = station_to_cast.get(r[COL_STATION], (None, False, ""))
        if not usable:
            continue
        keep_rows.append({
            "cruise_id": CRUISE_ID,
            "cast_id": cast_id,
            "station": r[COL_STATION],
            "depth_m": r[COL_DEPTH],
            "time_utc": r["_t"].isoformat() if pd.notna(r["_t"]) else "",
            "nitrate_uM": r[COL_NITRATE],
            "replicate_id": f"{r[COL_STATION]}-{r[COL_DEPTH]:g}m" if pd.notna(r[COL_DEPTH]) else r[COL_STATION],
            "flag": review_flag,
        })
    bottle = pd.DataFrame(keep_rows)
    bottle_path = BOTTLE_DIR / f"{CRUISE_ID}_bottle_nitrate.csv"
    bottle.to_csv(bottle_path, index=False, encoding="utf-8-sig")

    # --- report ---
    print(f"Cruise {CRUISE_ID}: {len(nut)} nutrient bottles across {nut[COL_STATION].nunique()} stations")
    print(f"Casts with SUNA: {sorted(suna)}\n")
    print("Pairing audit:")
    for _, r in audit.iterrows():
        mark = "USE " if r["usable_for_calibration"] else "skip"
        print(f"  [{mark}] {r['station']:>3} -> {r['matched_cast']}  "
              f"({r['time_gap_min']} min, {r['confidence']}, "
              f"SUNA={r['cast_has_suna']})")
    print(f"\nUsable stations: {audit['usable_for_calibration'].sum()} of {len(audit)}")
    print(f"Bottle table rows written: {len(bottle)}")
    print(f"\nAudit  -> {audit_path}")
    print(f"Bottle -> {bottle_path}")
    if len(bottle) < 20:
        print(f"\nNOTE: only {len(bottle)} bottles are usable. This is below the ~50-bottle,\n"
              "three-subrange guidance in Zheng (2024). Stages 2-3 will run, but the\n"
              "cruise-specific bias fit is under-constrained; treat the calibrated\n"
              "product as provisional / a pipeline test, not a definitive calibration.")


if __name__ == "__main__":
    main()

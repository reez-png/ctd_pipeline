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

A medium-confidence pairing is KEPT, with the reason in the bottle table's `pairing_note` column. Until Nov 2026 it
was written into `flag` (REVIEW_PAIRING_<gap>min); since step03_nitrate_qc, step02b and the SUNA bottle tests leave
out any row whose `flag` is non-empty, `flag` is now left blank here and means only "exclude this bottle" (BDL, NES,
a suspected misfire, set by hand). To exclude a medium pairing after checking the field log, put a reason in `flag`.
A flag set by hand is kept when the table is rebuilt (matched by station, depth and nitrate value); one that no
longer matches a row is printed as a warning.

USAGE (edit the settings, then run):
    python build_bottle_table.py
"""

from __future__ import annotations

import re
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
from project_paths import project_root, product_dir   # beside this script, at the repo root
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))   # suncat, also when it is not installed
from suncat.stops import cast_stops, match_depths, TOL_M, TOL_FRAC  # noqa: E402
CTD_ROOT   = project_root()
CRUISE_ID  = "P45_06"
CRUISE_DIR = CTD_ROOT / "cruises" / CRUISE_ID

# Raw multi-cruise nutrient file (Pierre's OMI export).
NUTRIENT_FILE = CRUISE_DIR / "metadata" / "bottle_nitrate" / "OMI_nutrients_all_raw.csv"
# How the raw file names the cruise (it uses hyphens: "P45-06").
NUTRIENT_CRUISE_LABEL = "P45-06"

L1_CTD_ROOT  = product_dir(CRUISE_ID, "ctd_datcnv")        # cast .cnv headers (start_time)
L2_SUNA_ROOT = product_dir(CRUISE_ID, "suna_l2")       # which casts have SUNA merged products
CTD_1S_ROOT  = product_dir(CRUISE_ID, "ctd_1s")        # fallback for the stops when the L1 .cnv has no up-cast

# Bottle depths: the sheet gives the target depth. With no firing record (no .afm, so no .btl), the depth where the
# rosette was held on the up-cast (suncat.stops) is the candidate for where each bottle closed:
#   "report"  list every stop and the stop each sheet depth would take; the bottle table keeps the sheet depths
#   "apply"   depth_m becomes the stop depth (sheet depth kept as depth_nominal_m, depth_source "stop" or "sheet");
#             a sheet depth with no stop within max(5 m, 15 %) keeps the sheet depth
#   "off"     neither
# P45_06 is "report": its up-casts held only at the bottom and at the surface before recovery (and once at 12 m),
# never at a mid-depth sheet depth, so the stops cannot be where those bottles closed (Nov 2026). That fits how they
# were taken: fired on the up-cast while the rosette was moving, with perhaps a few near-surface samples from a
# separate hand-deployed Niskin (Reez, Nov 2026). Its depths are the sheet's targets.
STOP_DEPTHS = "report"
BOTTLE_DIR   = CRUISE_DIR / "metadata" / "bottle_nitrate"

# Column names in the raw nutrient file.
COL_CRUISE, COL_STATION, COL_TIME = "cruise", "station", "date_time"
COL_DEPTH, COL_NITRATE, COL_REP   = "depth", "nitrate_nitrite", "replicate"
COL_LAT, COL_TEMP = "latitude", "temperature"     # optional: for pressure to depth, and a check on the stop

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


def _key(station, depth, nitrate) -> tuple:
    d = pd.to_numeric(pd.Series([depth, nitrate]), errors="coerce")
    return (str(station).strip(), round(float(d[0]), 3) if pd.notna(d[0]) else None,
            round(float(d[1]), 4) if pd.notna(d[1]) else None)


def carry_hand_flags(bottle: pd.DataFrame, bottle_path: Path) -> tuple[list[str], list[str]]:
    """A flag set by hand in the existing bottle table (a reason to exclude a bottle) survives a rebuild: it is copied
    onto the new row with the same station, sheet depth and nitrate value (the value tells replicates apart). A flag whose
    row is gone (a station re-paired, a value corrected in the nutrient file) is reported, not silently dropped."""
    if not bottle_path.exists() or bottle.empty:
        return [], []
    old = pd.read_csv(bottle_path, keep_default_na=False, dtype=str)
    if "flag" not in old.columns:
        return [], []
    old = old[old["flag"].str.strip() != ""]
    new_keys = {_key(r["station"], r["depth_m"], r["nitrate_uM"]): i for i, r in bottle.iterrows()}
    carried, lost = [], []
    for _, r in old.iterrows():
        k = _key(r.get("station", ""), r.get("depth_nominal_m") or r.get("depth_m", ""), r.get("nitrate_uM", ""))
        desc = f"{k[0]} {r.get('depth_m', '')} m, {r.get('nitrate_uM', '')} uM: {r['flag']}"
        if k in new_keys:
            bottle.loc[new_keys[k], "flag"] = r["flag"]
            carried.append(desc)
        else:
            lost.append(desc)
    return carried, lost


def apply_stop_depths(bottle: pd.DataFrame, lat: float) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """depth_m -> the up-cast stop nearest the sheet depth, per cast (suncat.stops). Adds depth_nominal_m,
    depth_source ("stop" or "sheet"), stop_dur_s, ctd_temp_stop_c and depth_note. Returns (bottle, every stop found,
    one summary line per cast)."""
    b = bottle.copy()
    b["depth_m"] = pd.to_numeric(b["depth_m"], errors="coerce").astype(float)
    b["depth_nominal_m"] = b["depth_m"]
    b["depth_source"] = "sheet"
    b["stop_dur_s"] = np.nan
    b["ctd_temp_stop_c"] = np.nan
    b["depth_note"] = ""
    all_stops, lines = [], []
    for cast in sorted(b["cast_id"].dropna().unique()):
        rows = b["cast_id"] == cast
        st, used = cast_stops([L1_CTD_ROOT / f"{cast}.cnv", CTD_1S_ROOT / f"{cast}_1s.cnv"], lat=lat)
        zs = sorted(b.loc[rows, "depth_nominal_m"].dropna().unique())
        m = match_depths(zs, st)
        st = st.assign(cast_id=cast, source=used, sample_depths="")
        for _, r in m.iterrows():
            sel = rows & (b["depth_nominal_m"] == r["nominal_m"])
            if int(r["stop"]) >= 0:
                b.loc[sel, ["depth_m", "stop_dur_s", "ctd_temp_stop_c"]] = [r["stop_depth_m"], r["dur_s"], r["temp_c"]]
                b.loc[sel, "depth_source"] = "stop"
                st.loc[int(r["stop"]), "sample_depths"] = f"{r['nominal_m']:g}"
            else:
                b.loc[sel, "depth_note"] = r["note"]
        all_stops.append(st)
        up = st[st["leg"] != "down"]
        hit = m[m["stop"] >= 0]
        legs = ", ".join(f"{g} {z:.1f} m ({d:.0f} s)" for g, z, d in zip(up["leg"], up["depth_m"], up["dur_s"]))
        s = (f"  {cast}: stops from the bottom up in {used}: {legs or 'none'}; "
             f"{len(hit)} of {len(m)} sheet depths at a stop")
        if len(hit):
            s += " (stop minus sheet: " + ", ".join(f"{z:g} m {d:+.1f}" for z, d in zip(hit["nominal_m"], hit["dz_m"])) + ")"
        miss = m[m["stop"] < 0]
        for _, r in miss.iterrows():
            s += f"\n      {r['nominal_m']:g} m kept at the sheet depth: {r['note']}"
        spare = up[up["sample_depths"] == ""]
        if len(spare):
            s += f"\n      stops no sample matched: " + ", ".join(f"{z:.1f} m ({d:.0f} s)" for z, d in
                                                              zip(spare["depth_m"], spare["dur_s"]))
        lines.append(s)
    stops = pd.concat(all_stops, ignore_index=True) if all_stops else pd.DataFrame()
    return b, stops, lines


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
        review_note = ""   # written into the bottle table's pairing_note column for medium-confidence matches
        if has_suna and conf == "high":
            reason = f"matched {best_cast} at {best_gap:.0f} min; cast has SUNA -> usable (high confidence)"
            usable = True
        elif has_suna and conf == "medium":
            reason = f"matched {best_cast} at {best_gap:.0f} min (medium confidence, > {HIGH_TOL_MIN} min); cast has SUNA -> KEPT, with a pairing_note for review"
            usable = True
            review_note = f"REVIEW_PAIRING_{best_gap:.0f}min: medium-confidence time pairing, check the field log"
        elif not has_suna:
            reason = f"matched {best_cast} at {best_gap:.0f} min, but that cast has NO SUNA product -> excluded (cannot calibrate against SUNA)"
            usable = False
        else:  # low confidence
            reason = f"nearest cast {best_cast} is {best_gap:.0f} min away (> {MED_TOL_MIN} min) -> too weak to trust -> excluded"
            usable = False
        station_to_cast[station] = (best_cast, usable, review_note)
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
        cast_id, usable, review_note = station_to_cast.get(r[COL_STATION], (None, False, ""))
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
            "flag": "",
            "pairing_note": review_note,
            "temp_logged_c": pd.to_numeric(pd.Series([r.get(COL_TEMP)]), errors="coerce").iloc[0],
        })
    bottle = pd.DataFrame(keep_rows)
    bottle_path = BOTTLE_DIR / f"{CRUISE_ID}_bottle_nitrate.csv"
    carried, lost = carry_hand_flags(bottle, bottle_path)
    stop_lines, with_stops = [], None
    if STOP_DEPTHS in ("report", "apply") and len(bottle):
        lat = float(pd.to_numeric(nut.get(COL_LAT), errors="coerce").median()) if COL_LAT in nut else np.nan
        with_stops, stops, stop_lines = apply_stop_depths(bottle, lat if np.isfinite(lat) else 5.0)
        stops_path = BOTTLE_DIR / f"{CRUISE_ID}_bottle_stops.csv"
        stops.to_csv(stops_path, index=False, encoding="utf-8-sig")
        if STOP_DEPTHS == "apply":
            bottle = with_stops
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
    for r in carried:
        print(f"  hand flag kept: {r}")
    for r in lost:
        print(f"  WARNING: hand flag NOT carried over (no matching row now; check the nutrient file): {r}")
    if with_stops is not None:
        n_stop = int((with_stops["depth_source"] == "stop").sum())
        print(f"\nUp-cast stops (STOP_DEPTHS = {STOP_DEPTHS!r}): {n_stop} of {len(with_stops)} rows have a stop within "
              f"max({TOL_M:g} m, {TOL_FRAC:.0%}) of the sheet depth; "
              + ("depth_m is the stop depth for those" if STOP_DEPTHS == "apply" else "the bottle table keeps the sheet depths"))
        for s in stop_lines:
            print(s)
        d = with_stops.loc[with_stops["depth_source"] == "stop"]
        if len(d):
            dz = d["depth_m"] - d["depth_nominal_m"]
            print(f"  stop minus sheet over those rows: median {dz.median():+.2f} m, range {dz.min():+.1f} to {dz.max():+.1f} m")
            tt = (d["temp_logged_c"] - d["ctd_temp_stop_c"]).dropna()
            if len(tt):
                print(f"  check: temperature logged with the sample minus CTD temperature at the stop: median "
                      f"{tt.median():+.3f} C, largest |difference| {tt.abs().max():.3f} C ({len(tt)} rows)")
        print(f"Stops  -> {stops_path}")
    print(f"\nAudit  -> {audit_path}")
    print(f"Bottle -> {bottle_path}")
    if len(bottle) < 20:
        print(f"\nNOTE: only {len(bottle)} bottles are usable. This is below the ~50-bottle,\n"
              "three-subrange guidance in Zheng (2024). Stages 2-3 will run, but the\n"
              "cruise-specific bias fit is under-constrained; treat the calibrated\n"
              "product as provisional / a pipeline test, not a definitive calibration.")


if __name__ == "__main__":
    main()

"""Build a cruise's position table for step06, cruises/<CRUISE>/metadata/stations.csv, from what the cruise already holds.

    python tools/make_stations_table.py P45_06
    python tools/make_stations_table.py P45_06 --write

One row per CTD cast (the L1 .cnv files): cast_id, station, date_time_utc_start, lat_deg, lon_deg, position_source.
  date_time_utc_start  the .cnv header's start_time (step06 groups casts into sections by this date)
  lat_deg, lon_deg     from, in this order:
     1. the .cnv header's NMEA position (a "NMEA Latitude" / "NMEA Longitude" line, written when the deck unit had a
        GPS feed): the fix at the start of the cast
     2. the discrete nutrient file (latitude, longitude per sample), through the station-to-cast pairing that
        build_bottle_table.py wrote (<CRUISE>_bottle_pairing_audit.csv): the median over the station's sample rows.
        These are the positions logged with the samples, so they are the cast's position only as far as the
        pairing is right; a medium-confidence pairing is named in position_source.
     A cast that neither gives is left blank and listed: fill it by hand from the cast sheet or field log.
When both give a position, the header fix is used and the distance between them is printed as a check.

Without --write the table is printed only. An existing stations.csv is not replaced without --force. The file is
cruise data and stays out of git (*.csv is ignored).
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src"))
import project_paths as pp  # noqa: E402

_START = re.compile(r"#\s*start_time\s*=\s*([A-Za-z]{3}\s+\d{1,2}\s+\d{4}\s+\d{2}:\d{2}:\d{2})")
_NMEA = re.compile(r"^\*+\s*NMEA\s+(Latitude|Longitude)\s*[=:]\s*(.+)$", re.I)


def parse_coord(text: str) -> float:
    """'05 29.37 N', '000 00.99 W', '5.4895', '-0.0165' -> decimal degrees (S and W negative); NaN if unreadable."""
    s = text.strip().upper()
    hemi = re.search(r"[NSEW]\s*$", s)
    sign = -1.0 if hemi and hemi.group(0).strip() in ("S", "W") else 1.0
    nums = [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?", s)]
    if not nums:
        return np.nan
    deg = abs(nums[0]) + (nums[1] / 60.0 if len(nums) > 1 else 0.0) + (nums[2] / 3600.0 if len(nums) > 2 else 0.0)
    if nums[0] < 0:
        sign = -1.0
    return sign * deg


def cast_headers(cnv_paths) -> pd.DataFrame:
    """cast_id, start (UTC), hdr_lat, hdr_lon from each .cnv header (up to *END*)."""
    rows = []
    for p in cnv_paths:
        start, lat, lon = pd.NaT, np.nan, np.nan
        for line in Path(p).read_text(encoding="latin-1", errors="replace").splitlines():
            if line.startswith("*END*"):
                break
            m = _START.search(line)
            if m and pd.isna(start):
                start = pd.Timestamp(datetime.strptime(m.group(1), "%b %d %Y %H:%M:%S"))
            m = _NMEA.match(line.strip())
            if m:
                v = parse_coord(m.group(2))
                if m.group(1).lower() == "latitude":
                    lat = v
                else:
                    lon = v
        if lat == 0.0 and lon == 0.0:                      # a GPS feed with no fix
            lat = lon = np.nan
        rows.append(dict(cast_id=Path(p).stem, start=start, hdr_lat=lat, hdr_lon=lon))
    return pd.DataFrame(rows, columns=["cast_id", "start", "hdr_lat", "hdr_lon"])


def _km(lat1, lon1, lat2, lon2):
    x = (np.asarray(lon2) - np.asarray(lon1)) * np.cos(np.radians((np.asarray(lat1) + np.asarray(lat2)) / 2)) * 111.320
    return np.hypot(x, (np.asarray(lat2) - np.asarray(lat1)) * 110.574)


def nutrient_positions(nut: pd.DataFrame, label: str, audit: pd.DataFrame) -> pd.DataFrame:
    """cast_id, station, nut_lat, nut_lon, n_rows, spread_km, confidence, gap_min: per paired station, the median of its
    sample rows' positions; spread_km is the largest distance of a row from that median (a typing error shows here)."""
    d = nut[nut["cruise"].astype(str).str.strip() == label].copy()
    d["lat"] = pd.to_numeric(d["latitude"], errors="coerce")
    d["lon"] = pd.to_numeric(d["longitude"], errors="coerce")
    d = d.dropna(subset=["lat", "lon"])
    a = audit.dropna(subset=["matched_cast"])
    a = a[a["confidence"].isin(["high", "medium"])]
    rows = []
    for _, r in a.iterrows():
        g = d[d["station"].astype(str).str.strip() == str(r["station"]).strip()]
        if g.empty:
            continue
        lat, lon = float(g["lat"].median()), float(g["lon"].median())
        rows.append(dict(cast_id=r["matched_cast"], station=r["station"], nut_lat=lat, nut_lon=lon, n_rows=len(g),
                         spread_km=float(np.max(_km(lat, lon, g["lat"], g["lon"]))),
                         confidence=r["confidence"], gap_min=r.get("time_gap_min")))
    out = pd.DataFrame(rows, columns=["cast_id", "station", "nut_lat", "nut_lon", "n_rows", "spread_km",
                                      "confidence", "gap_min"])
    dup = out["cast_id"].duplicated(keep=False)
    if dup.any():                                         # two stations on one cast: one row, no position
        one = out[dup].groupby("cast_id", as_index=False).agg(station=("station", lambda x: "+".join(map(str, x))))
        out = pd.concat([out[~dup], one.assign(nut_lat=np.nan, nut_lon=np.nan, n_rows=0, spread_km=np.nan,
                                               confidence="", gap_min=np.nan)], ignore_index=True)
    return out


def build(headers: pd.DataFrame, nutpos: pd.DataFrame) -> pd.DataFrame:
    t = headers.merge(nutpos, on="cast_id", how="left")
    lat, lon, src, check = [], [], [], []
    for _, r in t.iterrows():
        has_h = np.isfinite(r["hdr_lat"]) and np.isfinite(r["hdr_lon"])
        has_n = pd.notna(r.get("nut_lat")) and pd.notna(r.get("nut_lon"))
        if has_h:
            lat.append(r["hdr_lat"]); lon.append(r["hdr_lon"]); src.append("cnv header NMEA fix")
            check.append(float(_km(r["hdr_lat"], r["hdr_lon"], r["nut_lat"], r["nut_lon"])) if has_n else np.nan)
        elif has_n:
            lat.append(r["nut_lat"]); lon.append(r["nut_lon"]); check.append(np.nan)
            s = f"nutrient file, station {r['station']} ({int(r['n_rows'])} rows)"
            if r["confidence"] == "medium":
                s += f", medium-confidence pairing ({r['gap_min']} min): check the field log"
            src.append(s)
        else:
            lat.append(np.nan); lon.append(np.nan); check.append(np.nan)
            src.append("two stations paired to this cast: fill by hand" if pd.notna(r.get("station")) else "")
    t["lat_deg"], t["lon_deg"], t["position_source"], t["header_vs_nutrient_km"] = lat, lon, src, check
    t["station"] = t["station"].fillna("")
    t["date_time_utc_start"] = pd.to_datetime(t["start"]).dt.strftime("%Y-%m-%d %H:%M:%S")
    return t


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cruise")
    ap.add_argument("--label", help="how the nutrient file names the cruise (default: the cruise id with '-' for '_')")
    ap.add_argument("--nutrient-file", type=Path, help="default: metadata/bottle_nitrate/OMI_nutrients_all_raw.csv")
    ap.add_argument("--write", action="store_true", help="write metadata/stations.csv")
    ap.add_argument("--force", action="store_true", help="replace an existing stations.csv")
    a = ap.parse_args()

    meta = pp.cruise_dir(a.cruise) / "metadata"
    cnvs = sorted(pp.product_dir(a.cruise, "ctd_datcnv").glob(f"{a.cruise}_CTD_*.cnv"))
    if not cnvs:
        print(f"no {a.cruise}_CTD_*.cnv in {pp.product_dir(a.cruise, 'ctd_datcnv')}"); return 1
    headers = cast_headers(cnvs)

    nfile = a.nutrient_file or meta / "bottle_nitrate" / "OMI_nutrients_all_raw.csv"
    audit_p = meta / "bottle_nitrate" / f"{a.cruise}_bottle_pairing_audit.csv"
    nutpos = pd.DataFrame(columns=["cast_id", "station", "nut_lat", "nut_lon", "n_rows", "spread_km", "confidence",
                                   "gap_min"])
    if nfile.exists() and audit_p.exists():
        nut = pd.read_csv(nfile)
        if {"latitude", "longitude", "cruise", "station"} <= set(nut.columns):
            nutpos = nutrient_positions(nut, a.label or a.cruise.replace("_", "-"), pd.read_csv(audit_p))
        else:
            print(f"{nfile.name} has no latitude / longitude columns: header fixes only")
    else:
        print(f"no nutrient positions: needs {nfile.name} and {audit_p.name} (run build_bottle_table.py first)")

    t = build(headers, nutpos)
    show = t[["cast_id", "station", "date_time_utc_start", "lat_deg", "lon_deg", "position_source"]]
    with pd.option_context("display.width", 200, "display.max_colwidth", 90):
        print(show.round({"lat_deg": 5, "lon_deg": 5}).to_string(index=False))
    if (wide := nutpos[nutpos["spread_km"] > 1.0]).size:
        print("\nstations whose sample rows disagree by more than 1 km (check the nutrient file):")
        print(wide[["station", "cast_id", "n_rows", "spread_km"]].round(2).to_string(index=False))
    if t["header_vs_nutrient_km"].notna().any():
        c = t.dropna(subset=["header_vs_nutrient_km"])
        print(f"\nheader fix vs nutrient file, {len(c)} casts: median {c['header_vs_nutrient_km'].median():.2f} km, "
              f"largest {c['header_vs_nutrient_km'].max():.2f} km ({c.loc[c['header_vs_nutrient_km'].idxmax(), 'cast_id']})")
    miss = t.loc[t["lat_deg"].isna(), "cast_id"].tolist()
    print(f"\n{len(t) - len(miss)} of {len(t)} casts positioned" + (f"; no position: {', '.join(miss)} "
          "(fill lat_deg / lon_deg by hand; step06 leaves these casts out)" if miss else ""))

    if not a.write:
        print("\n(not written: add --write)"); return 0
    out = meta / "stations.csv"
    if out.exists() and not a.force:
        print(f"\n{out} exists: not replaced (add --force)"); return 1
    t.insert(0, "cruise_id", a.cruise)
    t[["cruise_id", "cast_id", "station", "date_time_utc_start", "lat_deg", "lon_deg", "position_source"]].to_csv(
        out, index=False, encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

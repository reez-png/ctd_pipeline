"""Write a synthetic cruise (no real data) built to differ from P45_06 / P45_07 wherever the pipeline could assume them.

    python tools/make_synthetic_cruise.py <root> [CRUISE_ID] [--lat0 -12.4] [--lon0 40.6] [--utc-offset 3]

Writes <root>/cruises/<CRUISE_ID>/ in the v2 layout. Point the pipeline at it with CTD_PIPELINE_ROOT=<root>; never
write it into the repo's own cruises folder. tests/test_synthetic_cruise.py uses it.

What differs from P45 on purpose:
  * cruise id not starting with P45; casts numbered 1..4 but named with three digits (<CRUISE>_CTD_001)
  * position set by --lat0/--lon0 (southern hemisphere east of Greenwich by default; near the equator for the
    step06 test), local clock UTC + --utc-offset
  * SBE 9plus pressure name prDM; ship heave makes pressure loops, flagged Loop Edit style; cast 3 has NO flag column
  * Bin Average files write the MEAN pressure of each bin (the P45_07 style), not the bin centre
  * bottles give the cast as a number for casts 1 and 3 and as a full cast id (hyphenated cruise) for casts 2 and 4
Known truth, so outputs can be checked:
  chlorophyll  chl(z) = 0.05 + 1.2 exp(-((z - 25) / 6)^2) mg/m3, fluorometer raw = 0.04 + 1.6 chl (+ noise sd 0.01)
  nitrate      no3(z) = 2 + 10 (1 + tanh((z - 40) / 8)) umol/L (+ noise sd 0.1)
  pressure     from depth by inverting UNESCO 1983 at the cast's latitude (the formula behind SBE depSM)
  cast starts  local clock 09, 12, 15 and 21 h
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

GAIN, OFFSET = 1.6, 0.04
# cast number, d_lat, d_lon (from lat0/lon0), max depth (m), local clock start hour, has Loop Edit flag column
CASTS = [(1, 0.00, 0.00, 60.0, 9, True), (2, -0.10, 0.15, 90.0, 12, True),
         (3, -0.22, 0.32, 130.0, 15, False), (4, -0.35, 0.50, 160.0, 21, True)]


def unesco_depth(p, lat):
    x = np.sin(np.radians(lat)) ** 2
    g = 9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * x) * x) + 1.092e-6 * p
    return ((((-1.82e-15 * p + 2.279e-10) * p - 2.2512e-5) * p + 9.72659) * p) / g


def unesco_pressure(z, lat):
    z = np.asarray(z, float); p = z.copy()
    for _ in range(6):
        p = p + (z - unesco_depth(p, lat)) / 0.99446
    return p


def chl_true(z):
    return 0.05 + 1.2 * np.exp(-((np.asarray(z, float) - 25.0) / 6.0) ** 2)


def no3_true(z):
    return 2.0 + 10.0 * (1.0 + np.tanh((np.asarray(z, float) - 40.0) / 8.0))


def _write_cnv(path: Path, df: pd.DataFrame, start: pd.Timestamp, extra: tuple[str, ...] = ()) -> None:
    desc = {"timeS": "Time, Elapsed [seconds]", "prDM": "Pressure, Digiquartz [db]", "tv290C": "Temperature [ITS-90, deg C]",
            "sal00": "Salinity, Practical [PSU]", "flECO-AFL": "Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]",
            "sigma-é00": "Density [sigma-theta, kg/m^3]", "depSM": "Depth [salt water, m]",
            "nbin": "number of scans per bin", "flag": "flag"}
    names = list(df.columns)
    head = ["* Sea-Bird SBE 9 Data File:", "* SYNTHETIC: tools/make_synthetic_cruise.py, not real data"]
    head += [f"# name {i} = {n}: {desc.get(n, n)}" for i, n in enumerate(names)]
    head += [f"# start_time = {start:%b %d %Y %H:%M:%S} [Instrument's time stamp, header]", "# interval = seconds: 0.25"]
    head += list(extra) + ["*END*"]
    body = [" ".join(("%11.3e" % v) if n == "flag" else ("%11.4f" % v) for n, v in zip(names, r))
            for r in df.itertuples(index=False)]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(head + body) + "\n", encoding="latin-1")


def make_cruise(root, cruise: str = "ZZ12_03", lat0: float = -12.4, lon0: float = 40.6, utc_offset: float = 3.0,
                seed: int = 7) -> Path:
    rng = np.random.default_rng(seed)
    C = Path(root) / "cruises" / cruise
    (C / "_audit").mkdir(parents=True, exist_ok=True)
    (C / "_audit" / "layout.json").write_text('{"layout": "v2"}')
    summary, chl_rows, no3_rows = [], [], []
    for num, dlat, dlon, zmax, hour, has_flag in CASTS:
        cast, lat, lon = f"{cruise}_CTD_{num:03d}", lat0 + dlat, lon0 + dlon
        dt = 0.25
        t_down = np.arange(0, zmax / 0.8, dt)                       # 0.8 m/s with 0.7 m, 8 s heave
        z_down = np.clip(2.0 + 0.8 * t_down + 0.7 * np.sin(2 * np.pi * t_down / 8.0), 0.5, None)
        z_down = z_down[z_down <= zmax]
        z_up = [np.full(int(20 / dt), z_down[-1])]                   # 20 s at the bottom
        for a, b in ((z_down[-1], 0.6 * zmax), (0.6 * zmax, 0.3 * zmax), (0.3 * zmax, 2.0)):
            z_up += [np.arange(a, b, -dt), np.full(int(30 / dt), b)]  # 1 m/s up, 30 s bottle stops
        z = np.concatenate([z_down] + z_up)
        t = np.arange(len(z)) * dt
        p = unesco_pressure(z, lat)
        nd = len(z_down)
        bad = np.zeros(len(z), bool)
        bad[:nd] = p[:nd] < np.maximum.accumulate(p[:nd]) - 1e-9      # Loop Edit: reversals on the down-cast
        df = pd.DataFrame({"timeS": t, "prDM": p, "tv290C": 26.0 - 0.08 * z, "sal00": 35.0 + 0.004 * z,
                           "flECO-AFL": OFFSET + GAIN * chl_true(z) + rng.normal(0, 0.01, len(z)),
                           "sigma-é00": 22.5 + 0.03 * z, "depSM": z})
        if has_flag:
            df["flag"] = np.where(bad, -9.99e-29, 0.0)
        start_utc = pd.Timestamp("2026-03-10") + pd.Timedelta(days=num - 1) + pd.Timedelta(hours=hour - utc_offset)
        _write_cnv(C / "L1/CTD/Derived" / f"{cast}_der.cnv", df, start_utc)
        kept = ~bad if has_flag else np.ones(len(z), bool)
        imax = int(np.argmax(np.where(kept, p, -np.inf)))
        k = df[kept & (np.arange(len(z)) <= imax)].drop(columns=["flag"], errors="ignore")
        b = np.floor(k["prDM"] + 0.5)
        g = k.groupby(b).mean(); g["nbin"] = k.groupby(b).size(); g["flag"] = 0.0
        _write_cnv(C / "L3/CTD/CTD_1m_down" / f"{cast}_1m_down.cnv", g.reset_index(drop=True), start_utc,
                   ("# binavg_bintype = decibars", "# binavg_binsize = 1", "# binavg_excl_bad_scans = yes"))
        idx = np.arange(0, len(z), 4)                                  # 1 Hz nitrate frames, both legs
        nf = pd.DataFrame({"utc_time": (start_utc + pd.to_timedelta(t[idx], unit="s")).strftime("%Y-%m-%d %H:%M:%S.%f"),
                           "temp_c": df["tv290C"].to_numpy()[idx], "salinity": df["sal00"].to_numpy()[idx],
                           "pressure_db": p[idx], "depth_m": z[idx], "leg": np.where(idx <= imax, "down", "up"),
                           "ctd_gap_s": 0.2, "ctd_scan_pressure_db": p[idx], "ctd_scan_kept": kept[idx],
                           "L1_clean": rng.random(len(idx)) > 0.1,
                           "no3_qc": no3_true(z[idx]) + rng.normal(0, 0.1, len(idx)), "lag_limited": False})
        (C / "L2/nitrate").mkdir(parents=True, exist_ok=True)
        nf.to_csv(C / "L2/nitrate" / f"{cast}_nitrate_QC.csv", index=False)
        summary.append(dict(cast=num, station=f"T{num}", date_time_utc_start=f"{start_utc:%Y-%m-%d %H:%M:%S}",
                            latitude_start=round(lat, 4), longitude_start=round(lon, 4)))
        local = start_utc + pd.Timedelta(hours=utc_offset)
        for zb in (5, 15, 25, 35, 50):
            if zb <= zmax - 3:
                chl_rows.append(dict(cruise=cruise.replace("_", "-"), station=f"T{num}",
                                     cast=(num if num % 2 else cast.replace("_", "-", 1)),
                                     date_time=f"{local:%m/%d/%Y %H:%M}", depth=zb, replicate="a",
                                     chl=round(float(chl_true(zb)) * (1 + rng.normal(0, 0.04)), 3)))
        for zb in (10, 40, 55):
            if zb <= zmax - 3:
                no3_rows.append(dict(cast_id=cast, depth_m=zb, nitrate_uM=round(float(no3_true(zb)), 2)))
    (C / "metadata/chla").mkdir(parents=True, exist_ok=True)
    (C / "metadata/bottle_nitrate").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary).to_csv(C / "metadata" / f"{cruise}_CTD_cast_summary.csv", index=False)
    pd.DataFrame(chl_rows).to_csv(C / "metadata/chla" / f"{cruise}_chla.csv", index=False)
    pd.DataFrame(no3_rows).to_csv(C / "metadata/bottle_nitrate" / f"{cruise}_bottle_nitrate.csv", index=False)
    return C


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("root"); ap.add_argument("cruise", nargs="?", default="ZZ12_03")
    ap.add_argument("--lat0", type=float, default=-12.4); ap.add_argument("--lon0", type=float, default=40.6)
    ap.add_argument("--utc-offset", type=float, default=3.0)
    a = ap.parse_args()
    root = Path(a.root).resolve()
    if (Path(__file__).resolve().parents[1] / "cruises") == root / "cruises":
        ap.error("refusing to write a synthetic cruise into the repo's own cruises folder; use another root")
    print("synthetic cruise ->", make_cruise(root, a.cruise, a.lat0, a.lon0, a.utc_offset))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

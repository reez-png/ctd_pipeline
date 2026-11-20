"""What is on the ECO's third channel (flECO-AFL1)? A per-cast check on the Derived .cnv files.

The ECO V2 on these cruises (S/N 10043) measures three things: chlorophyll fluorescence excited at 470 nm, at 435 nm,
and backscattering at 700 nm (the NTU sheet is a second calibration of the same 700 nm channel). It reaches the
SBE 19plus V2 on the RS-232 port (.xmlcon SerialRS232C_Sensor 10, WET Labs), which takes three channels:
    channel 1  FluoroWetlabECO_AFL_FL   scale 7.345e-3, dark 54   = the 470 nm sheet           -> flECO-AFL
    channel 2  TurbidityMeter           scale 1.676e-6, dark 55   = the 700 nm beta sheet      -> not exported
    channel 3  FluoroWetlabECO_AFL_FL   scale 10,       dark 55   = no sheet (placeholder)     -> flECO-AFL1
So flECO-AFL1 = 10 x (counts - 55), and the likeliest occupant is the 435 nm chlorophyll channel, the one measurement
the .xmlcon gives no coefficients for. If the ECO's channel order differs, it is the 700 nm backscatter instead.

The test: on kept down-cast scans below MIN_DBAR, rank correlations of the channel's counts with the 470 nm
chlorophyll (flECO-AFL) and with beam attenuation (CStarAt0, particles), each also with the other held fixed
(partial). A second chlorophyll channel follows flECO-AFL once particles are held fixed; a backscatter channel
follows CStarAt0 once chlorophyll is held fixed.

    python tools\\inspect_eco_channels.py P45_06 [P45_07 ...]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from project_paths import product_dir  # noqa: E402
from suncat.io.cnv import read_cnv  # noqa: E402

MIN_DBAR = 3.0
AFL = ("flECO-AFL", 7.345e-3, 54.0)       # column, scale, dark from the .xmlcon
AFL1 = ("flECO-AFL1", 10.0, 55.0)
CSTAR = "CStarAt0"


def counts(v, scale, dark):
    return np.asarray(v, float) / scale + dark


def _rank(x):
    return pd.Series(x).rank().to_numpy()


def partial_r(x, y, z):
    """Spearman correlation of x and y with z held fixed (ranks, then residuals of a straight line in z)."""
    rx, ry, rz = _rank(x), _rank(y), _rank(z)
    A = np.c_[np.ones_like(rz), rz]
    ex = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    ey = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1]) if ex.std() > 0 and ey.std() > 0 else np.nan


def one_cast(path: Path) -> dict:
    d = read_cnv(path, mask_badflag=True)
    r = dict(cast=path.stem, n=0)
    if d is None or d.empty or AFL1[0] not in d or AFL[0] not in d:
        r["note"] = "channel missing"
        return r
    pcol = next((c for c in d.columns if c.lower().startswith(("prdm", "prm"))), None)
    kept = d["flag"].eq(0) if "flag" in d else pd.Series(True, index=d.index)
    p = pd.to_numeric(d[pcol], errors="coerce")
    imax = int(np.nanargmax(p.where(kept).to_numpy()))
    m = kept & (d.index <= imax) & (p > MIN_DBAR)
    c1 = counts(d.loc[m, AFL1[0]], *AFL1[1:])
    c0 = counts(d.loc[m, AFL[0]], *AFL[1:])
    cs = pd.to_numeric(d.loc[m, CSTAR], errors="coerce").to_numpy() if CSTAR in d else np.full(c1.size, np.nan)
    ok = np.isfinite(c1) & np.isfinite(c0) & np.isfinite(cs)
    c1, c0, cs = c1[ok], c0[ok], cs[ok]
    r.update(n=int(ok.sum()), afl1_counts_min=float(np.min(c1)) if c1.size else np.nan,
             afl1_counts_med=float(np.median(c1)) if c1.size else np.nan,
             afl1_counts_max=float(np.max(c1)) if c1.size else np.nan,
             afl1_below_50=int((c1 < 50).sum()), afl1_unique=int(np.unique(np.round(c1, 3)).size))
    if r["afl1_unique"] <= 1:
        r["note"] = "flat"
        return r
    if c1.size < 20:
        r["note"] = "too few scans"
        return r
    r.update(r_chl470=float(pd.Series(c1).corr(pd.Series(c0), method="spearman")),
             r_cstar=float(pd.Series(c1).corr(pd.Series(cs), method="spearman")),
             r_chl470_given_cstar=partial_r(c1, c0, cs), r_cstar_given_chl470=partial_r(c1, cs, c0),
             note="")
    return r


def main(cruises):
    for cid in cruises:
        files = sorted(product_dir(cid, "ctd_derived").glob("*.cnv"))
        t = pd.DataFrame([one_cast(f) for f in files])
        print(f"\n{cid}: {len(files)} Derived .cnv; flECO-AFL1 as counts (value / 10 + 55), kept down-cast scans "
              f"below {MIN_DBAR:g} dbar; Spearman r")
        with pd.option_context("display.width", 220, "display.max_columns", 30):
            print(t.round(3).to_string(index=False))
        live = t[t.get("note", "").eq("")] if "note" in t else t
        if len(live) and "r_chl470_given_cstar" in live:
            a, b = live["r_chl470_given_cstar"].median(), live["r_cstar_given_chl470"].median()
            print(f"  median partial r: with 470 nm chlorophyll (particles fixed) {a:.2f}; "
                  f"with beam attenuation (chlorophyll fixed) {b:.2f}")


if __name__ == "__main__":
    main(sys.argv[1:] or ["P45_06"])

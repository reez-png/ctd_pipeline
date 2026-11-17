"""Write (or show) a cruise's SUNA clock offset, cruises/<CRUISE>/metadata/suna_time_offsets.csv.

    python tools/set_suna_time_offset.py P45_07 --offset -16.0 --source "water entry (step04): +16.1 s, ..."
    python tools/set_suna_time_offset.py P45_07 --show

The offset is ADDED to the SUNA time: minus the measured lag (SUNA 16 s late -> --offset -16.0). It is written for
every cast the cruise has merged SUNA data for (L2/SUNA/*_SUNA_1s.csv), or for --casts. step02_merge applies the
file whenever it exists; rerun step02_merge and step03 afterwards. An existing file is not replaced without --force.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src"))
import project_paths as pp  # noqa: E402
from suncat.suna_clock import describe, load_offsets, write_offsets  # noqa: E402


def cruise_casts(cruise: str) -> list[str]:
    """Casts with merged SUNA data, else casts in the co-located product."""
    l2 = pp.product_dir(cruise, "suna_l2")
    casts = [f.name.replace("_SUNA_1s.csv", "") for f in sorted(l2.glob("*_SUNA_1s.csv"))]
    if not casts:
        col = l2 / f"{cruise}__suna_L2_colocated.csv"
        if col.exists():
            import pandas as pd
            casts = sorted(pd.read_csv(col, usecols=["cast"])["cast"].dropna().astype(str).unique())
    return casts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cruise")
    ap.add_argument("--offset", type=float, help="seconds ADDED to the SUNA time (minus the measured lag)")
    ap.add_argument("--source", help="where the offset comes from (required)")
    ap.add_argument("--casts", nargs="+", help="casts to write (default: every cast with merged SUNA data)")
    ap.add_argument("--force", action="store_true", help="replace an existing file")
    ap.add_argument("--show", action="store_true", help="print the current file and stop")
    a = ap.parse_args()
    cdir = pp.cruise_dir(a.cruise)
    if a.show or a.offset is None:
        try:
            print(f"{a.cruise}: {describe(cdir)}")
            for c, o in sorted(load_offsets(cdir).items()):
                print(f"  {c}  {o:+g} s")
        except ValueError as e:
            print(e)
            return 1
        return 0
    casts = a.casts or cruise_casts(a.cruise)
    try:
        path = write_offsets(cdir, casts, a.offset, a.source, overwrite=a.force)
    except (FileExistsError, ValueError) as e:
        print(e)
        return 1
    print(f"{a.cruise}: {describe(cdir)}\n  wrote {path} ({len(casts)} casts)")
    print("  next: rerun step02_merge_ctd_suna and step03_nitrate_qc for this cruise (step02 applies the file whenever it exists)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

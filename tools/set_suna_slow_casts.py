"""Write (or show) a cruise's per-cast SUNA response flags, cruises/<CRUISE>/metadata/suna_response.csv.

    python tools/set_suna_slow_casts.py P45_07                  from step02b's per-cast tau, slow = tau >= 3 x median
    python tools/set_suna_slow_casts.py P45_07 --factor 2.5
    python tools/set_suna_slow_casts.py P45_07 --casts P45_07_CTD_01 P45_07_CTD_09 --source "..."   by hand
    python tools/set_suna_slow_casts.py P45_07 --show

Reads _audit/suna/<CRUISE>__suna_response_time.csv and the cruise median from __suna_response_summary.json (run
step02b first). Prints every cast's tau, sorted, against the threshold, so the gap the rule falls in is visible
before anything is written. step08 copies the flags into the L3 nitrate; rerun step08 and step06 afterwards. An
existing file is not replaced without --force.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src"))
import project_paths as pp  # noqa: E402
from suncat.suna_response import FACTOR, describe, flag_slow, load_response, write_response  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cruise")
    ap.add_argument("--factor", type=float, default=FACTOR, help=f"slow = tau >= factor x cruise median (default {FACTOR:g})")
    ap.add_argument("--casts", nargs="+", help="flag these casts as slow instead of applying the rule (needs --source)")
    ap.add_argument("--source", help="where the flags come from (default: composed from step02b and the rule)")
    ap.add_argument("--force", action="store_true", help="replace an existing file")
    ap.add_argument("--show", action="store_true", help="print the current file and stop")
    a = ap.parse_args(argv)
    cdir = pp.cruise_dir(a.cruise)
    if a.show:
        try:
            print(f"{a.cruise}: {describe(cdir)}")
            for c, (tau, slow) in sorted(load_response(cdir).items()):
                print(f"  {c}  tau {tau:5.1f} s  {'SLOW' if slow else ''}")
        except ValueError as e:
            print(e)
            return 1
        return 0
    audit = cdir / "_audit" / "suna"
    tab_f, sum_f = audit / f"{a.cruise}__suna_response_time.csv", audit / f"{a.cruise}__suna_response_summary.json"
    if not tab_f.exists() or not sum_f.exists():
        print(f"no {tab_f.name} or {sum_f.name} in {audit}: run step02b_suna_response_correction for {a.cruise} first")
        return 1
    tab = pd.read_csv(tab_f, encoding="utf-8-sig")
    med = json.loads(sum_f.read_text()).get("tau_median")
    if med is None:
        print(f"{sum_f.name} has no cruise median tau (fewer than 3 casts with an estimate): nothing to flag against")
        return 1
    flags = flag_slow(tab, float(med), a.factor)
    if a.casts:
        if not a.source:
            print("--casts needs --source: say why these casts are flagged")
            return 1
        unknown = sorted(set(a.casts) - set(flags["cast_id"]))
        if unknown:
            print(f"not in step02b's table: {unknown}")
            return 1
        flags["slow_response"] = flags["cast_id"].isin(a.casts)
    source = a.source or (f"step02b tau per cast; slow = tau >= {a.factor:g} x cruise median {float(med):.1f} s "
                          f"({a.factor * float(med):.1f} s) or at the grid's edge")
    thr = a.factor * float(med)
    print(f"{a.cruise}: tau per cast (step02b), cruise median {float(med):.1f} s, threshold {thr:.1f} s"
          + (" (not used: --casts given)" if a.casts else ""))
    for _, r in flags.sort_values("tau_s", ascending=False, na_position="last").iterrows():
        mark = "SLOW" if r["slow_response"] else ""
        print(f"  {r['cast_id']:16s} {r['tau_s']:6.1f} s  {mark}" if pd.notna(r["tau_s"]) else
              f"  {r['cast_id']:16s}      - s  {mark} (no tau from step02b)")
    s = flags.dropna(subset=["tau_s"])
    above, below = s.loc[s["slow_response"], "tau_s"], s.loc[~s["slow_response"], "tau_s"]
    if len(above) and len(below):
        print(f"  gap at the threshold: slowest unflagged {below.max():.1f} s, fastest flagged {above.min():.1f} s")
    try:
        path = write_response(cdir, flags, source, overwrite=a.force)
    except (FileExistsError, ValueError) as e:
        print(e)
        return 1
    print(f"{a.cruise}: {describe(cdir)}\n  wrote {path}")
    print("  next: rerun step08_grid_l3 and step06_hydrographic_sections for this cruise")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

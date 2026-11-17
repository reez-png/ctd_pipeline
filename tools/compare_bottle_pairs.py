"""Set step02b's and the lag diagnosis's bottle pairs side by side, bottle by bottle.

    python tools/compare_bottle_pairs.py P45_07

Both tools pair bottles by one rule (suncat.nitrate.bottle_match) but fit different frames: step02b the L2
co-located product, the diagnosis step03's frames. This reads the per-bottle tables they write
(_audit/suna/<CRUISE>__suna_response_bottles.csv from step02b, __suna_lag_bottles.csv from
tools/diagnose_suna_lag.py) and separates the two reasons their scores can differ:

  composition  which bottles each compared (and why the other left a bottle out)
  frames       on the SAME bottles, the scores of each tool's down-cast and lag-corrected values

then lists the bottles that carry the misfit (largest residuals in bottle units, from each tool's own fit) and
the bottles whose values differ most between the tools. Writes __bottle_pairs_compared.csv beside the inputs.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src"))
import project_paths as pp  # noqa: E402
import suncat.nitrate.bottle_match as BM  # noqa: E402

TOOLS = {"02b": ("__suna_response_bottles.csv", "N@tau"), "diag": ("__suna_lag_bottles.csv", "N@tau_cruise")}


def _read(path: Path, ntau: str) -> pd.DataFrame:
    t = pd.read_csv(path, encoding="utf-8-sig")
    t = t.rename(columns={ntau: "Ntau", "N@0": "N0"})
    t["reason"] = t["reason"].fillna("").astype(str)
    t["cast"] = t["cast"].astype(str)
    t["k"] = t.groupby(["cast", "depth_m", "bottle"]).cumcount()          # replicate rows stay apart
    return t[["cast", "depth_m", "bottle", "k", "down", "N0", "Ntau", "reason"]]


def _resid(t: pd.DataFrame, col: str) -> pd.Series:
    """(sensor - offset) / gain - bottle, from the tool's own OLS over its compared bottles."""
    ok = t["reason"].eq("") & t[col].notna()
    out = pd.Series(np.nan, index=t.index)
    if ok.sum() >= 3:
        g, o = np.polyfit(t.loc[ok, "bottle"], t.loc[ok, col], 1)
        out[ok] = (t.loc[ok, col] - o) / g - t.loc[ok, "bottle"]
    return out


def compare(cruise: str, top: int = 8, say=print) -> dict:
    a = pp.cruise_dir(cruise) / "_audit" / "suna"
    paths = {k: a / f"{cruise}{f}" for k, (f, _) in TOOLS.items()}
    missing = [p.name for p in paths.values() if not p.exists()]
    if missing:
        say(f"missing {missing} in {a}: run step02b and tools/diagnose_suna_lag.py for {cruise} first")
        return {}
    T = {k: _read(paths[k], TOOLS[k][1]) for k in TOOLS}
    for k, t in T.items():
        t["e_down"], t["e_Ntau"] = _resid(t, "down"), _resid(t, "Ntau")
    M = T["02b"].merge(T["diag"], on=["cast", "depth_m", "bottle", "k"], how="outer", suffixes=("_02b", "_diag"))
    for k in TOOLS:
        M[f"reason_{k}"] = M[f"reason_{k}"].fillna("not in this tool's table")
    in_a, in_b = M["reason_02b"].eq(""), M["reason_diag"].eq("")
    say(f"{cruise}: bottle pairs, step02b (L2 co-located frames) against the diagnosis (step03 frames)")
    say(f"  compared by both {int((in_a & in_b).sum())}; only by step02b {int((in_a & ~in_b).sum())}; only by the "
        f"diagnosis {int((~in_a & in_b).sum())}; by neither {int((~in_a & ~in_b).sum())}")
    for lab, m, other in (("only step02b", in_a & ~in_b, "reason_diag"), ("only the diagnosis", in_b & ~in_a, "reason_02b")):
        if m.any():
            say(f"  {lab}: the other tool's reason: " +
                "; ".join(f"{v} {r}" for r, v in M.loc[m, other].value_counts().items()))
    # frames: the same bottles, each tool's values
    both = M[in_a & in_b]
    res = {}
    if len(both) >= 5:
        say(f"\n  on the same {len(both)} bottles (composition removed; what differs is the frames):")
        say(f"  {'tool':6s} {'estimator':9s} {'gain (95% CI)':>27s} {'offset':>8s} {'R2':>6s}  LOOCV RMSE (umol/L)")
        for k in TOOLS:
            tb = pd.DataFrame({"cast": both["cast"], "depth_m": both["depth_m"], "bottle": both["bottle"],
                               "down": both[f"down_{k}"], "N@tau": both[f"Ntau_{k}"], "reason": ""})
            c = BM.compare(tb, ["N@tau"])
            res[k] = c
            for e, f in c["fits"].items():
                say(f"  {k:6s} {e:9s} {f['gain']:7.4f} ({f['lo']:.4f} to {f['hi']:.4f}) {f['offset']:+8.3f} "
                    f"{f['R2']:6.3f}  {f['loocv']:.3f}")
            f = c["fits"]["N@tau"]
            say(f"  {k:6s} N@tau vs down: better on {f['better']} of {f['n']}, sign test p = {f['p_sign']:.3f}")
        dd, dn = (both["down_02b"] - both["down_diag"]).abs(), (both["Ntau_02b"] - both["Ntau_diag"]).abs()
        say(f"  between the tools, same bottle: |down difference| median {dd.median():.2f}, max {dd.max():.2f}; "
            f"|N@tau difference| median {dn.median():.2f}, max {dn.max():.2f} (sensor units)")
    # the bottles that carry the misfit
    for k in TOOLS:
        t = M[M[f"reason_{k}"].eq("")].copy()
        t["abs"] = t[f"e_Ntau_{k}"].abs()
        t = t.sort_values("abs", ascending=False).head(top)
        say(f"\n  largest N@tau residuals in {k} (bottle units; the other tool's on the same bottle beside it):")
        say(f"    {'cast':16s} {'depth':>6s} {'bottle':>7s} {'e_down':>7s} {'e_Ntau':>7s}   other e_down e_Ntau / reason")
        o = "diag" if k == "02b" else "02b"
        for _, r in t.iterrows():
            oth = (f"{r[f'e_down_{o}']:+7.2f} {r[f'e_Ntau_{o}']:+7.2f}" if r[f"reason_{o}"] == ""
                   else r[f"reason_{o}"])
            say(f"    {r['cast']:16s} {r['depth_m']:6.1f} {r['bottle']:7.2f} {r[f'e_down_{k}']:+7.2f} "
                f"{r[f'e_Ntau_{k}']:+7.2f}   {oth}")
    out = a / f"{cruise}__bottle_pairs_compared.csv"
    M.to_csv(out, index=False, encoding="utf-8-sig")
    say(f"\nper-bottle table -> {out}")
    return dict(table=M, same_bottles=res)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cruise")
    ap.add_argument("--top", type=int, default=8, help="bottles listed per tool (default 8)")
    a = ap.parse_args(argv)
    return 0 if compare(a.cruise, a.top) else 1


if __name__ == "__main__":
    raise SystemExit(main())

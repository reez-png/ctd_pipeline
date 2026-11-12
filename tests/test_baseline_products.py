"""Regression test: current pipeline products vs the frozen baseline.

The baseline is created by tools/make_baseline.py on a machine that has the cruise data. It copies the
trusted L2/L3 product CSVs into tests/baseline/local/ and records a manifest. This test reloads each
frozen copy and the matching LIVE product from the cruise folder and asserts they still agree within
tolerance, so a refactor that silently moves a number is caught.

It SKIPS cleanly when:
  - there is no baseline (tests/baseline/local/manifest.json absent): e.g. on CI, or before you freeze;
  - no live product is present: you have not re-run the producing step in this working tree yet.

So the meaningful gate is:  refactor  ->  re-run the producing notebook/script  ->  pytest.

Tolerance: numeric columns must match within rtol=1e-6, atol=1e-8 (byte-faithful refactors reproduce
the CSV text exactly, so this only absorbs harmless float-formatting jitter). Columns, row count, NaN
pattern and every text column must match exactly. All mismatches across all products are reported at
once, so after a refactor you see every affected product, not just the first.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

RTOL = 1e-6
ATOL = 1e-8

_HERE = Path(__file__).resolve().parent
_LOCAL = _HERE / "baseline" / "local"
_MANIFEST = _LOCAL / "manifest.json"


def _cruise_dir(cruise: str) -> Path:
    # Resolve the cruise folder the same way make_baseline.py did (via suncat.config). Works whether
    # suncat is installed editable or only present under src/.
    sys.path.insert(0, str(_HERE.parent / "src"))
    from suncat import config
    return config.for_cruise(cruise).cruise_dir


def _compare(rel: str, frozen: Path, live: Path) -> list[str]:
    """Return a list of human-readable problems for one product (empty if it matches)."""
    problems: list[str] = []
    a = pd.read_csv(frozen)
    b = pd.read_csv(live)

    if list(a.columns) != list(b.columns):
        return [f"{rel}: columns changed\n    frozen={list(a.columns)}\n    live  ={list(b.columns)}"]
    if len(a) != len(b):
        return [f"{rel}: row count changed (frozen={len(a)}, live={len(b)})"]

    for col in a.columns:
        fa, fb = a[col], b[col]
        num = pd.api.types.is_numeric_dtype(fa) and pd.api.types.is_numeric_dtype(fb)
        if num:
            xa = fa.to_numpy(dtype=float)
            xb = fb.to_numpy(dtype=float)
            na, nb = np.isnan(xa), np.isnan(xb)
            if not np.array_equal(na, nb):
                problems.append(f"{rel}:{col}: NaN pattern changed")
                continue
            if not np.allclose(xa[~na], xb[~nb], rtol=RTOL, atol=ATOL):
                d = float(np.max(np.abs(xa[~na] - xb[~nb]))) if (~na).any() else float("nan")
                problems.append(f"{rel}:{col}: values moved (max abs diff {d:.3g}, rtol={RTOL}, atol={ATOL})")
        else:
            if not fa.astype("object").equals(fb.astype("object")):
                problems.append(f"{rel}:{col}: text column changed")
    return problems


def test_products_match_baseline():
    if not _MANIFEST.exists():
        pytest.skip("no regression baseline present (run tools/make_baseline.py where the cruise data lives)")

    manifest = json.loads(_MANIFEST.read_text())
    cruise_dir = _cruise_dir(manifest["cruise"])

    checked: list[str] = []
    skipped: list[str] = []
    problems: list[str] = []

    for f in manifest["files"]:
        rel = f["path"]
        frozen = _LOCAL / rel
        live = cruise_dir / rel

        # Integrity: the frozen reference copy must not have been edited since it was frozen.
        from suncat import baseline
        if frozen.exists() and baseline.file_sha256(frozen) != f["sha256"]:
            problems.append(f"{rel}: frozen baseline copy was modified (sha256 mismatch); re-freeze it")
            continue

        if not live.exists():
            skipped.append(rel)
            continue

        problems.extend(_compare(rel, frozen, live))
        checked.append(rel)

    if problems:
        raise AssertionError("baseline regressions (" + str(len(problems)) + "):\n  " + "\n  ".join(problems))
    if not checked:
        pytest.skip(
            f"no live products to compare (re-run the producing steps first); "
            f"{len(skipped)} frozen products have no live match"
        )

"""Freeze the current pipeline outputs as a regression baseline, BEFORE any refactor.

Run on a machine that has the cruise data:

    python tools/make_baseline.py

It copies the selected product files into tests/baseline/local/ (preserving their folder layout under
the cruise, so files that share a name in different folders do not clash) and records a sha256 of each,
plus any scalar calibration results you list, into manifest.json. Later, after a suncat module is
extracted and the producing step is re-run, tests/test_baseline_products.py reloads each frozen copy
and the freshly produced product and checks they still agree within tolerance.

Nothing here is committed: the whole tests/baseline/local/ folder is git-ignored (add the rule from
the README/handover if it is not already there). The frozen copies are REAL cruise products.

Edit CRUISE, ROOT, PRODUCTS and SCALARS to match your repo, then confirm the printed inventory
at the end lists every product you meant to freeze.
"""
from __future__ import annotations
import json, shutil, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from suncat import baseline, config
import project_paths as pp

# ---- EDIT THIS BLOCK -------------------------------------------------------
CRUISE = "P45_06"
ROOT   = None          # None -> $SUNCAT_CTD_ROOT or CWD; or set r"C:\Projects\ctd_pipeline"

# Finished products to freeze, as (product key, file glob). The key is looked up in the layout registry
# (project_paths.PRODUCTS), so the baseline follows ACTIVE_LAYOUT. If a line prints "(no match)", check
# the DISCOVERED inventory printed at the end and correct the glob here.
PRODUCTS = [
    ("nitrate_qc",           "*_nitrate_QC.csv"),          # step03: per-cast SUNA nitrate QC (primary product)
    ("chl_calibrated",       "*_chl_L3_corrected.csv"),    # step07: chlorophyll calibration (NPQ removed)
    ("chl_calibrated_noNPQ", "*_chl_L3_corrected.csv"),    # step07_noremoval: chlorophyll (NPQ kept)
    ("suna_calibrated",      "*.csv"),                     # bottle-corrected SUNA
    ("suna_l2",              "*_SUNA_1s.csv"),             # step02: merged CTD+SUNA 1 s (mid-chain checkpoint)
    ("suna_l2",              "*_suna_L2_colocated.csv"),   # step02 standalone co-location (cross-check)
]

# Scalar results to freeze (optional but recommended). Fill from your current calibration outputs,
# e.g. the step03 bottle fit and the step07 chlorophyll fit. These are checked by eye in the manifest
# and can be asserted later; they are a human-readable record of "what the numbers were tonight".
#   {"key": "suna_bottle_fit", "reference": [<intercept>, <slope>], "rtol": 1e-6}
SCALARS = [
    # {"key": "suna_bottle_fit", "reference": [-1.097, 0.871], "rtol": 1e-6},   # step03 P45_06
    # {"key": "chl_fit",         "reference": [-0.004, 1.138], "rtol": 1e-6},    # step07 P45_06
]
# ---------------------------------------------------------------------------


def _discover(base: Path) -> None:
    """Print the CSV products actually present, so wrong globs are easy to fix."""
    print("\nDISCOVERED product CSVs under this cruise (for cross-checking the globs above):")
    layout = pp.LAYOUTS[pp.active_layout()]
    for key in sorted({k for k, _ in PRODUCTS}):
        for h in sorted((base / layout[key]).glob("*.csv")):
            print("   ", f"[{key}]", h.relative_to(base).as_posix())


def main() -> int:
    cfg = config.for_cruise(CRUISE, root=ROOT)
    base = cfg.cruise_dir
    if not base.exists():
        print(f"ERROR: cruise folder not found: {base}")
        print("Set ROOT above, or $SUNCAT_CTD_ROOT, to your ctd_pipeline folder.")
        return 1
    out = Path(__file__).resolve().parents[1] / "tests" / "baseline" / "local"
    out.mkdir(parents=True, exist_ok=True)

    files = []
    seen = set()
    layout = pp.LAYOUTS[pp.active_layout()]
    print(f"layout {pp.active_layout()}")
    for prod, pat in PRODUCTS:
        matched = sorted((base / layout[prod]).glob(pat))
        if not matched:
            print(f"  (no match for {prod}: {layout[prod]}/{pat})")
        for src in matched:
            rel = src.relative_to(base)                 # keep the folder structure under the cruise
            key = rel.as_posix()
            if key in seen:                              # a file caught by two globs: freeze once
                continue
            seen.add(key)
            dest = out / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
            files.append({"path": key, "product": prod,
                          "sha256": baseline.file_sha256(dest),
                          "from": str(src)})
            print("froze", key)

    manifest = {"cruise": CRUISE, "layout": pp.active_layout(), "files": files, "products": SCALARS}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))

    _discover(base)
    if not files:
        print("\nWARNING: no product files matched. Fix PRODUCT_GLOBS using the DISCOVERED list above.")
    print(f"\nbaseline: {len(files)} files + {len(SCALARS)} scalars -> {out / 'manifest.json'}")
    print("Reminder: tests/baseline/local/ must be git-ignored. Never commit these copies.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

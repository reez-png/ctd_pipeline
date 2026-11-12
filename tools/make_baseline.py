"""Freeze the current pipeline outputs as a regression baseline, BEFORE any refactor.

Run on a machine that has the cruise data:

    python tools/make_baseline.py

It copies the selected product files into tests/baseline/local/ (preserving their L2/L3/... subfolder
layout, so files that share a name in different folders do not clash) and records a sha256 of each,
plus any scalar calibration results you list, into manifest.json. Later, after a suncat module is
extracted and the producing step is re-run, tests/test_baseline_products.py reloads each frozen copy
and the freshly produced product and checks they still agree within tolerance.

Nothing here is committed: the whole tests/baseline/local/ folder is git-ignored (add the rule from
the README/handover if it is not already there). The frozen copies are REAL cruise products.

Edit CRUISE, ROOT, PRODUCT_GLOBS and SCALARS to match your repo, then confirm the printed inventory
at the end lists every product you meant to freeze.
"""
from __future__ import annotations
import json, shutil, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from suncat import baseline, config

# ---- EDIT THIS BLOCK -------------------------------------------------------
CRUISE = "P45_06"
ROOT   = None          # None -> $SUNCAT_CTD_ROOT or CWD; or set r"C:\Projects\ctd_pipeline"

# Finished products to freeze, as globs RELATIVE TO the cruise folder (cruises/<CRUISE>/).
# The names below follow README section 4b. If a glob prints "(no match)", check the DISCOVERED
# inventory printed at the end and correct the name here, these are the only things you must get right.
PRODUCT_GLOBS = [
    "L3/*_nitrate_QC.csv",                       # step03: per-cast SUNA nitrate QC (primary product)
    "L3/productivity/*_productivity.csv",        # step07: chlorophyll calibration (NPQ removed)
    "L3/productivity_noNPQ/*_productivity.csv",  # step07_noremoval: chlorophyll (NPQ kept)
    "L3/SUNA/*.csv",                             # bottle-corrected SUNA L3, if your step writes one
    "L2/SUNA/*_SUNA_1s.csv",                     # step02: merged CTD+SUNA 1 s (mid-chain checkpoint)
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
    for pat in ("L3/*.csv", "L3/**/*.csv", "L2/SUNA/*.csv"):
        hits = sorted(base.glob(pat))
        for h in hits:
            print("   ", h.relative_to(base).as_posix())


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
    for pat in PRODUCT_GLOBS:
        matched = sorted(base.glob(pat))
        if not matched:
            print(f"  (no match for {pat})")
        for src in matched:
            rel = src.relative_to(base)                 # keep the L3/<sub>/<file> structure
            key = rel.as_posix()
            if key in seen:                              # a file caught by two globs: freeze once
                continue
            seen.add(key)
            dest = out / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
            files.append({"path": key,
                          "sha256": baseline.file_sha256(dest),
                          "from": str(src)})
            print("froze", key)

    manifest = {"cruise": CRUISE, "files": files, "products": SCALARS}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))

    _discover(base)
    if not files:
        print("\nWARNING: no product files matched. Fix PRODUCT_GLOBS using the DISCOVERED list above.")
    print(f"\nbaseline: {len(files)} files + {len(SCALARS)} scalars -> {out / 'manifest.json'}")
    print("Reminder: tests/baseline/local/ must be git-ignored. Never commit these copies.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

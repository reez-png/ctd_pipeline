"""Freeze the current pipeline outputs as a regression baseline, BEFORE any refactor.

Run on a machine that has the cruise data:

    python tools/make_baseline.py

It copies the selected product files into tests/baseline/local/ (preserving their L3/... subfolder
layout, so files that share a name in different folders do not clash) and records a sha256 of each,
plus any scalar calibration results you list, into manifest.json. Later, as suncat modules are
extracted, the regression test regenerates these and checks they still match within tolerance.

Edit CRUISE, ROOT, PRODUCT_GLOBS and SCALARS to match your repo. Nothing here is committed: the
local baseline folder is git-ignored.
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
# These match the real P45_06 layout. Add or remove lines as your products change.
PRODUCT_GLOBS = [
    "L3/*_nitrate_QC.csv",                        # per-cast SUNA nitrate QC
    "L3/SUNA/*_suna_L3_corrected.csv",            # calibrated SUNA L3
    "L3/productivity/*_chl_L3_corrected.csv",     # chlorophyll L3 (NPQ corrected)
    "L3/productivity_noNPQ/*_chl_L3_corrected.csv",  # chlorophyll L3 (no NPQ)
]

# Scalar results to freeze (optional but recommended). Fill from your current calibration outputs.
#   {"key": "suna_bottle_fit", "reference": [<intercept>, <slope>], "atol": 0.0, "rtol": 1e-6}
SCALARS = [
    # {"key": "suna_bottle_fit", "reference": [-0.39, 0.89], "rtol": 1e-6},
    # {"key": "chl_fit",         "reference": [0.457, 0.193], "rtol": 1e-6},
]
# ---------------------------------------------------------------------------


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
    for pat in PRODUCT_GLOBS:
        matched = sorted(base.glob(pat))
        if not matched:
            print(f"  (no match for {pat})")
        for src in matched:
            rel = src.relative_to(base)                 # keep the L3/<sub>/<file> structure
            dest = out / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
            files.append({"path": str(rel).replace("\\", "/"),
                          "sha256": baseline.file_sha256(dest),
                          "from": str(src)})
            print("froze", rel)

    if not files:
        print("\nWARNING: no product files matched. Check PRODUCT_GLOBS and ROOT.")
    manifest = {"cruise": CRUISE, "files": files, "products": SCALARS}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nbaseline: {len(files)} files + {len(SCALARS)} scalars -> {out / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

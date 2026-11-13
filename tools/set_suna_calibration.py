"""Write (or show) a cruise's SUNA Stage 3 calibration file, cruises/<CRUISE>/metadata/suna_calibration.json.

    python tools/set_suna_calibration.py P45_06 fit
    python tools/set_suna_calibration.py P45_07 adopted --a1 0.89 --a0 -0.39 --source "homogeneous-layer fit ..."
    python tools/set_suna_calibration.py P45_07 --show

step03 reads this file and stops without it, so each cruise states its calibration once, here, and no run
command or notebook default can change it. An existing file is not replaced unless --force is given.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "src"))
import project_paths as pp  # noqa: E402
from suncat.suna_cal import (SunaCalibrationError, calibration_path, load_suna_calibration,  # noqa: E402
                             write_suna_calibration)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cruise")
    ap.add_argument("mode", nargs="?", choices=["fit", "adopted"])
    ap.add_argument("--a1", type=float, help="gain (adopted only)")
    ap.add_argument("--a0", type=float, help="offset, umol/L (adopted only)")
    ap.add_argument("--source", help="where a1 and a0 come from (required for adopted)")
    ap.add_argument("--force", action="store_true", help="replace an existing file")
    ap.add_argument("--show", action="store_true", help="print the current file and stop")
    a = ap.parse_args()
    cdir = pp.cruise_dir(a.cruise)
    if a.show or a.mode is None:
        try:
            rec = load_suna_calibration(cdir)
        except (FileNotFoundError, SunaCalibrationError) as e:
            print(e)
            return 1
        print(f"{a.cruise}: {rec['describe']}\n  {rec['path']}")
        return 0
    try:
        path = write_suna_calibration(cdir, a.mode, a.a1, a.a0, a.source, overwrite=a.force)
    except (FileExistsError, SunaCalibrationError) as e:
        print(e)
        return 1
    print(f"{a.cruise}: {load_suna_calibration(cdir)['describe']}\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

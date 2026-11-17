"""tools/compare_bottle_pairs.py: step02b's and the diagnosis's per-bottle tables side by side."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("compare_bottle_pairs", REPO / "tools" / "compare_bottle_pairs.py")
CBP = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CBP)


def test_composition_and_frames_are_separated(tmp_path, monkeypatch):
    monkeypatch.setenv("CTD_PIPELINE_ROOT", str(tmp_path))
    a = tmp_path / "cruises" / "ZZ" / "_audit" / "suna"
    a.mkdir(parents=True)
    rng = np.random.default_rng(0)
    b = np.linspace(0, 20, 30)
    base = pd.DataFrame({"cast": [f"C_{i % 3:02d}" for i in range(30)], "depth_m": np.arange(30.0), "bottle": b})
    t1 = base.assign(down=0.9 * b + rng.normal(0, 0.5, 30), **{"N@0": 0.9 * b, "N@tau": 0.9 * b + 0.1}, reason="")
    t2 = base.assign(down=0.9 * b + rng.normal(0, 0.5, 30), **{"N@0": 0.9 * b, "N@tau_cruise": 0.9 * b}, reason="")
    t2.loc[0, "reason"] = "no moving down-cast frame within +-2 m"
    t2.loc[1, "N@tau_cruise"] = 0.9 * b[1] + 15.0                     # one bottle the diagnosis gets badly wrong
    t1.to_csv(a / "ZZ__suna_response_bottles.csv", index=False)
    t2.to_csv(a / "ZZ__suna_lag_bottles.csv", index=False)
    lines = []
    r = CBP.compare("ZZ", say=lines.append)
    text = "\n".join(lines)
    assert "compared by both 29; only by step02b 1; only by the diagnosis 0" in text
    assert "1 no moving down-cast frame within +-2 m" in text
    s = r["same_bottles"]
    assert s["02b"]["fits"]["N@tau"]["R2"] > 0.99 > s["diag"]["fits"]["N@tau"]["R2"]
    # the bad bottle heads the diagnosis's residual list
    i = next(k for k, ln in enumerate(lines) if "largest N@tau residuals in diag" in ln)
    assert lines[i + 2].split()[0] == "C_01" and float(lines[i + 2].split()[1]) == 1.0
    assert (a / "ZZ__bottle_pairs_compared.csv").exists()

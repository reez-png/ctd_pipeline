import json
from pathlib import Path
import pytest
from suncat import baseline

BASE = Path(__file__).parent / "baseline" / "local"
MANIFEST = BASE / "manifest.json"
pytestmark = pytest.mark.baseline


def _load():
    if not MANIFEST.exists():
        pytest.skip("no local baseline captured; run tools/make_baseline.py on a machine with the data")
    return json.loads(MANIFEST.read_text())


def test_reference_files_intact():
    """The frozen reference outputs are present and unchanged on disk."""
    m = _load()
    missing, changed = [], []
    for e in m.get("files", []):
        p = BASE / e["path"]
        if not p.exists():
            missing.append(e["path"]); continue
        if baseline.file_sha256(p) != e["sha256"]:
            changed.append(e["path"])
    assert not missing, f"baseline files missing: {missing}"
    assert not changed, f"baseline files changed on disk since capture: {changed}"


def test_producers_reproduce_reference():
    """Once a suncat module registers a producer, its regenerated output must match the frozen
    reference within tolerance. Keys with no producer yet are skipped (module not extracted)."""
    m = _load()
    checked = 0
    for e in m.get("products", []):
        fn = baseline.producer(e["key"])
        if fn is None:
            continue
        assert baseline.within_tol(e["reference"], fn(), e.get("atol", 0.0), e.get("rtol", 1e-6)), \
            f"producer '{e['key']}' no longer reproduces the frozen reference"
        checked += 1
    if checked == 0:
        pytest.skip("no suncat producers registered yet; reference is frozen and intact")

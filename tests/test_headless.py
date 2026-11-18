"""Headless notebook runs: tools/run_notebook.py marks the kernel SUNCAT_HEADLESS=1 and every notebook that builds
interactive widgets skips them then (no frontend; a step06 cell after its widget viewer hung twice on Windows)."""
from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def test_every_widget_notebook_checks_for_a_headless_run():
    missing = []
    for nb in sorted((REPO / "processing_scripts" / "notebook").rglob("*.ipynb")):
        if ".ipynb_checkpoints" in nb.parts:      # Jupyter's autosave copies: not the pipeline's notebooks
            continue
        src = "\n".join("".join(c["source"]) for c in json.loads(nb.read_text(encoding="utf-8"))["cells"]
                        if c["cell_type"] == "code")
        if "import ipywidgets" in src and "SUNCAT_HEADLESS" not in src:
            missing.append(nb.name)
    assert not missing, missing


def test_the_runner_sets_it():
    assert 'os.environ["SUNCAT_HEADLESS"] = "1"' in (REPO / "tools" / "run_notebook.py").read_text(encoding="utf-8")

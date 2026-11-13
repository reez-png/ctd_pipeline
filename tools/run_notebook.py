"""Run a pipeline notebook from the command line, optionally at another git revision, with settings overridden.

    python tools/run_notebook.py processing_scripts/notebook/step03_nitrate_qc.ipynb ^
        --set CRUISE_ID=P45_06 --set ADOPT_CALIBRATION=False
    python tools/run_notebook.py processing_scripts/notebook/step03_nitrate_qc.ipynb --rev main --set CRUISE_ID=P45_06

Why: an A/B comparison (or re-making the baseline) must run exactly the code in git, not whatever an editor
has open. This reads the notebook from the working tree, or from `git show <rev>:<path>` with --rev, replaces
each `NAME = ...` settings line given with --set (the line must occur exactly once in the code cells),
and executes it with the notebook's own folder as the working directory. The notebook file is never
written, so no outputs or edited settings end up in the repo. --save writes the executed copy elsewhere.
"""
from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path

import nbformat


def _value(text: str) -> str:
    """Python source for a --set value: a literal stays as written (False, 2.0, 'x'); anything else is a string."""
    try:
        ast.literal_eval(text)
        return text
    except (ValueError, SyntaxError):
        return repr(text)


def load(path: Path, rev: str | None):
    if rev is None:
        return nbformat.read(path, as_version=4)
    rel = path.resolve().relative_to(Path.cwd().resolve()).as_posix()
    raw = subprocess.run(["git", "show", f"{rev}:{rel}"], capture_output=True, check=True).stdout
    return nbformat.reads(raw.decode("utf-8"), as_version=4)


def apply_settings(nb, settings: list[str]) -> list[str]:
    done = []
    for item in settings:
        name, _, val = item.partition("=")
        name = name.strip()
        pat = re.compile(rf"^{re.escape(name)}\s*=.*$", re.M)
        cells = [c for c in nb.cells if c.cell_type == "code" and pat.search(c.source)]
        n = sum(len(pat.findall(c.source)) for c in cells)
        if n != 1:
            raise SystemExit(f"--set {name}: expected one '{name} = ...' line in the code cells, found {n}")
        line = f"{name} = {_value(val.strip())}   # set by run_notebook"
        cells[0].source = pat.sub(lambda _m: line, cells[0].source, count=1)
        done.append(line)
    return done


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("notebook", type=Path)
    ap.add_argument("--rev", help="git revision to take the notebook from (default: working tree)")
    ap.add_argument("--set", action="append", default=[], metavar="NAME=VALUE")
    ap.add_argument("--timeout", type=int, default=3600, help="seconds per cell")
    ap.add_argument("--save", type=Path, help="write the executed notebook here (keep it outside the repo)")
    ap.add_argument("--quiet", action="store_true", help="do not print cell output")
    a = ap.parse_args(argv)

    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError

    nb = load(a.notebook, a.rev)
    for line in apply_settings(nb, a.set):
        print("  " + line)
    where = a.notebook.resolve().parent
    print(f"running {a.notebook.name} ({'git ' + a.rev if a.rev else 'working tree'}) in {where}")
    client = NotebookClient(nb, timeout=a.timeout, kernel_name="python3",
                            resources={"metadata": {"path": str(where)}})
    status = 0
    try:
        client.execute()
    except CellExecutionError as exc:
        print(exc)
        status = 1
    if not a.quiet:
        for c in nb.cells:
            for o in c.get("outputs", []):
                if o.get("output_type") == "stream":
                    sys.stdout.write(o.get("text", ""))
                elif o.get("output_type") in ("execute_result", "display_data"):
                    txt = o.get("data", {}).get("text/plain")
                    if txt:
                        print(txt)
    if a.save:
        nbformat.write(nb, a.save)
        print("executed copy ->", a.save)
    print("FAILED" if status else "OK")
    return status


if __name__ == "__main__":
    raise SystemExit(main())

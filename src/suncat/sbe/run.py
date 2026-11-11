"""Thin runner for SBEDataProcessing modules (the optional, Windows-only backend).

Consolidates the per-notebook run_sbe / run_filter: build the SBE command, run it, and return the
output .cnv path (or None). Mirrors the step01 command pattern exactly:
    <exe> /p<psa> /i<in_cnv> /o<work_dir> /f<out_stem> /s
The portable core never imports this module; it is used only where the SBE software is present.
"""
from __future__ import annotations
import subprocess
from pathlib import Path


def run_module(exe, psa, in_cnv, out_stem, *, sbe_bin, work_dir, echo=True, warn=True):
    """Run one SBEDataProcessing module (FilterW.exe, LoopEditW.exe, ...) on a .cnv.

    Returns work_dir/<out_stem>.cnv if it was written, else None. echo prints the command line; warn
    prints a non-zero return code. Both default on (the step00b behaviour); pass echo=False, warn=False
    for the quiet step00a behaviour.
    """
    cmd = [str(Path(sbe_bin) / exe), f"/p{psa}", f"/i{in_cnv}", f"/o{work_dir}", f"/f{out_stem}", "/s"]
    if echo:
        print("   ", " ".join(str(x) for x in cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if warn and res.returncode != 0:
        print("    WARNING rc", res.returncode, (res.stderr or "")[:200])
    out = Path(work_dir) / f"{out_stem}.cnv"
    return out if out.exists() else None

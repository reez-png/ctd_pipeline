"""Equivalence + contract tests for suncat.sbe.command.build_command.

The _reference function below is the exact pure logic extracted from step01's build_sbe_command. The
fuzz test pins build_command to it across many random inputs, so any future drift is caught; the named
tests document the command contract (argument order, /f stem, optional /c, local-name rules, /s last).

Paths are POSIX-style so .resolve() behaves the same on Linux CI and on Windows.
"""
from __future__ import annotations
import random
from pathlib import Path
from typing import List, Optional

from suncat.sbe.command import build_command, command_to_text


def _reference(
    exe_path: Path,
    input_file: Path,
    output_dir: Path,
    output_file: str,
    psa_file: Path,
    config_file: Optional[Path] = None,
    extra_args: Optional[List[str]] = None,
    use_local_names: bool = False,
) -> List[str]:
    extra_args = extra_args or []
    input_arg = input_file.name if use_local_names else str(input_file)
    config_arg = None
    if config_file is not None:
        if use_local_names and config_file.parent.resolve() == input_file.parent.resolve():
            config_arg = config_file.name
        else:
            config_arg = str(config_file)
    stem = Path(output_file).stem
    cmd = [str(exe_path), f"/i{input_arg}", f"/p{psa_file}", f"/f{stem}", f"/o{output_dir}"]
    if config_arg is not None:
        cmd.append(f"/c{config_arg}")
    cmd.extend(extra_args)
    cmd.append("/s")
    return cmd


# --- contract tests ---------------------------------------------------------

EXE = Path("/sbe/DatCnvW.exe")
IN = Path("/work/raw/P45_06_CTD_03.hex")
OUT = Path("/work/out")
PSA = Path("/psa/01_datcnv.psa")


def test_order_and_stem_no_config():
    cmd = build_command(EXE, IN, OUT, "P45_06_CTD_03.cnv", PSA)
    assert cmd == [
        str(EXE),
        f"/i{IN}",            # full path when use_local_names is False (default)
        f"/p{PSA}",
        "/fP45_06_CTD_03",    # /f is the stem, extension dropped
        f"/o{OUT}",
        "/s",
    ]


def test_no_config_has_no_c_flag():
    cmd = build_command(EXE, IN, OUT, "x.cnv", PSA)
    assert not any(a.startswith("/c") for a in cmd)
    assert cmd[-1] == "/s"


def test_local_names_use_bare_input():
    cmd = build_command(EXE, IN, OUT, "x.cnv", PSA, use_local_names=True)
    assert cmd[1] == "/iP45_06_CTD_03.hex"


def test_config_local_only_when_beside_input():
    beside = IN.with_name("19-8460.xmlcon")
    cmd = build_command(EXE, IN, OUT, "x.cnv", PSA, config_file=beside, use_local_names=True)
    assert "/c19-8460.xmlcon" in cmd            # bare name, because it is beside the input

    elsewhere = Path("/cal/19-8460.xmlcon")
    cmd2 = build_command(EXE, IN, OUT, "x.cnv", PSA, config_file=elsewhere, use_local_names=True)
    assert f"/c{elsewhere}" in cmd2             # full path, because it is not beside the input


def test_config_full_path_when_not_local():
    beside = IN.with_name("19-8460.xmlcon")
    cmd = build_command(EXE, IN, OUT, "x.cnv", PSA, config_file=beside, use_local_names=False)
    assert f"/c{beside}" in cmd                 # use_local_names False -> always full path


def test_extra_args_before_terminal_s():
    cmd = build_command(EXE, IN, OUT, "x.cnv", PSA, extra_args=["/a1", "/a2"])
    assert cmd[-3:] == ["/a1", "/a2", "/s"]


def test_command_to_text_roundtrips():
    cmd = build_command(EXE, IN, OUT, "x.cnv", PSA)
    text = command_to_text(cmd)
    assert isinstance(text, str) and "/s" in text


# --- fuzz equivalence probe -------------------------------------------------

def test_matches_reference_over_random_inputs():
    rng = random.Random(20261112)
    exes = [Path("/sbe/DatCnvW.exe"), Path("/sbe/FilterW.exe"), Path("/sbe/LoopEditW.exe")]
    stems = ["P45_06_CTD_03", "cast01", "a.b", "x"]
    exts = [".cnv", ".asc", ""]
    for _ in range(2000):
        indir = Path(random.choice(["/work/raw", "/data/in", "/x"]))
        inp = indir / (random.choice(stems) + ".hex")
        out = Path(random.choice(["/work/out", "/o"]))
        psa = Path("/psa") / (random.choice(["01_datcnv", "05_filter"]) + ".psa")
        outfile = random.choice(stems) + random.choice(exts)
        use_local = rng.random() < 0.5

        if rng.random() < 0.5:
            cfg = None
        elif rng.random() < 0.5:
            cfg = inp.with_name("19-8460.xmlcon")          # beside the input
        else:
            cfg = Path("/cal/19-8460.xmlcon")              # elsewhere

        extra = None
        if rng.random() < 0.3:
            extra = random.choice([["/a"], ["/a", "/b"], []])

        got = build_command(exes[0], inp, out, outfile, psa,
                            config_file=cfg, extra_args=extra, use_local_names=use_local)
        ref = _reference(exes[0], inp, out, outfile, psa,
                        config_file=cfg, extra_args=extra, use_local_names=use_local)
        assert got == ref

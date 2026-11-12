"""Build the Sea-Bird command line for one processing module (optional Windows backend).

This is step01's production command builder, lifted verbatim: it supports a config file (/c), extra
arguments, and the "local name" convention (pass a bare filename when the runner works from the input
file's own folder). It is PURE: it only assembles the argument list, with no filesystem access and no
module globals, so it is fully testable without the SBE software. The side effects step01's
build_sbe_command also performed (creating the output folder, rejecting paths that contain spaces,
reading the RUN_FROM_INPUT_FILE_FOLDER / USE_LOCAL_INPUT_AND_CONFIG_NAMES globals) stay in step01's
thin wrapper, which computes use_local_names and calls this.

Argument order matches step01 exactly:

    <exe>  /i<input>  /p<psa>  /f<output_stem>  /o<output_dir>  [/c<config>]  [extra...]  /s

This is a different command from suncat.sbe.run.run_module (which is /p /i /o /f /s, no config); the
two are intentionally kept separate.
"""
from __future__ import annotations
import subprocess
from pathlib import Path
from typing import List, Optional


def build_command(
    exe_path: Path,
    input_file: Path,
    output_dir: Path,
    output_file: str,
    psa_file: Path,
    config_file: Optional[Path] = None,
    extra_args: Optional[List[str]] = None,
    *,
    use_local_names: bool = False,
) -> List[str]:
    """Assemble a Sea-Bird no-space command line. Pure: no I/O, no globals.

    use_local_names: when True, pass the input as its bare filename, and the config as its bare
    filename ONLY when it sits beside the input file. This reproduces step01's
    RUN_FROM_INPUT_FILE_FOLDER and USE_LOCAL_INPUT_AND_CONFIG_NAMES behaviour; the caller computes the
    combined flag and passes it in. Default False matches step01's globals().get(..., False) default.

    /f is given the output BASENAME (no extension) so modules do not produce doubled extensions.
    """
    extra_args = extra_args or []

    input_arg = input_file.name if use_local_names else str(input_file)
    config_arg: Optional[str] = None
    if config_file is not None:
        if use_local_names and config_file.parent.resolve() == input_file.parent.resolve():
            config_arg = config_file.name
        else:
            config_arg = str(config_file)

    output_stem = Path(output_file).stem
    command = [
        str(exe_path),
        f"/i{input_arg}",
        f"/p{psa_file}",
        f"/f{output_stem}",
        f"/o{output_dir}",
    ]
    if config_arg is not None:
        command.append(f"/c{config_arg}")
    command.extend(extra_args)
    command.append("/s")
    return command


def command_to_text(command: List[str]) -> str:
    """Render a command list as readable Windows command text."""
    return subprocess.list2cmdline(command)

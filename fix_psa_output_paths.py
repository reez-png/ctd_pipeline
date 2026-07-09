"""
fix_psa_output_paths.py

One-time fixer for the six v2 PSA files. It neutralizes the three fields that
were frozen into File Setup when the PSAs were saved in the SBE GUI:

    <OutputDir value="...">   ->  value=""   (let the /o flag decide)
    <OutputFile value="...">  ->  value=""   (let the /f flag decide)
    <NameAppend value="...">  ->  value=""   (no appended suffix; the runner names files)

Why: step01 already passes /o<dir>, /f<basename> and /s on the command line. When
the PSA also carries hardcoded output settings, Sea-Bird can prefer the PSA values
over the CLI flags, writing output to the wrong place / wrong name (return_code=0
but missing or mislocated output). Blanking these fields makes Sea-Bird fall back
to the command-line flags, which is what we want.

It also rewrites the stale InputDir (some PSAs still point at the OLD ctd_pipeline
project) to blank, since /i supplies the input on the command line too.

A timestamped .bak copy of each PSA is kept next to the original.

USAGE (from the PyCharm terminal, with the v2 venv active):
    python fix_psa_output_paths.py
Then re-run step01 (still RUN_SBE_COMMANDS = False) and confirm the PSA warnings
are gone.
"""

import re
import shutil
from datetime import datetime
from pathlib import Path

# Folder holding the six PSAs. Edit if yours differ.
PSA_DIR = Path(r"C:\Users\OA_2023-03\Projects\ctd_pipeline_v2\ctd\psa")

PSA_FILES = [
    "01_datcnv.psa",
    "02_alignctd.psa",
    "06_loopedit.psa",
    "07_derive.psa",
    "08b_binavg_1m_down.psa",
    "08d_binavg_1s_full.psa",
]

# Tags whose value attribute should be blanked.
BLANK_TAGS = ["OutputDir", "OutputFile", "NameAppend", "InputDir"]


def blank_tag(text: str, tag: str) -> tuple[str, str | None]:
    """Blank value="..." for a <Tag value="..." /> element. Returns (new_text, old_value)."""
    pattern = re.compile(rf'(<{tag}\s+value=")([^"]*)("\s*/?>)')
    m = pattern.search(text)
    if not m:
        return text, None
    old = m.group(2)
    new_text = pattern.sub(rf'\g<1>\g<3>', text, count=1)
    return new_text, old


def fix_one(path: Path) -> dict:
    raw = path.read_text(encoding="latin-1")
    changes = {}
    text = raw
    for tag in BLANK_TAGS:
        text, old = blank_tag(text, tag)
        if old is not None and old != "":
            changes[tag] = old
    if text != raw:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = path.with_suffix(path.suffix + f".bak_{stamp}")
        shutil.copy2(path, backup)
        path.write_text(text, encoding="latin-1")
        changes["_backup"] = backup.name
    return changes


def main():
    if not PSA_DIR.exists():
        raise SystemExit(f"PSA folder not found: {PSA_DIR}")
    print(f"Fixing PSAs in: {PSA_DIR}\n")
    for name in PSA_FILES:
        p = PSA_DIR / name
        if not p.exists():
            print(f"  SKIP (missing): {name}")
            continue
        changes = fix_one(p)
        if changes:
            cleared = {k: v for k, v in changes.items() if k != "_backup"}
            print(f"  {name}")
            for tag, old in cleared.items():
                print(f"      blanked {tag}  (was: {old})")
            print(f"      backup: {changes.get('_backup')}")
        else:
            print(f"  {name}: nothing to change (already clean)")
    print("\nDone. Re-run step01 (RUN_SBE_COMMANDS = False) and check the PSA diagnostic.")


if __name__ == "__main__":
    main()

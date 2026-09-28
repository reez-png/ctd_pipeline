"""
cnv_to_csv.py  -  Convert Sea-Bird .cnv files in a folder to .csv

Standalone utility. Run it on its own:
  - Double-click it, or run:  python cnv_to_csv.py
    -> a folder picker opens; choose a folder that contains .cnv files.
  - Or pass the folder on the command line:  python cnv_to_csv.py "C:\\path\\to\\folder"

What it does:
  - Finds .cnv files in the chosen folder (set RECURSIVE = True to include subfolders).
  - Writes a .csv next to each .cnv with the same name (or into a "csv" subfolder;
    see WRITE_BESIDE below).
  - Column headers come from the .cnv "# name N = ..." lines.
  - Values are copied through exactly as text (no rounding), so the numbers match the .cnv.

Uses only the Python standard library - no pandas or numpy needed.
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

# ============================ options =====================================
RECURSIVE       = False   # True also converts .cnv inside subfolders of the chosen folder
WRITE_BESIDE    = True    # True: write each .csv next to its .cnv. False: into a "csv" subfolder
OUTPUT_SUBDIR   = "csv"   # folder name used when WRITE_BESIDE is False
OVERWRITE       = True    # overwrite an existing .csv of the same name
USE_SHORT_NAMES = True    # header = short name (e.g. prdM, tv290C). False = full label
# ==========================================================================

# Matches:  # name 1 = prDM: Pressure, Strain Gauge [db]
NAME_RE = re.compile(r"#\s*name\s+(\d+)\s*=\s*([^:]+?)\s*:\s*(.*?)\s*$", re.I)


def parse_cnv(path: Path):
    """Return (column_names, data_rows) from an SBE .cnv, or (None, None) if not readable."""
    lines = Path(path).read_text(encoding="latin-1", errors="replace").splitlines()
    names = {}
    end = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("#") or s.startswith("*"):
            m = NAME_RE.search(ln)
            if m:
                idx = int(m.group(1))
                short, full = m.group(2).strip(), m.group(3).strip()
                names[idx] = short if USE_SHORT_NAMES else (full or short)
            if s.upper().startswith("*END*"):
                end = i + 1
                break
    if end is None:
        return None, None
    cols = [names[k] for k in sorted(names)]
    rows = []
    for ln in lines[end:]:
        if ln.strip() == "":
            continue
        parts = ln.split()
        if cols:
            if len(parts) >= len(cols):
                rows.append(parts[:len(cols)])
        else:
            rows.append(parts)
    if not cols and rows:
        cols = [f"col{i + 1}" for i in range(len(rows[0]))]
    return cols, rows


def convert_one(cnv_path: Path, out_dir: Path | None) -> str:
    cols, rows = parse_cnv(cnv_path)
    if cols is None:
        return "skipped (no *END* marker; not a readable .cnv)"
    if out_dir is None:
        out_path = cnv_path.with_suffix(".csv")
    else:
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / (cnv_path.stem + ".csv")
    if out_path.exists() and not OVERWRITE:
        return f"exists, skipped ({out_path.name})"
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        w.writerows(rows)
    return f"-> {out_path.name}  ({len(rows)} rows, {len(cols)} columns)"


def pick_folder() -> Path | None:
    """Command-line argument if given, otherwise a GUI folder picker."""
    if len(sys.argv) > 1:
        return Path(sys.argv[1])
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        chosen = filedialog.askdirectory(title="Select a folder containing .cnv files")
        root.destroy()
        return Path(chosen) if chosen else None
    except Exception:
        print("No folder given and no folder picker available.")
        print("Run again with a path, e.g.:  python cnv_to_csv.py \"C:\\path\\to\\folder\"")
        return None


def main() -> int:
    folder = pick_folder()
    if not folder or not folder.exists() or not folder.is_dir():
        print("No valid folder selected.")
        return 1

    walker = folder.rglob("*") if RECURSIVE else folder.glob("*")
    cnvs = sorted(p for p in walker if p.is_file() and p.suffix.lower() == ".cnv")
    if not cnvs:
        where = "folder or its subfolders" if RECURSIVE else "folder"
        print(f"No .cnv files found in the selected {where}:\n  {folder}")
        return 1

    out_dir = None if WRITE_BESIDE else (folder / OUTPUT_SUBDIR)
    print(f"Converting {len(cnvs)} .cnv file(s) in: {folder}")
    if out_dir:
        print(f"Writing .csv files into: {out_dir}")
    print("-" * 60)

    converted = skipped = errors = 0
    for p in cnvs:
        try:
            msg = convert_one(p, out_dir)
            print(f"  {p.name}: {msg}")
            if msg.startswith("->"):
                converted += 1
            else:
                skipped += 1
        except Exception as exc:
            print(f"  {p.name}: ERROR {exc!r}")
            errors += 1

    print("-" * 60)
    print(f"Done. converted={converted}  skipped={skipped}  errors={errors}")
    return 0


if __name__ == "__main__":
    code = main()
    try:
        input("\nPress Enter to close...")   # keeps the window open when double-clicked
    except (EOFError, KeyboardInterrupt):
        pass
    sys.exit(code)

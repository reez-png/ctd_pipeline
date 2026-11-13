"""Move a cruise's files from one data-level layout to the other (phase 3, docs/DATA_LEVELS.md).

    python tools/migrate_layout.py                          dry run: every cruise, to v2; nothing moves
    python tools/migrate_layout.py --apply                  move every cruise to v2
    python tools/migrate_layout.py --cruise P45_07 --apply  one cruise
    python tools/migrate_layout.py --to v1 --apply          roll back

Each file is classified by the layout registry in project_paths (folder + migrate glob, longest folder
wins), so the rule for what belongs to which product lives in one place. For every file whose product
changes folder:

  1. the whole cruise is planned first; if any destination already holds DIFFERENT content, or two files
     would land on one destination, the cruise is refused and nothing in it moves;
  2. the file is copied to a temporary name beside its destination, the copy's SHA-256 is checked against
     the source, the copy is renamed into place, and only then is the source deleted;
  3. a destination that already holds the SAME content (an interrupted earlier run) just has its source
     removed, so a run can always be repeated;
  4. every move is appended to cruises/<id>/_audit/layout_migration_<UTC time>.csv as it happens.

Files no product claims are never touched; they are listed so you can decide. When a cruise finishes
cleanly, empty source folders are removed and cruises/<id>/_audit/layout.json records the layout, which
project_paths.product_dir() then checks before handing out any folder.

The frozen regression baseline (tests/baseline/local/) is remapped the same way when its cruise moves:
copies move with the same SHA-256 check (against the hashes in the manifest), the paths in manifest.json
are rewritten, and the old manifest is kept as manifest.<from>.json.
"""
from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import project_paths as pp  # noqa: E402

LOG_FIELDS = ["product", "source", "destination", "bytes", "sha256", "action"]


# ------------------------------------------------------------------ classification
def classify(rel: str, layout: str) -> str | None:
    """Product key that owns a cruise-relative posix path in a layout, or None. Longest folder wins."""
    best, best_len = None, -1
    for key, folder in pp.LAYOUTS[layout].items():
        if not rel.startswith(folder + "/"):
            continue
        rest = rel[len(folder) + 1:]
        glob = pp.MIGRATE_GLOB[key]
        if glob != "*" and ("/" in rest or not fnmatch.fnmatch(rest, glob)):
            continue
        if len(folder) > best_len:
            best, best_len = key, len(folder)
    return best


def relocate(rel: str, src: str, dst: str) -> tuple[str | None, str | None]:
    """(product, new relative path) for a path moving from layout src to dst; new path None if it stays."""
    key = classify(rel, src)
    if key is None:
        return None, None
    a, b = pp.LAYOUTS[src][key], pp.LAYOUTS[dst][key]
    return key, (None if a == b else b + rel[len(a):])


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def level_files(cdir: Path):
    """Every file under the top-level folders either layout uses (L0..L3, _figures)."""
    tops = sorted({f.split("/")[0] for lay in pp.LAYOUTS.values() for f in lay.values()})
    for top in tops:
        root = cdir / top
        if root.is_dir():
            for p in sorted(root.rglob("*")):
                if p.is_file():
                    yield p.relative_to(cdir).as_posix()


# ------------------------------------------------------------------ planning
def plan(base: Path, rels, src: str, dst: str):
    """Moves for files under base: list of (key, rel, new_rel, kind). kind is "move" or "dedupe".
    Returns (moves, unclaimed, problems)."""
    moves, unclaimed, problems, targets = [], [], [], {}
    for rel in rels:
        key, new = relocate(rel, src, dst)
        if key is None:
            if classify(rel, dst) is None:          # not already sitting in a destination-layout folder
                unclaimed.append(rel)
            continue
        if new is None:
            continue
        if new in targets:
            problems.append(f"two files would land on {new}: {targets[new]} and {rel}")
            continue
        targets[new] = rel
        kind = "move"
        if (base / new).exists():
            if sha256(base / new) == sha256(base / rel):
                kind = "dedupe"
            else:
                problems.append(f"{new} already exists with different content (would overwrite {rel})")
                continue
        moves.append((key, rel, new, kind))
    return moves, unclaimed, problems


# ------------------------------------------------------------------ moving
def move_one(base: Path, rel: str, new: str, kind: str, expect_sha: str | None = None) -> tuple[str, int, str]:
    """Copy, verify, rename into place, delete source. Returns (action, bytes, sha256)."""
    s, d = base / rel, base / new
    h = sha256(s)
    if expect_sha is not None and h != expect_sha:
        raise RuntimeError(f"{rel}: SHA-256 differs from the recorded one before moving; not touched")
    size = s.stat().st_size
    if kind == "move":
        d.parent.mkdir(parents=True, exist_ok=True)
        tmp = d.with_name(d.name + ".migrating")
        shutil.copy2(s, tmp)
        if sha256(tmp) != h:
            tmp.unlink()
            raise RuntimeError(f"{rel}: copy failed its SHA-256 check; source kept, nothing changed")
        os.replace(tmp, d)
    elif sha256(d) != h:                             # dedupe: re-check right before deleting the source
        raise RuntimeError(f"{new} changed since planning; {rel} kept")
    try:
        s.unlink()
    except OSError as exc:                           # e.g. a file held open on Windows
        return f"copied; source NOT removed ({exc.__class__.__name__})", size, h
    return ("moved" if kind == "move" else "source removed (identical copy already in place)"), size, h


def prune_empty(base: Path, rels) -> int:
    """Remove now-empty folders that held moved files, deepest first, never above the cruise folder."""
    dirs = set()
    for rel in rels:
        p = (base / rel).parent
        while p != base and base in p.parents:
            dirs.add(p)
            p = p.parent
    n = 0
    for d in sorted(dirs, key=lambda x: len(x.parts), reverse=True):
        try:
            d.rmdir()
            n += 1
        except OSError:
            pass                                     # not empty (or gone): leave it
    return n


# ------------------------------------------------------------------ one cruise
def migrate_cruise(cruise: str, src: str, dst: str, apply: bool, stamp: str) -> bool:
    cdir = pp.cruise_dir(cruise)
    stored = pp.cruise_layout(cruise)
    print(f"\n=== {cruise}  (stored as {stored}; {src} -> {dst})")
    if stored == dst:
        print("  already in the target layout; nothing to do.")
        return True
    rels = list(level_files(cdir))
    moves, unclaimed, problems = plan(cdir, rels, src, dst)
    by_key: dict[str, int] = {}
    for key, *_ in moves:
        by_key[key] = by_key.get(key, 0) + 1
    for key in sorted(by_key):
        print(f"  {key:22s} {by_key[key]:5d} files   {pp.LAYOUTS[src][key]}  ->  {pp.LAYOUTS[dst][key]}")
    if not moves:
        print("  no files to move.")
    if unclaimed:
        print(f"  {len(unclaimed)} file(s) belong to no product and stay where they are:")
        for rel in unclaimed[:25]:
            print(f"      {rel}")
        if len(unclaimed) > 25:
            print(f"      ... and {len(unclaimed) - 25} more")
    if problems:
        print(f"  REFUSED: {len(problems)} conflict(s); nothing in this cruise was moved:")
        for p in problems:
            print(f"      {p}")
        return False
    if not apply:
        return True

    log = cdir / "_audit" / f"layout_migration_{stamp}.csv"
    log.parent.mkdir(parents=True, exist_ok=True)
    clean = True
    with open(log, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=LOG_FIELDS)
        w.writeheader()
        for key, rel, new, kind in moves:
            try:
                action, size, h = move_one(cdir, rel, new, kind)
            except Exception as exc:                 # record and carry on; this cruise will not be marked
                action, size, h, clean = f"FAILED: {exc}", "", "", False
            if not action.startswith(("moved", "source removed")):
                clean = False
            w.writerow(dict(product=key, source=rel, destination=new, bytes=size, sha256=h, action=action))
            fh.flush()
    pruned = prune_empty(cdir, [r for _, r, _, _ in moves])
    print(f"  moved {len(moves)} files, removed {pruned} empty folders; log -> {log.relative_to(pp.project_root())}")
    if clean:
        (cdir / pp.LAYOUT_MARKER).write_text(json.dumps(
            {"layout": dst, "migrated_utc": stamp, "files_moved": len(moves), "log": log.name}, indent=2),
            encoding="utf-8")
        print(f"  marked {cruise} as {dst}")
    else:
        print("  NOT marked: some files did not move cleanly (see the log); fix and re-run, it resumes.")
    return clean


# ------------------------------------------------------------------ frozen baseline
def migrate_baseline(bdir: Path, cruises, src: str, dst: str, apply: bool) -> bool:
    mf = bdir / "manifest.json"
    if not mf.is_file():
        return True
    man = json.loads(mf.read_text(encoding="utf-8"))
    if man.get("cruise") not in cruises:
        return True
    if man.get("layout", "v1") == dst:
        print(f"\n=== baseline ({man['cruise']}): already {dst}")
        return True
    moves, problems = [], []
    for f in man["files"]:
        key, new = relocate(f["path"], src, dst)
        if key is None:
            problems.append(f"{f['path']}: belongs to no product in {src}")
        elif new is not None:
            moves.append((f, key, new))
    print(f"\n=== baseline ({man['cruise']}): {len(moves)} of {len(man['files'])} frozen files move")
    for f, key, new in moves:
        print(f"  {f['path']}  ->  {new}")
    if problems:
        print("  REFUSED:", *problems, sep="\n      ")
        return False
    if not apply:
        return True
    for f, key, new in moves:
        kind = "dedupe" if (bdir / new).exists() else "move"
        move_one(bdir, f["path"], new, kind, expect_sha=f["sha256"])
        f["path"], f["product"] = new, key
    prune_empty(bdir, [f["path"] for f, _, _ in moves])
    shutil.copy2(mf, bdir / f"manifest.{src}.json")
    man["layout"] = dst
    mf.write_text(json.dumps(man, indent=2), encoding="utf-8")
    print(f"  manifest rewritten (previous kept as manifest.{src}.json)")
    return True


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--to", default="v2", choices=sorted(pp.LAYOUTS))
    ap.add_argument("--cruise", action="append", help="cruise id (repeatable); default: every cruise")
    ap.add_argument("--apply", action="store_true", help="move files (default is a dry run)")
    ap.add_argument("--baseline-dir", type=Path, default=REPO / "tests" / "baseline" / "local")
    a = ap.parse_args(argv)
    dst = a.to
    src = next(lay for lay in pp.LAYOUTS if lay != dst)
    croot = pp.project_root() / "cruises"
    cruises = a.cruise or (sorted(p.name for p in croot.iterdir() if p.is_dir()) if croot.is_dir() else [])
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    print(("APPLY" if a.apply else "DRY RUN (nothing moves; add --apply)") + f"  {src} -> {dst}  in {croot}")
    ok = True
    for c in cruises:
        ok &= migrate_cruise(c, src, dst, a.apply, stamp)
    ok &= migrate_baseline(a.baseline_dir, cruises, src, dst, a.apply)
    if a.apply and ok and pp.active_layout() != dst:
        print(f"\nData is now {dst}; the code still says ACTIVE_LAYOUT = {pp.active_layout()!r}. "
              f"Steps will refuse to run on migrated cruises until the layout is switched in project_paths.py.")
    print("\nDONE" if ok else "\nFINISHED WITH PROBLEMS (see above)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

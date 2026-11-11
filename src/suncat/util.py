"""Small shared helpers used across the pipeline notebooks."""
from __future__ import annotations


def pick(df, *cands):
    """First column name in df matching a candidate: exact match first, then case-insensitive
    substring. Candidates may be passed as varargs (pick(df, "prdM", "prDM")) or as a single
    tuple/list (pick(df, ("prdM", "prDM"))). Returns None if df is None or nothing matches.

    Consolidates the three notebook copies (step00a's list form, step00b/step00c's varargs form).
    """
    if df is None:
        return None
    if len(cands) == 1 and isinstance(cands[0], (list, tuple)):
        cands = tuple(cands[0])
    for c in cands:
        if c in df.columns:
            return c
    for c in df.columns:
        if any(x.lower() in c.lower() for x in cands):
            return c
    return None

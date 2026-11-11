from __future__ import annotations
import hashlib
from pathlib import Path
from typing import Callable
import numpy as np

# Producers let the regression test regenerate a frozen product from source and compare it to the
# reference. Each science module calls register_producer(...) as it is extracted, so the baseline
# test becomes a real gate the moment a module lands, and simply skips that key until then.
_PRODUCERS: "dict[str, Callable]" = {}


def register_producer(key: str, fn: Callable) -> None:
    _PRODUCERS[key] = fn


def producer(key: str):
    return _PRODUCERS.get(key)


def file_sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def within_tol(ref, got, atol: float = 0.0, rtol: float = 1e-6) -> bool:
    return bool(np.allclose(np.asarray(ref, float), np.asarray(got, float),
                            atol=atol, rtol=rtol, equal_nan=True))

from __future__ import annotations
import os
from dataclasses import dataclass
from pathlib import Path

_ENV_ROOT = "SUNCAT_CTD_ROOT"


@dataclass(frozen=True)
class CruiseConfig:
    """All pipeline paths for one cruise, derived from a single root. Replaces the hardcoded
    C:\\Projects\\ctd_pipeline scattered through the notebooks, so nothing is tied to one machine."""
    ctd_root: Path
    cruise_id: str

    @property
    def cruise_dir(self) -> Path: return self.ctd_root / "cruises" / self.cruise_id
    @property
    def L0(self) -> Path: return self.cruise_dir / "L0" / "CTD"
    @property
    def L1(self) -> Path: return self.cruise_dir / "L1" / "CTD"
    @property
    def L2(self) -> Path: return self.cruise_dir / "L2" / "CTD"
    @property
    def L3(self) -> Path: return self.cruise_dir / "L3" / "CTD"
    @property
    def align_dir(self) -> Path: return self.L2 / "Align_CTD"
    @property
    def qc_plots(self) -> Path: return self.L2 / "qc_plots"
    @property
    def psa_dir(self) -> Path: return self.ctd_root / "psa"
    @property
    def work(self) -> Path: return self.ctd_root / "_sbe_work"


def ctd_root(explicit: "str | os.PathLike | None" = None) -> Path:
    """Resolve the pipeline root. Priority: explicit argument, then $SUNCAT_CTD_ROOT, then the
    current working directory. No path is hardcoded, so the package is portable."""
    if explicit is not None:
        return Path(explicit).expanduser().resolve()
    env = os.environ.get(_ENV_ROOT)
    if env:
        return Path(env).expanduser().resolve()
    return Path.cwd().resolve()


def for_cruise(cruise_id: str, root: "str | os.PathLike | None" = None) -> CruiseConfig:
    return CruiseConfig(ctd_root=ctd_root(root), cruise_id=cruise_id)

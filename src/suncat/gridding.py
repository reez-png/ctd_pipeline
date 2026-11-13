"""L3 gridding: calibrated, flagged L2 products onto the SBE 1 m down-cast pressure bins (docs/DATA_LEVELS.md).

One table per cast. The rows ARE the bins Sea-Bird Bin Average wrote for that cast (`<cast>_1m_down.cnv`,
pressure bins, down-cast, bad scans excluded), so every CTD column is carried exactly as SBE averaged it and
the nitrate and chlorophyll columns are averaged onto the same bins the same way:

  bin      b = floor(p / binsize + 1/2): bin b spans [b - 1/2, b + 1/2) * binsize around the SBE bin centre
           b * binsize. A pressure exactly on an edge (x.500 dbar at the 0.001 dbar resolution, about 1 frame
           in 1000) goes up; Sea-Bird's own edge rule is not documented, so such a frame may sit one bin away
           from where Bin Average would put it
  value    the arithmetic mean of the frames in the bin (what Bin Average computes for the CTD columns)
  n        the number of frames averaged; reported for every bin so thin bins are visible
  sd       the sample standard deviation (ddof = 1) where n >= 2, else empty; with n of 2 or 3 it is a range
           indicator, not a reliable estimate of the spread

Nitrate frames are used only if they are on the down-cast leg (up to the scan of maximum pressure among the
scans Loop Edit kept, as for the CTD bins), coincide with a CTD scan Bin Average used (step03's ctd_scan_kept:
the 4 Hz Derived-file scan at the same moment was kept by Loop Edit), passed SUNA spectral QC (L1_clean), and
have a finite calibrated value (no3_qc). Each frame is binned by the pressure of that same scan
(ctd_scan_pressure_db), so a nitrate frame enters a bin exactly when the CTD's own scan did. Frames whose bin
is not among the CTD's bins are counted and left out. For a cast without per-scan status (no Derived file) the
fallback is the 1 s record: frames within max_ctd_gap_s of a CTD second, binned by its pressure; this catches
gaps in the 1 s record but not Loop Edit's reversal flags, and a note says so.

The CTD rows are keyed by their bin: some Bin Average runs write the bin centre as the pressure (P45_06), others
the mean pressure of the scans in the bin (P45_07). The bin centre goes in pressure_db and the CTD's own value
in ctd_pressure_db; a CTD pressure half a bin or more from its centre, or two rows in one bin, is refused. Flags travel with the values: the
nitrate's lag_limited flag, the SUNA clock offset step02 applied, and a non-photochemical-quenching risk
flag for chlorophyll in daylit near-surface bins (the criterion step07 uses to exclude bottles from its fit).

Pure functions on DataFrames plus thin readers; the pipeline notebook resolves the folders.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from suncat.io.cnv import read_cnv

_META_RE = {
    "bintype": re.compile(r"#\s*binavg_bintype\s*=\s*(\S+)", re.I),
    "binsize": re.compile(r"#\s*binavg_binsize\s*=\s*([0-9.eE+-]+)", re.I),
    "excl_bad": re.compile(r"#\s*binavg_excl_bad_scans\s*=\s*(\S+)", re.I),
    "start_time": re.compile(r"#\s*start_time\s*=\s*([A-Za-z]{3}\s+\d{1,2}\s+\d{4}\s+\d{2}:\d{2}:\d{2})", re.I),
}
_P_COL = re.compile(r"pr[dDsS]?M", re.I)


@dataclass
class GridConfig:
    npq_depth_dbar: float = 10.0     # step07 NPQ_DEPTH
    day_start_hour: int = 7          # step07 DAY_START, local hours (Ghana: local = UTC)
    day_end_hour: int = 17           # step07 DAY_END (exclusive)
    utc_offset_hours: float = 0.0    # local time = UTC + this
    min_frames: int = 1              # minimum frames for a bin value (SBE Bin Average MinScansPerBin)
    max_ctd_gap_s: float = 1.0       # fallback only: a frame must lie within this of a CTD 1 s record


@dataclass
class CastStats:
    cast_id: str
    ctd_bins: int = 0
    nitrate_frames: int = 0
    nitrate_used: int = 0
    nitrate_dropped: dict = field(default_factory=dict)
    nitrate_bins: int = 0
    nitrate_n_median: float = np.nan
    nitrate_n_max: int = 0
    ctd_bin_offset_max: float = np.nan
    chl_rows: int = 0
    chl_used: int = 0
    chl_bins: int = 0
    notes: list = field(default_factory=list)

    def row(self) -> dict:
        d = {k: v for k, v in self.__dict__.items() if k not in ("nitrate_dropped", "notes")}
        for k, v in self.nitrate_dropped.items():
            d[f"nitrate_dropped_{k}"] = v
        d["notes"] = "; ".join(self.notes)
        return d


# ------------------------------------------------------------------ readers
def read_sbe_bins(path) -> tuple[pd.DataFrame, dict]:
    """A Bin Average .cnv as a DataFrame (bad-flag sentinel masked) plus its binning metadata."""
    df = read_cnv(path, mask_badflag=True)
    meta = {"bintype": None, "binsize": np.nan, "excl_bad": None, "start_time": pd.NaT}
    for ln in Path(path).read_text(encoding="latin-1", errors="replace").splitlines():
        if ln.strip().upper().startswith("*END*"):
            break
        for k, rx in _META_RE.items():
            m = rx.search(ln)
            if m:
                meta[k] = m.group(1)
    meta["binsize"] = float(meta["binsize"]) if meta["binsize"] not in (None, np.nan) else np.nan
    meta["start_time"] = pd.to_datetime(meta["start_time"], format="%b %d %Y %H:%M:%S", errors="coerce", utc=True)
    return df, meta


def pressure_column(df: pd.DataFrame) -> str:
    for c in df.columns:
        if _P_COL.fullmatch(str(c).strip()):
            return c
    raise ValueError(f"no pressure column (prdM/prDM/prSM) among {list(df.columns)}")


def _flag(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])


# ------------------------------------------------------------------ binning
def bin_ids(p, binsize: float) -> np.ndarray:
    """Pressure-bin index of each pressure: floor(p / binsize + 1/2) (float array; NaN where p is NaN)."""
    p = np.asarray(p, dtype=float)
    return np.floor(p / binsize + 0.5)


def bin_stats(p, values, binsize: float, min_frames: int = 1) -> pd.DataFrame:
    """Per-bin mean, n and sample SD (ddof 1, NaN where n < 2) of finite values. Index: bin id."""
    d = pd.DataFrame({"b": bin_ids(p, binsize), "v": pd.to_numeric(pd.Series(values), errors="coerce").to_numpy()})
    d = d.dropna()
    if d.empty:
        return pd.DataFrame(columns=["mean", "n", "sd"], index=pd.Index([], name="b", dtype=float))
    g = d.groupby("b")["v"].agg(["mean", "count", "std"]).rename(columns={"count": "n", "std": "sd"})
    g.loc[g["n"] < 2, "sd"] = np.nan
    return g[g["n"] >= min_frames]


def npq_risk(bin_p: np.ndarray, start_time, cfg: GridConfig) -> pd.Series:
    """True where a bin is at or above npq_depth on a cast that started in local daylight; NA if the time
    is unknown. Same criterion step07 uses to drop bottles from its chlorophyll fit."""
    if start_time is None or pd.isna(start_time):
        return pd.Series(pd.NA, index=range(len(bin_p)), dtype="boolean")
    hour = (start_time + pd.Timedelta(hours=cfg.utc_offset_hours)).hour
    daylit = cfg.day_start_hour <= hour < cfg.day_end_hour
    return pd.Series(daylit & (np.asarray(bin_p, float) <= cfg.npq_depth_dbar), dtype="boolean")


# ------------------------------------------------------------------ one cast
def grid_cast(cast_id: str, ctd: pd.DataFrame, meta: dict, nitrate: pd.DataFrame | None = None,
              chl: pd.DataFrame | None = None, suna_offset_s: float | None = None,
              cfg: GridConfig | None = None) -> tuple[pd.DataFrame, CastStats]:
    cfg = cfg or GridConfig()
    st = CastStats(cast_id)
    if str(meta.get("bintype", "")).lower() != "decibars":
        raise ValueError(f"{cast_id}: the 1 m file is binned by {meta.get('bintype')!r}, not decibars; "
                         "the gridder only supports SBE pressure bins")
    bs = float(meta["binsize"])
    pcol = pressure_column(ctd)
    p = pd.to_numeric(ctd[pcol], errors="coerce").to_numpy(float)
    b_ctd = bin_ids(p, bs)
    offset = np.abs(p - b_ctd * bs)
    st.ctd_bin_offset_max = float(np.nanmax(offset)) if offset.size else np.nan
    if np.nanmax(offset) >= 0.5 * bs or pd.Series(b_ctd).duplicated().any():
        raise ValueError(f"{cast_id}: the 1 m file's pressures do not map one-to-one onto {bs:g} dbar bins "
                         f"(largest distance from a bin centre {np.nanmax(offset):.3f} dbar)")
    out = ctd.drop(columns=[c for c in ctd.columns if c in (pcol, "flag")]).copy()
    if "nbin" in out:
        out = out.rename(columns={"nbin": "ctd_n_scans"})
    # SBE writes sigma-theta as "sigma-" + a theta glyph in its own code page, which arrives as mojibake;
    # give it a plain name so the column survives every encoding.
    out = out.rename(columns={c: re.sub(r"^sigma-.+?(\d\d)$", r"sigma-theta\1", c)
                              for c in out.columns if c.startswith("sigma-") and not c.isascii()})
    out.insert(0, "ctd_pressure_db", p)
    out.insert(0, "pressure_db", b_ctd * bs)
    out.insert(0, "cast_id", cast_id)
    st.ctd_bins = len(out)

    # ---- nitrate (L2 nitrate_qc, one row per SUNA frame)
    cols = ["nitrate_uM", "nitrate_n", "nitrate_sd", "nitrate_lag_limited", "suna_time_offset_s"]
    if nitrate is None or nitrate.empty:
        st.notes.append("no nitrate product")
        for c in cols:
            out[c] = np.nan
        out["nitrate_n"] = 0
    else:
        need = {"pressure_db", "leg", "ctd_gap_s", "L1_clean", "no3_qc"}
        if not need <= set(nitrate.columns):
            raise ValueError(f"{cast_id}: nitrate product lacks {sorted(need - set(nitrate.columns))}; "
                             "re-run step03 (it writes pressure_db, leg, ctd_gap_s and the ctd_scan columns "
                             "since Nov 2026)")
        nf = nitrate.copy()
        st.nitrate_frames = len(nf)
        down = nf["leg"].astype(str).str.lower().eq("down")
        if "ctd_scan_kept" in nf and nf["ctd_scan_kept"].notna().any():
            kept = _flag(nf["ctd_scan_kept"])
            pz = pd.to_numeric(nf["ctd_scan_pressure_db"], errors="coerce")
        else:
            kept = pd.to_numeric(nf["ctd_gap_s"], errors="coerce") <= cfg.max_ctd_gap_s
            pz = pd.to_numeric(nf["pressure_db"], errors="coerce")
            st.notes.append("no per-scan Loop Edit status; 1 s record gap rule used (misses reversal flags)")
        clean = _flag(nf["L1_clean"])
        finite = pd.to_numeric(nf["no3_qc"], errors="coerce").notna() & pz.notna()
        nb = bin_ids(pz, bs)
        inside = pd.Series(np.isin(nb, b_ctd), index=nf.index)
        st.nitrate_dropped = {
            "upcast": int((~down).sum()),
            "ctd_not_kept": int((down & ~kept).sum()),
            "not_clean": int((down & kept & ~clean).sum()),
            "no_value": int((down & kept & clean & ~finite).sum()),
            "outside_ctd_bins": int((down & kept & clean & finite & ~inside).sum()),
        }
        use = down & kept & clean & finite & inside
        st.nitrate_used = int(use.sum())
        s = bin_stats(pz[use], nf.loc[use, "no3_qc"], bs, cfg.min_frames)
        idx = pd.Index(b_ctd)
        out["nitrate_uM"] = s["mean"].reindex(idx).to_numpy()
        out["nitrate_n"] = s["n"].reindex(idx).fillna(0).astype(int).to_numpy()
        out["nitrate_sd"] = s["sd"].reindex(idx).to_numpy()
        lag = bool(_flag(nf.loc[use, "lag_limited"]).any()) if "lag_limited" in nf and use.any() else False
        out["nitrate_lag_limited"] = pd.Series(np.where(out["nitrate_n"] > 0, lag, pd.NA), dtype="boolean").to_numpy()
        out["suna_time_offset_s"] = float(suna_offset_s) if suna_offset_s is not None else 0.0
        st.nitrate_bins = int((out["nitrate_n"] > 0).sum())
        st.nitrate_n_median = float(np.median(out.loc[out["nitrate_n"] > 0, "nitrate_n"])) if st.nitrate_bins else np.nan
        st.nitrate_n_max = int(out["nitrate_n"].max()) if len(out) else 0
        if st.nitrate_frames and not st.nitrate_used:
            st.notes.append("no nitrate frame passed (down-cast, CTD kept, spectral QC, finite)")

    # ---- chlorophyll (L2 chl_calibrated, rows of this cast)
    if chl is None or chl.empty:
        out["chl_mg_m3"], out["chl_n"] = np.nan, 0
        st.notes.append("no chlorophyll product")
    else:
        c = chl[chl["cast"].astype(str) == cast_id]
        st.chl_rows = len(c)
        ok = pd.to_numeric(c["chl_corrected"], errors="coerce").notna()
        if "sensor_blank" in c:
            ok &= ~_flag(c["sensor_blank"])
        if "leg" in c:                                   # native-resolution product (planned step07 change)
            ok &= c["leg"].astype(str).str.lower().eq("down")
        if "flag" in c:                                  # loop-edit flag, as Bin Average excludes it
            ok &= pd.to_numeric(c["flag"], errors="coerce").fillna(0) == 0
        st.chl_used = int(ok.sum())
        s = bin_stats(c.loc[ok, "pressure_db"], c.loc[ok, "chl_corrected"], bs, cfg.min_frames)
        out["chl_mg_m3"] = s["mean"].reindex(pd.Index(b_ctd)).to_numpy()
        out["chl_n"] = s["n"].reindex(pd.Index(b_ctd)).fillna(0).astype(int).to_numpy()
        st.chl_bins = int((out["chl_n"] > 0).sum())
    out["chl_npq_risk"] = npq_risk(out["pressure_db"].to_numpy(), meta.get("start_time"), cfg).to_numpy()
    return out, st


# ------------------------------------------------------------------ a cruise
def grid_cruise(cruise_id: str, ctd_dir: Path, nitrate_dir: Path | None, chl_file: Path | None,
                suna_l2_dir: Path | None, out_dir: Path | None, cfg: GridConfig | None = None):
    """Grid every cast that has a `<cast>_1m_down.cnv`. Returns (tables by cast, stats DataFrame).
    Writes `<cast>_L3_1m_down.csv` into out_dir when out_dir is given."""
    cfg = cfg or GridConfig()
    chl = pd.read_csv(chl_file) if chl_file is not None and Path(chl_file).exists() else None
    tables, stats = {}, []
    for f in sorted(Path(ctd_dir).glob("*_1m_down.cnv")):
        cast = f.name[: -len("_1m_down.cnv")]
        ctd, meta = read_sbe_bins(f)
        nf = None
        if nitrate_dir is not None and (Path(nitrate_dir) / f"{cast}_nitrate_QC.csv").exists():
            nf = pd.read_csv(Path(nitrate_dir) / f"{cast}_nitrate_QC.csv")
        offset = None
        if suna_l2_dir is not None and (Path(suna_l2_dir) / f"{cast}_SUNA_1s.csv").exists():
            m = pd.read_csv(Path(suna_l2_dir) / f"{cast}_SUNA_1s.csv", nrows=1)
            if "suna_time_offset_s" in m and m["suna_time_offset_s"].notna().any():
                offset = float(m["suna_time_offset_s"].iloc[0])
        table, st = grid_cast(cast, ctd, meta, nf, chl, offset, cfg)
        table.insert(0, "cruise_id", cruise_id)
        tables[cast] = table
        stats.append(st.row())
        if out_dir is not None:
            Path(out_dir).mkdir(parents=True, exist_ok=True)
            table.to_csv(Path(out_dir) / f"{cast}_L3_1m_down.csv", index=False, encoding="utf-8-sig")
    return tables, pd.DataFrame(stats)

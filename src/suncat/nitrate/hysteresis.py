"""Sensor response time (tau) from down/up hysteresis: a first-order lag fitted to both legs of a cast.

Moved from tools/estimate_response_time.py so that the pipeline (step02b) and the tool use one implementation.
The method, its tests and its assumptions are described in full in that tool's docstring; in short:

A first-order sensor obeys dy/dt = (x(t) - y(t)) / tau. The water is one profile N(z) for both legs, so
x(t) = N(z(t)) along the package's path. For a given tau the reading is linear in N; integrated exactly between
frames (x linear over each interval):
    y_i = a y_(i-1) + x_i [1 - r] + x_(i-1) [r - a],   a = exp(-dt/tau),   r = (tau/dt)(1 - a).
For each tau, N is fitted by penalised least squares on a depth grid, and tau is the minimum of the residual sum
of squares RSS(tau) (coarse-to-fine grid search, parabolic refinement). RSS(0) is the lag-free model ("both legs
see the same profile"), so 1 - RSS(tau)/RSS(0) is the share of the down/up mismatch the lag explains.
Cruise value: median over casts (frames are serially correlated, so each cast is one independent estimate), with
the interquartile range, a bootstrap over casts (90% interval), and a pooled estimate (RSS summed over casts).
Known bias: water displaced by d metres between the legs biases tau by about d / (v_down + v_up).

The first-order response is the model commonly used for sensor response-time correction (for oxygen optodes,
e.g. Bittig et al., 2014, Front. Mar. Sci.). The difference here is that it is run FORWARD (profile -> reading)
and fitted, not inverted: inverting it (y + tau dy/dt) amplifies frame noise by about tau/dt, which at tau = 30 s
and 1 s frames swamps the signal.
"""
from __future__ import annotations

import numpy as np

SMOOTH = 0.1          # curvature penalty weight, in (sensor units per m^2)^-2 relative to one frame's squared error
MIN_FRAMES = 30       # per leg, with a value


def _interp(grid: np.ndarray, z: np.ndarray):
    k = np.clip(np.searchsorted(grid, z, side="right") - 1, 0, len(grid) - 2)
    f = np.clip((z - grid[k]) / (grid[k + 1] - grid[k]), 0.0, 1.0)
    return k, f


def lag_design(t: np.ndarray, z: np.ndarray, grid: np.ndarray, tau: float) -> np.ndarray:
    """Matrix A with y = A @ N: the first-order sensor's reading at each frame for the profile N on grid."""
    n, G = len(t), len(grid)
    k, f = _interp(grid, z)
    A = np.zeros((n, G))
    rows = np.arange(n)
    A[rows, k] += 1 - f
    A[rows, k + 1] += f
    if tau <= 0:
        return A
    M = A                                                  # x_i = M[i] @ N
    dt = np.diff(t)
    a = np.exp(-dt / tau)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(dt > 0, tau / dt * (1 - a), 1.0)
    w1, w0 = 1 - r, r - a                                  # weights on x_i and x_(i-1), both >= 0
    U = np.empty_like(M)                                   # y_i = a_(i-1) y_(i-1) + U_i,  y_0 = x_0
    U[0] = M[0]
    U[1:] = w1[:, None] * M[1:] + w0[:, None] * M[:-1]
    # the recursion solved in closed form, chunk by chunk: y_i = sum_j exp(-(t_i - t_j) / tau) U_j, a cumulative
    # sum once each term is scaled by exp(+(t_j - t_s) / tau). Chunks span at most 30 tau so the scale stays below
    # e^30; all terms are non-negative, so the sums lose no precision. (A row-by-row loop does the same, slowly.)
    Y = np.empty_like(M)
    carry, s0 = None, 0
    while s0 < n:
        e = int(min(max(np.searchsorted(t, t[s0] + 30.0 * tau, side="right"), s0 + 1), n))
        E = np.exp(-(t[s0:e] - t[s0]) / tau)[:, None]
        blk = U[s0:e].copy()
        if carry is not None:
            blk[0] += a[s0 - 1] * carry
        Y[s0:e] = E * np.cumsum(blk / E, axis=0)
        carry, s0 = Y[e - 1], e
    return Y


def fit_profile(t, z, y, grid, tau, smooth=SMOOTH):
    """Least-squares profile N for one tau; returns (rss, N, residuals at the frames with a value)."""
    A = lag_design(t, z, grid, tau)
    ok = np.isfinite(y)
    G = len(grid)
    D = np.zeros((G - 2, G))
    i = np.arange(G - 2)
    D[i, i], D[i, i + 1], D[i, i + 2] = 1.0, -2.0, 1.0
    dz = float(np.median(np.diff(grid)))
    Ao = A[ok]
    # normal equations (A'A + smooth D'D / dz^4) N = A'y; the penalty makes the matrix positive definite as long as
    # two cells hold frames, and a tiny ridge guards the rest
    H = Ao.T @ Ao + (smooth / dz ** 4) * (D.T @ D)
    H[np.diag_indices(G)] += 1e-9 * max(float(np.trace(H)) / G, 1e-12)
    N = np.linalg.solve(H, Ao.T @ y[ok])
    res = y[ok] - Ao @ N
    return float(np.sum(res ** 2)), N, res


def refine_min(taus: np.ndarray, rss: np.ndarray):
    """Grid minimum refined by a parabola through it and its neighbours. edge=True when the minimum is the LAST
    grid value (no minimum within the grid). A minimum at tau = 0 is an estimate (no lag resolved), not an edge."""
    j = int(np.nanargmin(rss))
    if j == len(taus) - 1:
        return float(taus[j]), True
    if j == 0:
        return float(taus[0]), False
    x, yv = taus[j - 1:j + 2], rss[j - 1:j + 2]
    if np.isnan(yv).any():                        # a neighbour not evaluated (coarse-to-fine search)
        return float(taus[j]), False
    den = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
    A_ = (x[2] * (yv[1] - yv[0]) + x[1] * (yv[0] - yv[2]) + x[0] * (yv[2] - yv[1])) / den
    B_ = (x[2] ** 2 * (yv[0] - yv[1]) + x[1] ** 2 * (yv[2] - yv[0]) + x[0] ** 2 * (yv[1] - yv[2])) / den
    if A_ <= 0:
        return float(taus[j]), False
    return float(np.clip(-B_ / (2 * A_), x[0], x[2])), False


def search_tau(t, z, y, grid, taus, smooth=SMOOTH):
    """RSS over the tau grid, coarse to fine (every 4th tau and the last, then every tau within 4 steps of the
    coarse minimum; the RSS curve is smooth in tau). Returns (rss with NaN where not evaluated, tau_hat, edge)."""
    taus = np.asarray(taus, float)
    rss = np.full(len(taus), np.nan)
    for j in sorted(set(range(0, len(taus), 4)) | {len(taus) - 1}):
        rss[j] = fit_profile(t, z, y, grid, taus[j], smooth)[0]
    jc = int(np.nanargmin(rss))
    for j in range(max(0, jc - 4), min(len(taus), jc + 5)):
        if np.isnan(rss[j]):
            rss[j] = fit_profile(t, z, y, grid, taus[j], smooth)[0]
    tau_hat, edge = refine_min(taus, rss)
    return rss, tau_hat, edge


def cast_tau(t, z, y, is_down, is_up, taus, dz, smooth=SMOOTH, min_frames=MIN_FRAMES) -> dict:
    """tau for one cast. t: seconds (any origin), z: depth or pressure, y: sensor values with NaN where a frame is
    not to be fitted (not clean, not on a leg); every frame with t and z carries the sensor state forward.
    is_down / is_up: leg masks. Returns n, n_down, n_up, tau (NaN at the grid's upper edge), tau_grid, edge,
    rss0 / rss_min and rms0 / rms (lag-free and fitted), hyst0 / hyst (mean misfit up minus down, lag-free and
    fitted), curve (RSS over taus), data (t, z, y, grid) for refits, or skip (a reason) when it cannot."""
    t, z, y = (np.asarray(a, float) for a in (t, z, y))
    is_down, is_up = np.asarray(is_down, bool), np.asarray(is_up, bool)
    keep = np.isfinite(t) & np.isfinite(z)
    o = np.argsort(t[keep], kind="stable")
    t, z, y, is_down, is_up = (a[keep][o] for a in (t, z, y, is_down, is_up))
    y = np.where(is_down | is_up, y, np.nan)
    ok = np.isfinite(y)
    out = dict(n=int(ok.sum()), n_down=int((ok & is_down).sum()), n_up=int((ok & is_up).sum()), tau=np.nan,
               edge=False, rss0=np.nan, rss_min=np.nan, rms0=np.nan, rms=np.nan, hyst0=np.nan, hyst=np.nan,
               curve=None)
    if min(out["n_down"], out["n_up"]) < min_frames:
        out["skip"] = f"too few frames with a value on a leg (down {out['n_down']}, up {out['n_up']})"
        return out
    t = t - t[0]
    grid = np.arange(np.floor(z.min() / dz) * dz, np.ceil(z.max() / dz) * dz + dz / 2, dz)
    if len(grid) < 4:
        out["skip"] = "profile shorter than four grid cells"
        return out
    rss, tau_hat, edge = search_tau(t, z, y, grid, taus, smooth)
    up, dn = is_up[ok], is_down[ok]

    def hyst(res):                                # mean misfit up minus down: the hysteresis left unexplained
        return float(np.mean(res[up]) - np.mean(res[dn])) if up.any() and dn.any() else np.nan

    r0 = fit_profile(t, z, y, grid, 0.0, smooth)
    r1 = fit_profile(t, z, y, grid, tau_hat, smooth)
    n = int(ok.sum())
    out.update(tau=np.nan if edge else tau_hat, tau_grid=tau_hat, edge=edge, rss0=r0[0], rss_min=r1[0],
               rms0=float(np.sqrt(r0[0] / n)), rms=float(np.sqrt(r1[0] / n)), hyst0=hyst(r0[2]), hyst=hyst(r1[2]),
               curve=rss, data=(t, z, y, grid))
    return out


def summarise(per_cast: dict, taus, smooth=SMOOTH, n_boot: int = 2000, seed: int = 0) -> dict:
    """Cruise value from {cast: cast_tau result}: median over casts with an estimate, IQR, bootstrap over casts
    (90% interval of the median), pooled tau (RSS summed over casts; missing coarse-to-fine points are filled in
    around the pooled minimum) and the share of the lag-free RSS that the lag explains. Empty dict below 3 casts."""
    taus = np.asarray(taus, float)
    good = {c: r for c, r in per_cast.items() if np.isfinite(r.get("tau", np.nan))}
    if len(good) < 3:
        return {}
    v = np.array([r["tau"] for r in good.values()], float)
    rng = np.random.default_rng(seed)
    boot = np.median(rng.choice(v, (n_boot, len(v)), replace=True), axis=1)
    W = np.array([r["curve"] for r in good.values()], float)          # casts x taus
    tot = W.sum(axis=0)                                                 # NaN where any cast lacks the point
    jc = int(np.nanargmin(tot))
    for j in range(max(0, jc - 4), min(len(taus), jc + 5)):
        for i, r in enumerate(good.values()):
            if np.isnan(W[i, j]):
                W[i, j] = fit_profile(*r["data"], taus[j], smooth)[0]
    pooled, _ = refine_min(taus, W.sum(axis=0))
    return dict(tau_median=float(np.median(v)), tau_iqr=(float(np.percentile(v, 25)), float(np.percentile(v, 75))),
                tau_ci90=(float(np.percentile(boot, 5)), float(np.percentile(boot, 95))), n_casts=len(v),
                tau_pooled=float(pooled),
                rss_share_explained=float(1 - sum(r["rss_min"] for r in good.values())
                                          / sum(r["rss0"] for r in good.values())))

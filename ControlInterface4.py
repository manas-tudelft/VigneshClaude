# Authors: Maxim Beckers & Twan ten Have
# Load-cell redesign (FUTEK LSB200 + HX711): force F as the 4th state, elastic + inelastic
# material mechanics, measured-model LQR with gain scheduling, built-in simulator and self-test.
"""
Straining-stage control GUI with load-cell force feedback.

State x = [omega, I, d, F], input u = V_dc * D_cycle (V). Everything lives in this one file:
plant model, LQR design, material analysis, a simulated stage and a built-in self-test.

Run:
    python ControlInterface4.py                 real stage on SERIAL_PORT
    python ControlInterface4.py --port COM5     another serial port
    python ControlInterface4.py --sim           simulated stage, no hardware needed
    python ControlInterface4.py --sim --demo    scripted test in the simulator
    python ControlInterface4.py --selftest      built-in checks, no window

The only sample inputs are A_cs and L0. The displacement d comes from the encoder and the
force F from the load cell; E, the stiffness k, yield, hardening and relaxation are calculated
from them, never guessed. See LOADCELL_WIRING.md and MECHANICS_AND_CONTROL.md.
"""

import sys
import os
import math
import time
import argparse
import textwrap
import unittest
from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit
import control as ct

try:
    from PyQt5 import QtWidgets, QtCore, QtGui
    import pyqtgraph as pg
    QT_AVAILABLE = True
    _MainWindowBase = QtWidgets.QMainWindow
except ImportError:              # the self-test does not need Qt
    QtWidgets = QtCore = QtGui = pg = None
    QT_AVAILABLE = False
    _MainWindowBase = object


# -----------------------
# Configuration
# -----------------------
SERIAL_PORT = 'COM9'
BAUD_RATE = 115200

# System parameters (stage and motor, unchanged)
b = 1e-5
J = 6.20646e-6
k_t = 243.8 * 53.601e-3
L = 1e-4
R = 55.041
p = 0.0007
V_dc = 24.0
G_SCREW = p / np.pi          # sample elongation per radian of the output shaft (m/rad)

# LQR weights (existing tuning), augmented state order [omega, I, d, F, z]
Q_OMEGA = 1e-3
Q_CURRENT = 1e-3
Q_DISP = 5e8
Q_INT = 3e8
R_INPUT = 1e-3
SECONDARY_WEIGHT = 1e-3      # weight of the non-tracked output relative to the tracked one
TAU_PLACEHOLDER = 1.0        # s, force-state time constant until a hold has been fitted

# Load cell (FUTEK LSB200 + HX711) defaults
FORCE_LIMIT_N = 5.0          # LSB200 maximum load
PRELOAD_N = 0.01
HX711_RATE_HZ = 80.0

MODE_DISP = 0
MODE_FORCE = 1
MODE_NAMES = {MODE_DISP: "Displacement", MODE_FORCE: "Force", -1: "off"}
C_TRACK = {MODE_DISP: np.array([[0.0, 0.0, 1.0, 0.0]]),
           MODE_FORCE: np.array([[0.0, 0.0, 0.0, 1.0]])}
# Physically measured outputs: encoder displacement d and load-cell force F
C_MEASURED = np.array([[0.0, 0.0, 1.0, 0.0],
                       [0.0, 0.0, 0.0, 1.0]])


# -----------------------
# Plant dynamics and controller design
# -----------------------
@dataclass
class PlantModel:
    """Sample part of the plant (force row of A). All values are measured; before the
    stiffness has been measured the model is stiffness-free (k_i = k_r = 0)."""
    k_i: float = 0.0              # instantaneous stiffness (N/m)
    k_r: float = 0.0              # relaxed / tangent stiffness (N/m)
    tau: float = TAU_PLACEHOLDER  # relaxation time of the force state (s)
    regime: str = "unmeasured"
    measured: bool = False
    softening: bool = False
    k_e: float = 0.0              # measured elastic stiffness (N/m), scales the displacement-mode F weight


def plant_matrices(k_i, k_r, tau):
    """State x = [omega, I, d, F], input u = V_dc*D_cycle.
        J*domega = -b*omega + k_t*I - g*F
        L*dI     = -k_t*omega - R*I + u
        dd       = g*omega
        dF       = k_i*g*omega + (k_r/tau)*d - F/tau   (+ constant, removed by the integrator)
    With F = k*d this reduces to the 3-state model of the report."""
    g = G_SCREW
    A = np.array([
        [-b / J,    k_t / J,  0.0,        -g / J],
        [-k_t / L,  -R / L,   0.0,         0.0],
        [g,          0.0,     0.0,         0.0],
        [k_i * g,    0.0,     k_r / tau,  -1.0 / tau],
    ])
    B = np.array([[0.0], [1.0 / L], [0.0], [0.0]])
    C = np.eye(4)        # omega, d from the encoder, F from the load cell, I estimated
    D = np.zeros((4, 1))
    return A, B, C, D


def regime_parameters(regime, k_e=None, k_tan=None, k_v=0.0, tau_v=None):
    """Map the measured material state onto (k_i, k_r, tau):
        elastic, unloading, reloading:  k_i = k_e + k_v,   k_r = k_e
        plastic (on the envelope):      k_i = k_tan + k_v, k_r = k_tan
    k_e: elastic stiffness (slope of F-d), k_tan: tangent stiffness, k_v and tau_v: Maxwell arm
    from a hold fit. Returns the stiffness-free model when k_e has not been measured."""
    if k_e is None or not np.isfinite(k_e) or k_e <= 0:
        return PlantModel()
    k_v = max(float(k_v or 0.0), 0.0)
    tau = float(tau_v) if (tau_v is not None and np.isfinite(tau_v) and tau_v > 0) else TAU_PLACEHOLDER
    if regime in ("plastic", "softening") and k_tan is not None and np.isfinite(k_tan):
        floor = 0.01 * k_e
        softening = k_tan <= 0 or regime == "softening"
        k_tan = max(k_tan, floor)
        return PlantModel(k_tan + k_v, k_tan, tau, regime, True, softening, k_e)
    return PlantModel(k_e + k_v, k_e, tau, regime, True, False, k_e)


def _state_scales(k_i):
    """Typical magnitudes of omega (rad/s), I (A), d (m) and F (N), used to scale the
    PBH tests (raw rank tests fail because cond(W_r) ~ 1e20)."""
    return np.array([1.0, 0.1, 1e-4, max(k_i, 1.0) * 1e-4])


def reachability(A, B, model, tol=1e-12):
    """det W_r = g^2 k_t^3 (k_i - k_r) / (J^3 L^4 tau) and a scaled PBH test per mode."""
    det_analytic = G_SCREW**2 * k_t**3 * (model.k_i - model.k_r) / (J**3 * L**4 * model.tau)
    s = _state_scales(model.k_i)
    As = A * s[None, :] / s[:, None]
    Bs = B / s[:, None]
    modes = []
    for lam in np.linalg.eigvals(As):
        sv = np.linalg.svd(np.hstack([lam * np.eye(4) - As, Bs]), compute_uv=False)
        modes.append((lam, sv[-1] / sv[0]))
    unc = [lam for lam, ratio in modes if ratio < tol]
    return dict(det=det_analytic, rank=4 - len(unc), uncontrollable=unc,
                stabilizable=all(lam.real < 0 for lam in unc), modes=modes)


def observability(A, C, model, tol=1e-12):
    """Scaled PBH observability test for the outputs C (default: encoder d and load-cell F)."""
    s = _state_scales(model.k_i)
    As = A * s[None, :] / s[:, None]
    Cs = C * s[None, :]
    Cs = Cs / np.linalg.norm(Cs, axis=1, keepdims=True)
    unobs = []
    for lam in np.linalg.eigvals(As):
        sv = np.linalg.svd(np.vstack([lam * np.eye(4) - As, Cs]), compute_uv=False)
        if sv[-1] / sv[0] < tol:
            unobs.append(lam)
    return dict(rank=4 - len(unobs), unobservable=unobs,
                detectable=all(lam.real < 0 for lam in unobs))


def lqr_weights(mode, k_i, k_e=0.0):
    """Q for [omega, I, d, F, z]. In force mode the force weights scale with 1/k_i^2, which
    keeps the force-loop bandwidth (about 12 rad/s) the same in every material regime. In
    displacement mode the small F weight scales with the measured elastic stiffness k_e
    (scaling it with a softening tangent stiffness would make the F gain explode)."""
    if mode == MODE_DISP:
        k_ref = k_e if k_e > 0 else k_i
        q_F = SECONDARY_WEIGHT * Q_DISP / k_ref**2 if k_ref > 0 else 0.0
        return np.diag([Q_OMEGA, Q_CURRENT, Q_DISP, q_F, Q_INT])
    if k_i <= 0:
        raise ValueError("Force control needs the measured stiffness k, which has not been "
                         "measured yet. Load the sample in displacement control first.")
    return np.diag([Q_OMEGA, Q_CURRENT, SECONDARY_WEIGHT * Q_DISP, Q_DISP / k_i**2, Q_INT / k_i**2])


@dataclass
class ControllerDesign:
    mode: int
    model: PlantModel
    K: np.ndarray            # [k_omega, k_I, k_d, k_F, k_int]
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    poles_open: np.ndarray
    poles_closed: np.ndarray
    reach: dict
    obs: dict


def design_lqr(mode, model):
    """LQR with integral action on the tracked output (d in displacement mode, F in force mode)."""
    A, B, C, D = plant_matrices(model.k_i, model.k_r, model.tau)
    Ct = C_TRACK[mode]
    A_aug = np.block([[A, np.zeros((4, 1))], [-Ct, np.zeros((1, 1))]])
    B_aug = np.vstack((B, [[0.0]]))
    Q = lqr_weights(mode, model.k_i, model.k_e)
    K, _, E = ct.lqr(A_aug, B_aug, Q, np.array([[R_INPUT]]))
    return ControllerDesign(mode, model, np.asarray(K, dtype=float).flatten(), A, B, C, D,
                            np.linalg.eigvals(A), np.asarray(E),
                            reachability(A, B, model), observability(A, C_MEASURED, model))


_K0_CACHE = {}


def stiffness_free_design():
    """Displacement controller designed without any sample knowledge (k_i = k_r = 0).
    Its gains are practically those of every measured design, because the drive behaves
    as a velocity source (sample force barely changes the motor speed)."""
    if "design" not in _K0_CACHE:
        _K0_CACHE["design"] = design_lqr(MODE_DISP, PlantModel())
    return _K0_CACHE["design"]


def format_gain_command(r, K, mode):
    """Protocol v2 gain message: r,k_omega,k_I,k_d,k_F,k_int,mode"""
    return ",".join(f"{v:.8g}" for v in [r, *K]) + f",{int(mode)}\n"


def dominant_poles(poles, limit=1e3):
    """Poles slower than |Re| = limit (rad/s), slowest first."""
    poles = np.asarray(poles)
    slow = poles[np.abs(poles.real) < limit]
    return slow[np.argsort(np.abs(slow.real))]


def motor_quantities(omega, current, voltage):
    """Motor input power, mechanical power and produced torque (drive-train cross-check).
    voltage may be NaN; it is then reconstructed from V = R*I + k_t*omega."""
    omega = np.asarray(omega, dtype=float)
    current = np.asarray(current, dtype=float)
    voltage = np.where(np.isnan(voltage), R * current + k_t * omega, voltage)
    P_in = voltage * current
    P_mech = P_in - R * current**2
    tau_m = k_t * current
    return P_in, P_mech, tau_m


def sample_force(omega, tau_m):
    """Force on the sample implied by the motor torque, (pi/p)*(tau_m - b*omega).
    Only a drive-train estimate: gearbox and lead-screw friction are not modelled."""
    return (np.pi / p) * (np.asarray(tau_m) - b * np.asarray(omega))


def describe_design(design, A_cs, L0, K_sent=None):
    """Text for the Dynamics tab: model, matrices, gains, poles, reachability, observability."""
    m = design.model

    def mat(M):
        return "\n".join("  [" + "  ".join(f"{v: .3e}" for v in row) + " ]" for row in np.atleast_2d(M))

    def poles(P):
        return ", ".join(f"{z.real:.4g}{z.imag:+.4g}j" if abs(z.imag) > 1e-9 else f"{z.real:.4g}"
                         for z in P)

    def wrap(text):
        return textwrap.fill(text, 74, initial_indent="  ", subsequent_indent="  ")

    lines = [f"Controller mode : {MODE_NAMES[design.mode]}",
             f"Material regime : {m.regime}",
             "Plant model     : " + ("measured" if m.measured else "stiffness-free (nothing measured yet)"),
             f"k_i = {m.k_i:.4g} N/m   k_r = {m.k_r:.4g} N/m   tau = {m.tau:.4g} s"]
    if m.measured and A_cs > 0 and L0 > 0:
        lines.append(f"E_i = {m.k_i * L0 / A_cs / 1e9:.4g} GPa   E_r = {m.k_r * L0 / A_cs / 1e9:.4g} GPa")
    if m.softening:
        lines.append("Softening (necking): force control is ill-posed")
    lines += ["", "x = [omega, I, d, F]   u = V_dc*D_cycle",
              "A =", mat(design.A), "B^T =", mat(design.B.T),
              "C = I4 (omega, d: encoder; F: load cell; I: estimate)", "D = 0",
              "", "K = [k_omega, k_I, k_d, k_F, k_int]",
              "  designed: " + "  ".join(f"{v:.5g}" for v in design.K)]
    if K_sent is not None:
        lines.append("  sent    : " + "  ".join(f"{v:.5g}" for v in K_sent))
    r = design.reach
    lines += ["", "Open-loop poles (rad/s):", wrap(poles(design.poles_open)),
              "Closed-loop poles (rad/s):", wrap(poles(design.poles_closed)),
              "Dominant closed-loop poles (|Re| < 1000 rad/s):", wrap(poles(dominant_poles(design.poles_closed))),
              "", "det W_r = g^2 k_t^3 (k_i - k_r)/(J^3 L^4 tau)", f"        = {r['det']:.4g}",
              f"Reachable rank (scaled PBH) = {r['rank']} of 4"]
    if r["uncontrollable"]:
        lines.append("Uncontrollable mode(s): " + poles(r["uncontrollable"]))
        lines.append("  -> stable, system is stabilizable, LQR valid" if r["stabilizable"]
                     else "  -> NOT stabilizable")
        lines.append("  (a spring: force and displacement cannot be set independently)")
    o = design.obs
    lines.append(f"Observable rank from [d, F] (encoder + load cell) = {o['rank']} of 4")
    return "\n".join(lines)


# -----------------------
# Material mechanics and property extraction (inputs: d, F, A_cs, L0 only)
# -----------------------
HOLD, LOAD, RELOAD, UNLOAD = 0, 1, 2, 3
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def strain_stress(d, F, A_cs, L0, d0=0.0):
    """Engineering strain (d - d0)/L0 and stress F/A_cs, plus true strain and stress."""
    eps = (np.asarray(d, dtype=float) - d0) / L0
    sig = np.asarray(F, dtype=float) / A_cs
    one_plus = np.clip(1.0 + eps, 1e-9, None)
    return eps, sig, np.log(one_plus), sig * one_plus


def windowed_rate(t, y, window=0.3):
    """dy/dt from the difference over a centred time window (robust to encoder steps)."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(t)
    rate = np.zeros(n)
    if n < 2:
        return rate
    lo = np.clip(np.searchsorted(t, t - window / 2, side="left"), 0, n - 1)
    hi = np.clip(np.searchsorted(t, t + window / 2, side="right") - 1, 0, n - 1)
    dt = t[hi] - t[lo]
    ok = dt > 0
    rate[ok] = (y[hi[ok]] - y[lo[ok]]) / dt[ok]
    return rate


def displacement_quantum(d):
    """Smallest displacement step in the data (about one encoder count)."""
    steps = np.abs(np.diff(np.asarray(d, dtype=float)))
    steps = steps[steps > 1e-15]
    return float(np.percentile(steps, 5)) if len(steps) else 0.0


def segment_branches(t, d, rate_tol, d_tol, window=0.3):
    """Label each sample: LOAD (straining on the envelope), RELOAD (straining below the
    previous maximum), UNLOAD (d decreasing) or HOLD (|d rate| <= rate_tol)."""
    d = np.asarray(d, dtype=float)
    rate = windowed_rate(t, d, window)
    envelope = np.maximum.accumulate(d)
    labels = np.full(len(d), HOLD, dtype=np.int8)
    up = rate > rate_tol
    labels[up & (d >= envelope - d_tol)] = LOAD
    labels[up & (d < envelope - d_tol)] = RELOAD
    labels[rate < -rate_tol] = UNLOAD
    return labels, rate


def _linfit(x, y):
    slope, icpt = np.polyfit(x, y, 1)
    resid = y - (slope * x + icpt)
    ss_tot = np.sum((y - y.mean())**2)
    r2 = 1.0 - np.sum(resid**2) / ss_tot if ss_tot > 0 else 0.0
    return slope, icpt, r2, resid


def fit_elastic(eps, sig, d, F, load_mask, window_eps=None):
    """Young modulus E (slope of sigma-eps) and stiffness k (slope of F-d) on the virgin
    loading branch. Default window: sigma between 10 % and 40 % of sigma_max (first rise)."""
    idx = np.flatnonzero(load_mask)
    if len(idx) < 10:
        return None
    e, s = eps[idx], sig[idx]
    if window_eps is not None and window_eps[1] > window_eps[0]:
        sel = (e >= window_eps[0]) & (e <= window_eps[1])
        method = "user window"
    else:
        smax = s.max()
        if smax <= 0:
            return None
        sel = (s >= 0.10 * smax) & (s <= 0.40 * smax)
        above = np.flatnonzero(s > 0.40 * smax)
        if len(above):
            sel[above[0]:] = False
        method = "auto, 10-40 % of sigma_max"
    if sel.sum() < 8 or np.ptp(e[sel]) <= 0:
        return None
    E, s0, r2, resid = _linfit(e[sel], s[sel])
    if E <= 0:
        return None
    k, F_icpt, _, _ = _linfit(d[idx][sel], F[idx][sel])
    return dict(E=E, sigma0=s0, k=k, F0=F_icpt, r2=r2, n=int(sel.sum()), method=method,
                eps_lo=float(e[sel].min()), eps_hi=float(e[sel].max()),
                resid_std=float(resid.std()))


def _smooth(y, n=5):
    if len(y) < n:
        return y
    ys = np.convolve(y, np.ones(n) / n, mode="same")
    ys[:n // 2] = y[:n // 2]
    ys[-(n // 2):] = y[-(n // 2):]
    return ys


def offset_yield(eps, sig, load_mask, fit, offset=0.002):
    """Offset-method yield point: first crossing of the loading curve with
    sigma = E*(eps - offset) + sigma0 beyond the elastic fit window."""
    idx = np.flatnonzero(load_mask)
    if fit is None or len(idx) < 10:
        return None
    e, s = eps[idx], _smooth(sig[idx])
    h = s - (fit["E"] * (e - offset) + fit["sigma0"])
    cand = np.flatnonzero((e > fit["eps_hi"]) & (h <= 0))
    if len(cand) == 0 or cand[0] == 0:
        return None
    j = cand[0]
    w = h[j - 1] / (h[j - 1] - h[j]) if h[j - 1] != h[j] else 0.0
    return dict(eps=float(e[j - 1] + w * (e[j] - e[j - 1])),
                sig=float(s[j - 1] + w * (s[j] - s[j - 1])), offset=offset, index=int(idx[j]))


def proportional_limit(eps, sig, load_mask, fit):
    """First sustained deviation (5 points) from the elastic line beyond the fit window."""
    idx = np.flatnonzero(load_mask)
    if fit is None or len(idx) < 10:
        return None
    e, s = eps[idx], sig[idx]
    dev = np.abs(s - (fit["E"] * e + fit["sigma0"]))
    thr = np.maximum(0.01 * np.abs(s), 4 * fit["resid_std"])
    bad = (dev > thr) & (e > fit["eps_hi"])
    if len(bad) < 5:
        return None
    sustained = np.convolve(bad.astype(int), np.ones(5, dtype=int), mode="valid") == 5
    hits = np.flatnonzero(sustained)
    if len(hits) == 0:
        return None
    j = hits[0]
    return dict(eps=float(e[j]), sig=float(s[j]))


def ultimate(eps, sig):
    """UTS, uniform elongation and fracture: after the maximum, stress drops below 10 %
    of the UTS while the strain is larger than the strain at UTS."""
    if len(sig) < 5:
        return None
    i = int(np.argmax(sig))
    uts = float(sig[i])
    if uts <= 0:
        return None
    after = np.flatnonzero((np.arange(len(sig)) > i) & (sig < 0.1 * uts) & (eps > eps[i]))
    res = dict(uts=uts, eps_u=float(eps[i]), index=i, broken=False)
    if len(after):
        j = after[0]
        res.update(broken=True, eps_break=float(eps[j - 1]), sig_break=float(sig[j - 1]))
    return res


def tangent_modulus(eps, sig, load_mask, n_eval=150, rel_width=0.03):
    """Tangent modulus d(sigma)/d(eps) along the loading envelope (sliding regression)."""
    idx = np.flatnonzero(load_mask)
    if len(idx) < 20:
        return None
    x, y = eps[idx], sig[idx]
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    if np.ptp(x) <= 0:
        return None
    hw = rel_width * (x[-1] - x[0])

    def cum(a):
        return np.concatenate(([0.0], np.cumsum(a)))

    cx, cy, cxx, cxy = cum(x), cum(y), cum(x * x), cum(x * y)
    centers = np.linspace(x[0] + hw, x[-1] - hw, n_eval)
    lo = np.searchsorted(x, centers - hw)
    hi = np.searchsorted(x, centers + hw)
    n = (hi - lo).astype(float)
    Sx, Sy = cx[hi] - cx[lo], cy[hi] - cy[lo]
    Sxx, Sxy = cxx[hi] - cxx[lo], cxy[hi] - cxy[lo]
    den = n * Sxx - Sx**2
    ok = (n >= 5) & (den > 0)
    slope = np.full(n_eval, np.nan)
    slope[ok] = (n[ok] * Sxy[ok] - Sx[ok] * Sy[ok]) / den[ok]
    return centers, slope


def plastic_strain(eps, sig, fit):
    """Elastic-plastic decomposition eps_p = eps - (sigma - sigma0)/E."""
    return eps - (sig - fit["sigma0"]) / fit["E"]


def hardening(eps, sig, load_mask, fit, yld, ult):
    """Linear hardening modulus H = d(sigma)/d(eps_p) and mean tangent modulus E_t
    on the loading envelope between yield and UTS."""
    if fit is None or yld is None:
        return None
    hi = ult["eps_u"] if ult else np.inf
    sel = load_mask & (eps > yld["eps"]) & (eps <= hi)
    if sel.sum() < 10 or np.ptp(eps[sel]) <= 0:
        return None
    ep = plastic_strain(eps[sel], sig[sel], fit)
    if np.ptp(ep) <= 0:
        return None
    H, _, r2, _ = _linfit(ep, sig[sel])
    E_t, _, _, _ = _linfit(eps[sel], sig[sel])
    return dict(H=H, E_t=E_t, r2=r2, n=int(sel.sum()))


def hollomon_fit(eps, sig, load_mask, yld, ult):
    """Hollomon law sigma_true = K * eps_true^n between yield and UTS (log-log fit)."""
    if yld is None:
        return None
    hi = ult["eps_u"] if ult else np.inf
    sel = load_mask & (eps > yld["eps"]) & (eps <= hi) & (eps > 0) & (sig > 0)
    if sel.sum() < 10:
        return None
    et = np.log1p(eps[sel])
    st = sig[sel] * (1.0 + eps[sel])
    if np.ptp(np.log(et)) <= 0:
        return None
    n_exp, lnK, r2, _ = _linfit(np.log(et), np.log(st))
    return dict(K=float(np.exp(lnK)), n=float(n_exp), r2=r2)


def energies(eps, sig, load_mask, fit, yld, ult):
    """Modulus of resilience sigma_y^2/(2E) and toughness (area under the loading curve)."""
    res = {}
    if fit is not None and yld is not None:
        res["resilience"] = yld["sig"]**2 / (2.0 * fit["E"])
    idx = np.flatnonzero(load_mask)
    if len(idx) >= 10:
        e, s = eps[idx], sig[idx]
        if ult and ult.get("broken"):
            keep = e <= ult["eps_break"]
            e, s = e[keep], s[keep]
        order = np.argsort(e, kind="stable")
        if len(e) >= 2:
            res["toughness"] = float(_trapezoid(s[order], e[order]))
    return res


def detect_holds(t, rate, rate_tol, min_duration=2.0, mask=None):
    """Index ranges [i0, i1) where |rate| <= rate_tol for at least min_duration seconds."""
    still = np.abs(rate) <= rate_tol
    if mask is not None:
        still &= mask
    if len(still) == 0:
        return []
    edges = np.diff(still.astype(np.int8))
    starts = list(np.flatnonzero(edges == 1) + 1)
    ends = list(np.flatnonzero(edges == -1) + 1)
    if still[0]:
        starts.insert(0, 0)
    if still[-1]:
        ends.append(len(still))
    return [(s, e) for s, e in zip(starts, ends) if t[e - 1] - t[s] >= min_duration]


def _exp_model(tt, y_inf, delta, tau):
    return y_inf + delta * np.exp(-tt / tau)


def fit_hold(t, y):
    """Fit y(t) = y_inf + delta*exp(-(t - t0)/tau) to a hold segment
    (force at constant displacement = relaxation, displacement at constant force = creep)."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(t)
    if n < 15 or t[-1] <= t[0]:
        return None
    tt = t - t[0]
    y0 = y[:max(3, n // 50)].mean()
    y_end = y[-max(3, n // 10):].mean()
    delta0 = y0 - y_end
    tau0 = tt[-1] / 3.0
    if delta0 != 0:
        target = y_end + delta0 * math.exp(-1.0)
        crossed = np.flatnonzero((y - target) * np.sign(delta0) <= 0)
        if len(crossed):
            tau0 = max(tt[crossed[0]], tt[-1] / 50.0)
    try:
        popt, _ = curve_fit(_exp_model, tt, y, p0=[y_end, delta0, tau0],
                            bounds=([-np.inf, -np.inf, 1e-3], [np.inf, np.inf, 1e5]), maxfev=20000)
    except (RuntimeError, ValueError):
        return None
    y_inf, delta, tau = (float(v) for v in popt)
    resid = y - _exp_model(tt, *popt)
    noise = float(resid.std())
    ss_tot = float(np.sum((y - y.mean())**2))
    r2 = 1.0 - float(np.sum(resid**2)) / ss_tot if ss_tot > 0 else 0.0
    significant = abs(delta) > 3.0 * max(noise, 1e-12) and tau < 5.0 * tt[-1]
    return dict(y0=y_inf + delta, y0_measured=float(y0), y_inf=y_inf, delta=delta, tau=tau,
                r2=r2, noise=noise, significant=bool(significant),
                t_start=float(t[0]), duration=float(tt[-1]))


def relaxation_results(fit, d_hold, d0, A_cs, L0, v_load=None, t_load=None, t_ref=None):
    """Relaxation force and SLS parameters from a force-hold fit at constant displacement.
    Maxwell arm stiffness from the loading ramp before the hold (rate v for t_r seconds):
    dF = k_v * v * tau * (1 - exp(-t_r/tau)). t_ref: true start of the hold; the fitted
    exponential is extrapolated back to it (hold detection starts a little late)."""
    F_inf = fit["y_inf"]
    F0 = fit["y0"]
    if t_ref is not None and t_ref < fit["t_start"]:
        F0 = F_inf + fit["delta"] * math.exp((fit["t_start"] - t_ref) / fit["tau"])
    dF = F0 - F_inf
    res = dict(F0=F0, F_inf=F_inf, dF=dF, tau=fit["tau"], r2=fit["r2"],
               significant=fit["significant"], pct=100.0 * dF / F0 if F0 > 0 else float("nan"),
               duration=fit["duration"], t_start=fit["t_start"])
    eps_h = (d_hold - d0) / L0
    res["eps_hold"] = eps_h
    if eps_h > 0:
        res["E0"] = F0 / A_cs / eps_h
        res["E_inf"] = F_inf / A_cs / eps_h
    if v_load and v_load > 0 and t_load and t_load > 0 and fit["tau"] > 0:
        k_v = dF / (v_load * fit["tau"] * (1.0 - math.exp(-t_load / fit["tau"])))
        res["k_v"] = k_v
        res["E_v"] = k_v * L0 / A_cs
    return res


def classify_regime(t, d, F, results, eps_now, slack_force, rate_tol, d_tol, window=0.5):
    """Live material regime from the most recent data and the identified properties."""
    if results.get("fractured"):
        return "fractured"
    if len(t) < 5:
        return "unmeasured"
    recent = t >= t[-1] - window
    if recent.sum() < 3:
        recent = np.zeros(len(t), dtype=bool)
        recent[-3:] = True
    if np.nanmean(F[recent]) < slack_force:
        return "slack"
    tr = t[recent]
    rate = np.polyfit(tr - tr[0], d[recent], 1)[0] if np.ptp(tr) > 0 else 0.0
    if abs(rate) <= rate_tol:
        return "hold"
    if rate < 0:
        return "unloading"
    on_envelope = d[-1] >= np.max(d) - d_tol
    yld = results.get("yield")
    if yld is None or not on_envelope:
        return "elastic"
    k_now = results.get("k_tan_now")
    if k_now is not None and k_now < 0:
        return "softening"
    return "plastic" if eps_now > yld["eps"] else "elastic"


# -----------------------
# Material model: elastoplastic backbone + Maxwell arm (simulator ground truth and overlay)
# -----------------------
@dataclass
class MaterialParams:
    """Sample parameters. The defaults are only the simulator ground truth (a soft polymer
    film like the GUI default); the GUI never uses them for a real sample."""
    E: float = 1.0e9              # Pa, elastic modulus of the backbone
    sigma_y0: float = 30e6        # Pa, initial yield stress
    Q_voce: float = 15e6          # Pa, Voce hardening saturation
    beta_voce: float = 30.0       # -, Voce hardening rate
    H: float = 60e6               # Pa, linear hardening modulus
    E_v: float = 0.3e9            # Pa, Maxwell arm modulus (viscoelastic relaxation)
    tau_v: float = 8.0            # s, Maxwell arm relaxation time
    eps_neck: float = 0.12        # -, onset of necking (damage)
    eps_fracture: float = 0.15    # -, fracture strain
    A_cs: float = 0.01e-6         # m^2
    L0: float = 10e-3             # m
    slack: float = 20e-6          # m of slack before the sample is taut


class MaterialModel:
    """1D elasto-visco-plastic sample. Rate-independent elastoplastic backbone (Voce plus
    linear isotropic hardening, radial return) in parallel with a Maxwell arm (E_v, tau_v).
    Necking is modelled with a damage factor up to fracture. A thin sample carries no
    compression: below its taut length it is slack."""

    def __init__(self, P):
        self.P = P
        self.reset()

    def reset(self):
        self.eps_p = 0.0
        self.alpha = 0.0
        self.sig_v = 0.0
        self.eps_prev = 0.0
        self.broken = False
        self.sigma = 0.0

    def yield_stress(self, a):
        P = self.P
        return P.sigma_y0 + P.Q_voce * (1.0 - math.exp(-P.beta_voce * a)) + P.H * a

    def step(self, eps, dt):
        P = self.P
        if self.broken:
            self.eps_prev = eps
            self.sigma = 0.0
            return 0.0
        sig_tr = P.E * (eps - self.eps_p)
        if sig_tr <= 0.0:
            sig_ep = 0.0                    # slack or buckled
            self.sig_v = 0.0
        else:
            if sig_tr > self.yield_stress(self.alpha):
                da = 0.0                    # radial return, Newton on the consistency condition
                for _ in range(30):
                    a = self.alpha + da
                    g = sig_tr - P.E * da - self.yield_stress(a)
                    dg = -P.E - (P.Q_voce * P.beta_voce * math.exp(-P.beta_voce * a) + P.H)
                    inc = -g / dg
                    da += inc
                    if abs(inc) < 1e-15:
                        break
                self.eps_p += da
                self.alpha += da
                sig_ep = sig_tr - P.E * da
            else:
                sig_ep = sig_tr
            if dt > 0.0 and P.tau_v > 0.0 and P.E_v > 0.0:
                ex = math.exp(-dt / P.tau_v)
                self.sig_v = self.sig_v * ex + P.E_v * P.tau_v * (1.0 - ex) * (eps - self.eps_prev) / dt
        self.eps_prev = eps
        sig = sig_ep + self.sig_v
        if eps > P.eps_neck:
            damage = min(1.0, (eps - P.eps_neck) / (P.eps_fracture - P.eps_neck)) ** 2
            sig *= 1.0 - damage
            if damage >= 1.0:
                self.broken = True
                sig = 0.0
        self.sigma = max(sig, 0.0)
        return self.sigma

    def force(self, d, dt):
        """Force (N) for the absolute grip displacement d (m)."""
        return self.P.A_cs * self.step((d - self.P.slack) / self.P.L0, dt)


def simulate_material_force(t, d, P):
    """Force predicted by the material model P for a measured displacement history."""
    m = MaterialModel(P)
    out = np.empty(len(t))
    prev = t[0] if len(t) else 0.0
    for i in range(len(t)):
        out[i] = m.force(d[i], t[i] - prev)
        prev = t[i]
    return out


def identified_material(res, A_cs, L0, d0):
    """Material model built only from the measured properties (for the prediction overlay)."""
    fit = res.get("fit")
    if not fit:
        return None
    yld, hard, relax = res.get("yield"), res.get("hard"), res.get("relax")
    H = hard["H"] if hard else 0.0
    sigma_y0 = yld["sig"] - H * yld["offset"] if yld else 1e18
    E_v = relax.get("E_v", 0.0) if (relax and relax.get("significant") and relax.get("kind") == "relaxation") else 0.0
    E_v = max(E_v or 0.0, 0.0)
    return MaterialParams(E=fit["E"], sigma_y0=sigma_y0, Q_voce=0.0, beta_voce=1.0, H=H,
                          E_v=E_v, tau_v=relax["tau"] if E_v > 0 else 1.0,
                          eps_neck=1e9, eps_fracture=2e9, A_cs=A_cs, L0=L0,
                          slack=d0 - fit["sigma0"] * L0 / fit["E"])


# -----------------------
# Simulated Arduino (stage + sample + firmware), serial protocol v2
# -----------------------
class SimulatedArduino:
    """Arduino, stage and sample in software. Speaks serial protocol v2 exactly like
    firmware/StrainingStage/StrainingStage.ino and offers the serial.Serial methods the GUI
    uses (in_waiting, readline, write, close)."""
    TS = 0.002                 # s, simulated control period (the sketch runs at 1 kHz)
    ENCODER_RES = 1e-7         # m per count (placeholder for the real encoder)
    FORCE_NOISE = 1e-4         # N rms, LSB200 + HX711 at 80 SPS
    HX711_OFFSET = 0.015       # N, load-cell zero offset before TARE
    FRICTION_CURRENT = 0.005   # A, Coulomb friction of gearbox and lead screw as no-load current
    PRELOAD_SPEED = 1e-5       # m/s
    D_LIMIT = 0.02             # m, software travel limit
    OMEGA_FILTER = 0.005       # s, encoder speed filter

    def __init__(self, params=None, speed=1.0, realtime=True, seed=1):
        from collections import deque
        self.P = params or MaterialParams()
        self.material = MaterialModel(self.P)
        self.rng = np.random.default_rng(seed)
        self.speed = float(speed)
        self.realtime = realtime
        self.t = 0.0
        self._wall_ref, self._sim_ref = time.time(), 0.0
        self.theta = self.enc_zero = 0.0
        self.omega = self.omega_f = self.u = 0.0
        self.F_true = self.F_meas = 0.0
        self.tare_offset, self.cal_gain = 0.0, 1.0
        self.flim = FORCE_LIMIT_N
        self.K0 = stiffness_free_design().K.copy()
        self.K, self.K_disp = self.K0.copy(), self.K0.copy()
        self.mode, self.r, self.z = MODE_DISP, 0.0, 0.0
        self.state = "idle"        # idle | run | jog | preload | manual
        self.ramp = None           # [target, rate, kind]
        self.jog_settle = self.jog_timer = 0.0
        self.preload_force = None
        self.manual_u = 0.0
        self.collect = None        # [kind, value, samples] for TARE / CAL
        self.next_hx = 0.0
        self._out = deque()
        self._inbuf = ""
        self.last_line_time = None
        self._print("Straining stage FW v2 (load cell, protocol v2) [simulated]")

    # ---- serial.Serial interface ----
    @property
    def in_waiting(self):
        if self.realtime:
            self._advance_realtime()
        return len(self._out)

    def readline(self):
        if not self._out:
            return b""
        self.last_line_time, line = self._out.popleft()
        return line

    def write(self, data):
        self._inbuf += data.decode() if isinstance(data, (bytes, bytearray)) else str(data)
        while "\n" in self._inbuf:
            line, self._inbuf = self._inbuf.split("\n", 1)
            if line.strip():
                self._command(line.strip())
        return len(data)

    def close(self):
        pass

    # ---- time ----
    def now(self):
        return self.t

    def advance(self, duration):
        for _ in range(int(round(duration / self.TS))):
            self._step()

    def _advance_realtime(self, max_chunk=1.0):
        target = self._sim_ref + self.speed * (time.time() - self._wall_ref)
        if target - self.t > max_chunk:      # GUI fell behind: slow the simulation down
            target = self.t + max_chunk
            self._sim_ref, self._wall_ref = target, time.time()
        while self.t < target:
            self._step()

    def drain(self):
        """All pending output lines as text (used by the self-test)."""
        lines = [ln.decode().strip() for _, ln in self._out]
        self._out.clear()
        return lines

    def _print(self, msg):
        self._out.append((self.t, (msg + "\n").encode()))

    # ---- measurements and control law (same as the sketch) ----
    def _d_meas(self):
        d = G_SCREW * (self.theta - self.enc_zero)
        return round(d / self.ENCODER_RES) * self.ENCODER_RES

    def _x_meas(self):
        w = self.omega_f
        return w, (self.u - k_t * w) / R, self._d_meas(), self.F_meas

    @staticmethod
    def _law(K, mode, r, z, x):
        w, I, d, F = x
        if mode == MODE_DISP:
            e_d, e_F, y = d - r, F, d
        else:
            e_d, e_F, y = d, F - r, F
        return -(K[0] * w + K[1] * I + K[2] * e_d + K[3] * e_F + K[4] * z), y

    def _bumpless(self, K, mode, r, u_now):
        """Choose the integrator so the new law gives the present output u_now."""
        w, I, d, F = self._x_meas()
        e_d, e_F = (d - r, F) if mode == MODE_DISP else (d, F - r)
        rest = K[0] * w + K[1] * I + K[2] * e_d + K[3] * e_F
        self.z = -(u_now + rest) / K[4] if K[4] != 0 else 0.0

    def _stop(self):
        self.state = "idle"
        self.ramp = None
        self.manual_u = 0.0

    def _control(self):
        Ts = self.TS
        st = self.state
        if st == "idle":
            return 0.0
        if st == "manual":
            return self.manual_u
        if st == "preload":
            self.r += self.PRELOAD_SPEED * Ts
        elif self.ramp is not None:
            target, rate, kind = self.ramp
            step = rate * Ts
            if abs(target - self.r) <= step:
                self.r = target
                self.ramp = None
                if kind == "ramp":
                    self._print("RAMP done")
                elif kind == "jog":
                    self._print("JOG done")
                else:
                    self.jog_settle = self.jog_timer = 0.0
            else:
                self.r += math.copysign(step, target - self.r)
        if st == "jog" and self.ramp is None:
            self.jog_timer += Ts
            near = abs(self._d_meas() - self.r) <= 2 * self.ENCODER_RES
            self.jog_settle = self.jog_settle + Ts if near else 0.0
            if self.jog_settle >= 0.2 or self.jog_timer >= 5.0:
                self._stop()
                self._print("JOG done")
                return 0.0
        K, mode = (self.K_disp, MODE_DISP) if st in ("jog", "preload") else (self.K, self.mode)
        x = self._x_meas()
        u, y = self._law(K, mode, self.r, self.z, x)
        z_new = self.z + Ts * (self.r - y)
        u_new, _ = self._law(K, mode, self.r, z_new, x)
        if abs(u_new) <= V_dc or abs(u_new) < abs(u):      # conditional integration (anti-windup)
            self.z, u = z_new, u_new
        return max(-V_dc, min(V_dc, u))

    def _hx711_sample(self):
        raw = self.F_true + self.HX711_OFFSET + self.rng.normal(0.0, self.FORCE_NOISE)
        self.F_meas = self.cal_gain * (raw - self.tare_offset)
        if self.collect is not None:
            kind, value, samples = self.collect
            samples.append(raw)
            if len(samples) >= 16:
                mean = float(np.mean(samples))
                self.collect = None
                if kind == "tare":
                    self.tare_offset = mean
                    self._print(f"TARE done, offset = {mean:.6g} N")
                elif abs(mean - self.tare_offset) > 1e-6:
                    self.cal_gain = value / (mean - self.tare_offset)
                    self._print(f"CAL done, scale = {self.cal_gain:.6g}")
                else:
                    self._print("ERR CAL: no load on the load cell")
        if self.state != "idle" and abs(self.F_meas) > self.flim:
            self._stop()
            self._print(f"FORCE LIMIT: F = {self.F_meas:.4f} N > {self.flim:.4f} N, motor stopped")
        if self.state == "preload" and self.F_meas >= self.preload_force:
            self._stop()
            self._print("Pre-tensioning complete, controller stopped")
        if self.state != "idle" and abs(self._d_meas()) > self.D_LIMIT:
            self._stop()
            self._print("DISPLACEMENT LIMIT, motor stopped")
        w, I, d, F = self._x_meas()
        self._out.append((self.t, f"{w:.6f},{I:.6e},{d:.9f},{F:.6e},{self.u:.4f}\n".encode()))

    def _step(self):
        Ts = self.TS
        if self.t >= self.next_hx - 1e-12:
            self.next_hx += 1.0 / HX711_RATE_HZ
            self._hx711_sample()
        u = self._control()
        drive = k_t * u - R * G_SCREW * self.F_true          # quasi-static motor (L/R ~ 2 us)
        fric = R * k_t * self.FRICTION_CURRENT
        omega = 0.0 if abs(drive) <= fric else (drive - math.copysign(fric, drive)) / (k_t**2 + R * b)
        self.omega = omega
        self.theta += omega * Ts
        self.F_true = self.material.force(G_SCREW * self.theta, Ts)
        self.omega_f += (omega - self.omega_f) * min(1.0, Ts / self.OMEGA_FILTER)
        self.u = u
        self.t += Ts

    # ---- commands ----
    def _command(self, line):
        up = line.upper()
        try:
            if up == "START":
                if self.state != "run":
                    self._stop()
                    self.state = "run"
                    self._bumpless(self.K, self.mode, self.r, 0.0)
                    self._print("Controller started")
            elif up == "STOP":
                self._stop()
                self._print("Controller stopped")
            elif up == "RESET":
                self._stop()
                self.enc_zero, self.z, self.r = self.theta, 0.0, 0.0
                self._print("RESET done")
            elif up == "TARE":
                self.collect = ["tare", 0.0, []]
            elif up.startswith("CAL:"):
                self.collect = ["cal", float(line[4:]), []]
            elif up.startswith("FLIM:"):
                self.flim = abs(float(line[5:]))
                self._print(f"FLIM = {self.flim:.4f} N")
            elif up.startswith("RAMP:"):
                target, rate = (float(v) for v in line[5:].split(","))
                if self.state != "run":
                    self._print("ERR RAMP: controller not running")
                else:
                    self.ramp = [target, abs(rate), "ramp"]
            elif up.startswith("R:"):
                self.r = float(line[2:])
                self.ramp = None
            elif up.startswith("JOG:"):
                direction, distance, speed = (float(v) for v in line[4:].split(","))
                self._jog(direction, abs(distance), abs(speed))
            elif up.startswith("PRELOAD:"):
                self._preload(float(line[8:].split(",")[0]))
            elif up.startswith("MANUAL:"):
                direction, mag = (float(v) for v in line[7:].split(","))
                self._stop()
                if direction != 0 and mag > 0:
                    self.state = "manual"
                    self.manual_u = math.copysign(min(mag, 100.0) / 100.0 * V_dc, direction)
            else:
                self._gains(line)
        except (ValueError, IndexError):
            self._print(f"ERR bad command: {line}")

    def _gains(self, line):
        vals = [float(v) for v in line.split(",")]
        if len(vals) != 7:
            raise ValueError
        r, K, mode = vals[0], np.array(vals[1:6]), int(round(vals[6]))
        if mode not in (MODE_DISP, MODE_FORCE):
            raise ValueError
        if mode != self.mode:
            self.ramp = None
        elif self.ramp is not None:
            r = self.r                 # gain schedule during a ramp: the ramp keeps its own reference
        if self.state == "run":
            self._bumpless(K, mode, r, self.u)
        self.K, self.mode, self.r = K, mode, r
        if mode == MODE_DISP:
            self.K_disp = K.copy()
        self._print(f"GAINS OK mode={mode}")

    def _jog(self, direction, distance, speed):
        if speed <= 0:
            raise ValueError
        sgn = 1.0 if direction > 0 else -1.0
        if self.state == "run":
            if self.mode != MODE_DISP:
                self._print("ERR JOG: stop force control first")
                return
            target = self.r + sgn * distance if distance > 0 else sgn * self.D_LIMIT
            self.ramp = [target, speed, "jog"]
            return
        self._stop()
        self.state = "jog"
        self.r = self._d_meas()
        self._bumpless(self.K_disp, MODE_DISP, self.r, 0.0)
        target = self.r + sgn * distance if distance > 0 else sgn * self.D_LIMIT
        self.ramp = [target, speed, "jog_idle"]
        self.jog_settle = self.jog_timer = 0.0

    def _preload(self, F_target):
        if self.F_meas >= F_target:
            self._print("Pre-tensioning complete, controller stopped")
            return
        self._stop()
        self.state = "preload"
        self.preload_force = F_target
        self.r = self._d_meas()
        self._bumpless(self.K_disp, MODE_DISP, self.r, 0.0)


# -----------------------
# Data store
# -----------------------
class DataStore:
    """Growable column buffers; arrays() returns views without copying."""
    FIELDS = ("t", "omega", "I", "d", "F", "u", "r", "mode")

    def __init__(self, capacity=8192):
        self._cap = capacity
        self.clear()

    def clear(self):
        self._buf = np.full((len(self.FIELDS), self._cap), np.nan)
        self._n = 0

    def __len__(self):
        return self._n

    def append(self, row):
        if self._n == self._buf.shape[1]:
            grown = np.full((len(self.FIELDS), 2 * self._buf.shape[1]), np.nan)
            grown[:, :self._n] = self._buf[:, :self._n]
            self._buf = grown
        self._buf[:, self._n] = row
        self._n += 1

    def arrays(self):
        return {f: self._buf[i, :self._n] for i, f in enumerate(self.FIELDS)}

    def last(self, field):
        return self._buf[self.FIELDS.index(field), self._n - 1] if self._n else float("nan")


def _decimate(n, max_points):
    return slice(None, None, max(1, int(math.ceil(n / max_points))))


def _rel_change(new, old):
    return abs(new - old) / max(abs(old), 1e-12)


def _spin(lo, hi, value, step, decimals):
    s = QtWidgets.QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setDecimals(decimals)
    s.setSingleStep(step)
    s.setValue(value)
    return s


# -----------------------
# GUI
# -----------------------
class StageGUI(_MainWindowBase):
    """Controls on the left; Signals, Material and Dynamics tabs on the right."""
    RATE_TOL = 5e-7        # m/s, |d rate| at or below this is a hold (raised to 4 encoder counts per window)
    RATE_WINDOW = 0.5      # s, window of the displacement-rate estimate
    D_TOL = 1e-6           # m, tolerance for being on the loading envelope
    HOLD_MIN_S = 2.0       # s, shortest hold that is fitted
    CREEP_RATE_TOL = 2e-3  # N/s, |F rate| at or below this is a force hold
    PLOT_POINTS = 4000

    def __init__(self, link, sim=None, link_name=""):
        super().__init__()
        self.link, self.sim, self.link_name = link, sim, link_name
        self.t_wall0 = time.time()
        self.data = DataStore()
        self.events = []
        self.running = False
        self.jogging = False
        self.ctrl_mode = MODE_DISP
        self.protocol = None
        self.d0, self.t_origin = 0.0, None
        self.model = PlantModel()
        self.design = stiffness_free_design()
        self.sent_K = None
        self.r_gui, self.ramp_gui = 0.0, None
        self.results = {}
        self.permanent_set = None
        self.test = None
        self.last_loading = ("elastic", None)
        self.model_hist = []
        self.pole_hist = []
        self.last_schedule = -1e9
        self.schedule_pending = False
        self.fracture_handled = False
        self.last_overlay_wall = 0.0
        self.overlay = None
        self._dyn_text_cache = ""
        self._rx = b""
        self._t_last = -1e9
        self.setWindowTitle("Motor Control Live - load cell (FUTEK LSB200 + HX711)")
        self.resize(1550, 980)
        self._build_ui()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.on_tick)
        self.timer.start(20)                 # 50 Hz serial and plots
        self.analysis_timer = QtCore.QTimer(self)
        self.analysis_timer.timeout.connect(self.run_analysis)
        self.analysis_timer.start(500)       # 2 Hz material analysis
        self.update_status_display()
        self.update_dynamics_view()
        self.send(f"FLIM:{self.flim_input.value():.4f}")

    # ---- layout ----
    def _group(self, title):
        box = QtWidgets.QGroupBox(title)
        return box, QtWidgets.QFormLayout(box)

    def _build_ui(self):
        W = QtWidgets
        left = W.QWidget()
        lay = W.QVBoxLayout(left)

        self.status_label = W.QLabel("Status: OFF")
        font = self.status_label.font()
        font.setPointSize(12)
        self.status_label.setFont(font)
        self.regime_label = W.QLabel("Regime: unmeasured")
        self.regime_label.setStyleSheet("QLabel { font-weight: bold; }")
        self.link_label = W.QLabel(f"Link: {self.link_name}")
        self.fw_label = W.QLabel("Firmware protocol: waiting for data")
        for w in (self.status_label, self.regime_label, self.link_label, self.fw_label):
            lay.addWidget(w)

        box, form = self._group("Sample (only inputs: geometry)")
        self.area_input = _spin(1e-5, 100.0, 0.01, 0.01, 5)
        self.length_input = _spin(0.001, 200.0, 10.0, 0.1, 4)
        form.addRow("Cross-sectional area A_cs (mm^2):", self.area_input)
        form.addRow("Gauge length L0 (mm):", self.length_input)
        lay.addWidget(box)

        box, form = self._group("Control")
        self.mode_combo = W.QComboBox()
        self.mode_combo.addItems(["Displacement", "Force (stress)"])
        self.mode_combo.model().item(MODE_FORCE).setEnabled(False)
        self.mode_combo.setToolTip("Force mode unlocks once the stiffness k has been measured")
        self.disp_ref_input = _spin(-50.0, 50.0, 0.0, 0.01, 5)
        self.disp_ref_strain = W.QLabel("")
        self.stress_ref_input = _spin(0.0, 1000.0, 0.0, 1.0, 4)
        self.stress_ref_input.setEnabled(False)
        self.stress_ref_force = W.QLabel("")
        form.addRow("Mode:", self.mode_combo)
        form.addRow("Displacement d (mm, from preload zero):", self.disp_ref_input)
        form.addRow("", self.disp_ref_strain)
        form.addRow("Stress (MPa, force mode):", self.stress_ref_input)
        form.addRow("", self.stress_ref_force)
        self.send_button = W.QPushButton("Send Parameters to Arduino")
        self.start_button = W.QPushButton("Start Controller")
        self.update_ref_button = W.QPushButton("Update Reference")
        self.auto_sched = W.QCheckBox("Auto gain scheduling (measured model)")
        self.auto_sched.setChecked(True)
        for w in (self.send_button, self.start_button, self.update_ref_button, self.auto_sched):
            form.addRow(w)
        lay.addWidget(box)

        box, form = self._group("Test program (displacement)")
        self.target_input = _spin(-50.0, 50.0, 0.5, 0.05, 4)
        self.rate_input = _spin(1e-5, 5.0, 0.01, 0.005, 5)
        self.rate_strain = W.QLabel("")
        self.hold_input = _spin(0.0, 1e5, 30.0, 5.0, 1)
        form.addRow("Target d (mm, from preload zero):", self.target_input)
        form.addRow("Rate (mm/s):", self.rate_input)
        form.addRow("", self.rate_strain)
        form.addRow("Hold time (s):", self.hold_input)
        self.ramp_hold_button = W.QPushButton("Ramp -> Hold (relaxation)")
        self.unload_button = W.QPushButton("Unload (to preload force)")
        self.test_label = W.QLabel("Test: idle")
        for w in (self.ramp_hold_button, self.unload_button, self.test_label):
            form.addRow(w)
        lay.addWidget(box)

        box, form = self._group("Load cell (FUTEK LSB200 + HX711)")
        self.force_label = W.QLabel("F = - N")
        font = self.force_label.font()
        font.setPointSize(11)
        font.setBold(True)
        self.force_label.setFont(font)
        self.tare_button = W.QPushButton("Tare (sample unloaded)")
        self.cal_mass_input = _spin(0.1, 1000.0, 200.0, 10.0, 2)
        self.cal_button = W.QPushButton("Calibrate with hanging mass")
        self.preload_input = _spin(0.0001, 5.0, PRELOAD_N, 0.005, 4)
        self.pretension_button = W.QPushButton("Pre-tension (preload, sets strain zero)")
        self.flim_input = _spin(0.01, 50.0, FORCE_LIMIT_N, 0.1, 3)
        self.flim_button = W.QPushButton("Set force limit")
        form.addRow(self.force_label)
        form.addRow(self.tare_button)
        form.addRow("Calibration mass (g):", self.cal_mass_input)
        form.addRow(self.cal_button)
        form.addRow("Preload force (N):", self.preload_input)
        form.addRow(self.pretension_button)
        form.addRow("Force limit (N):", self.flim_input)
        form.addRow(self.flim_button)
        lay.addWidget(box)

        box, form = self._group("Manual control (distance-based)")
        self.jog_distance_input = _spin(0.0, 50.0, 0.1, 0.01, 4)
        self.jog_speed_input = _spin(1e-4, 5.0, 0.05, 0.01, 4)
        form.addRow("Distance (mm, 0 = continuous):", self.jog_distance_input)
        form.addRow("Speed (mm/s):", self.jog_speed_input)
        row = W.QHBoxLayout()
        self.jog_out_button = W.QPushButton("Outwards (+d)")
        self.jog_in_button = W.QPushButton("Inwards (-d)")
        self.jog_stop_button = W.QPushButton("Stop")
        for btn in (self.jog_out_button, self.jog_in_button, self.jog_stop_button):
            row.addWidget(btn)
        form.addRow(row)
        lay.addWidget(box)

        box, form = self._group("Analysis settings")
        self.fit_start_input = _spin(0.0, 50.0, 0.0, 0.01, 4)
        self.fit_end_input = _spin(0.0, 50.0, 0.0, 0.01, 4)
        self.offset_input = _spin(0.0, 5.0, 0.2, 0.05, 3)
        self.true_curve_check = W.QCheckBox("Show true stress-strain")
        self.drive_check = W.QCheckBox("Show drive-train force estimate")
        form.addRow("E fit window from d (mm):", self.fit_start_input)
        form.addRow("to d (mm), equal = auto:", self.fit_end_input)
        form.addRow("Yield offset (% strain):", self.offset_input)
        form.addRow(self.true_curve_check)
        form.addRow(self.drive_check)
        lay.addWidget(box)

        box = W.QGroupBox("Results (calculated from d, F, A_cs and L0)")
        v = W.QVBoxLayout(box)
        self.results_label = W.QLabel("Not measured yet.")
        self.results_label.setFont(QtGui.QFont("Consolas", 9))
        self.results_label.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        self.results_label.setWordWrap(True)
        v.addWidget(self.results_label)
        lay.addWidget(box)

        box, form = self._group("Data")
        self.reset_button = W.QPushButton("Reset Measurement (Clear Data)")
        self.save_button = W.QPushButton("Save Data to File")
        form.addRow(self.reset_button)
        form.addRow(self.save_button)
        lay.addWidget(box)

        box = W.QGroupBox("Arduino log")
        v = W.QVBoxLayout(box)
        self.log_box = W.QPlainTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumBlockCount(500)
        self.log_box.setFixedHeight(140)
        v.addWidget(self.log_box)
        lay.addWidget(box)
        lay.addStretch()

        scroll = W.QScrollArea()
        scroll.setWidget(left)
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(460)

        self.tabs = W.QTabWidget()
        self._build_signal_tab()
        self._build_material_tab()
        self._build_dynamics_tab()

        central = W.QWidget()
        h = W.QHBoxLayout(central)
        h.addWidget(scroll, 0)
        h.addWidget(self.tabs, 1)
        self.setCentralWidget(central)

        self.send_button.clicked.connect(self.send_parameters)
        self.start_button.clicked.connect(self.toggle_controller)
        self.update_ref_button.clicked.connect(self.update_reference)
        self.mode_combo.currentIndexChanged.connect(self._mode_changed)
        self.ramp_hold_button.clicked.connect(self.start_ramp_hold)
        self.unload_button.clicked.connect(self.start_unload)
        self.tare_button.clicked.connect(self.tare)
        self.cal_button.clicked.connect(self.calibrate)
        self.pretension_button.clicked.connect(self.pre_tension)
        self.flim_button.clicked.connect(self.set_force_limit)
        self.jog_out_button.clicked.connect(lambda: self.jog(+1))
        self.jog_in_button.clicked.connect(lambda: self.jog(-1))
        self.jog_stop_button.clicked.connect(self.stop_manual_control)
        self.reset_button.clicked.connect(self.reset_data)
        self.save_button.clicked.connect(lambda: self.save_data_to_file())
        for s in (self.disp_ref_input, self.stress_ref_input, self.area_input,
                  self.length_input, self.rate_input):
            s.valueChanged.connect(self._update_hint_labels)
        self._update_hint_labels()

    def _build_signal_tab(self):
        win = pg.GraphicsLayoutWidget()
        dash = QtCore.Qt.DashLine
        self.p_d = win.addPlot(title="Displacement d (from preload zero)")
        self.p_d.addLegend()
        self.c_d = self.p_d.plot(pen='y', name="d (encoder)")
        self.c_r = self.p_d.plot(pen=pg.mkPen('w', style=dash), name="reference", connect="finite")
        self.p_d.setLabel('left', 'd (mm)')
        self.p_d.setLabel('bottom', 't (s)')
        self.p_F = win.addPlot(title="Force F (load cell)")
        self.p_F.addLegend()
        self.c_F = self.p_F.plot(pen='c', name="F (load cell)")
        self.c_F_model = self.p_F.plot(pen=pg.mkPen('m', width=2, style=dash), name="identified model")
        self.c_F_drive = self.p_F.plot(pen=pg.mkPen((150, 150, 150)), name="drive-train estimate")
        self.p_F.setLabel('left', 'F (N)')
        self.p_F.setLabel('bottom', 't (s)')
        win.nextRow()
        self.p_w = win.addPlot(title="Omega (output shaft)")
        self.c_w = self.p_w.plot(pen='g')
        self.p_w.setLabel('left', 'omega (rad/s)')
        self.p_w.setLabel('bottom', 't (s)')
        self.p_u = win.addPlot(title="Drive: voltage u and estimated current")
        self.p_u.addLegend(offset=(-10, 10))
        self.c_u = self.p_u.plot(pen='r', name="u (V)")
        self.c_I = self.p_u.plot(pen=pg.mkPen((255, 165, 0)), name="I_est (mA)")
        self.p_u.setLabel('bottom', 't (s)')
        self.tabs.addTab(win, "Signals")

    def _build_material_tab(self):
        win = pg.GraphicsLayoutWidget()
        dash, dot = QtCore.Qt.DashLine, QtCore.Qt.DotLine
        self.p_ss = win.addPlot(title="Stress-strain")
        self.p_ss.addLegend(offset=(-10, -10))
        self.c_ss = self.p_ss.plot(pen='r', name="engineering")
        self.c_ss_true = self.p_ss.plot(pen=pg.mkPen((255, 165, 0)), name="true")
        self.c_ss_fit = self.p_ss.plot(pen=pg.mkPen('w', width=2, style=dash), name="elastic fit (E)")
        self.c_ss_off = self.p_ss.plot(pen=pg.mkPen((170, 170, 170), style=dot), name="offset line")
        self.c_ss_pts = pg.ScatterPlotItem(pen=None)
        self.p_ss.addItem(self.c_ss_pts)
        self.p_ss.setLabel('left', 'Stress (MPa)')
        self.p_ss.setLabel('bottom', 'Strain (%)')
        self.p_fd = win.addPlot(title="Force-displacement")
        self.p_fd.addLegend(offset=(-10, -10))
        self.c_fd = self.p_fd.plot(pen='y', name="measured")
        self.c_fd_fit = self.p_fd.plot(pen=pg.mkPen('w', width=2, style=dash), name="stiffness k")
        self.p_fd.setLabel('left', 'F (N)')
        self.p_fd.setLabel('bottom', 'd (mm, from preload zero)')
        win.nextRow()
        self.p_et = win.addPlot(title="Tangent modulus along the loading envelope")
        self.c_et = self.p_et.plot(pen='g', connect="finite")
        self.p_et.setLabel('left', 'E_t (GPa)')
        self.p_et.setLabel('bottom', 'Strain (%)')
        self.p_hold = win.addPlot(title="Last hold: relaxation (F) or creep (d)")
        self.p_hold.addLegend(offset=(-10, 10))
        self.c_hold = self.p_hold.plot(pen='c', name="measured")
        self.c_hold_fit = self.p_hold.plot(pen=pg.mkPen('w', width=2, style=dash), name="exponential fit")
        self.p_hold.setLabel('bottom', 'time in hold (s)')
        self.tabs.addTab(win, "Material")

    def _build_dynamics_tab(self):
        split = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.dyn_text = QtWidgets.QPlainTextEdit()
        self.dyn_text.setReadOnly(True)
        self.dyn_text.setFont(QtGui.QFont("Consolas", 9))
        self.dyn_text.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        split.addWidget(self.dyn_text)
        win = pg.GraphicsLayoutWidget()
        self.p_k = win.addPlot(title="Measured plant stiffness in the force row of A")
        self.p_k.addLegend(offset=(-10, 10))
        self.c_ki = self.p_k.plot(pen='m', symbol='o', symbolSize=6, symbolBrush='m', name="k_i (instantaneous)")
        self.c_kr = self.p_k.plot(pen='c', symbol='t', symbolSize=6, symbolBrush='c', name="k_r (relaxed / tangent)")
        self.p_k.setLabel('left', 'N/m')
        self.p_k.setLabel('bottom', 't (s)')
        win.nextRow()
        self.p_pole = win.addPlot(title="Dominant closed-loop poles (red: current design)")
        self.c_pole_hist = pg.ScatterPlotItem(size=7, pen=None, brush=pg.mkBrush(120, 120, 255, 90))
        self.c_pole_now = pg.ScatterPlotItem(size=13, pen=None, brush=pg.mkBrush(255, 60, 60), symbol='x')
        self.p_pole.addItem(self.c_pole_hist)
        self.p_pole.addItem(self.c_pole_now)
        self.p_pole.setLabel('left', 'Im (rad/s)')
        self.p_pole.setLabel('bottom', 'Re (rad/s)')
        self.p_pole.showGrid(x=True, y=True, alpha=0.3)
        split.addWidget(win)
        split.setSizes([640, 600])
        self.tabs.addTab(split, "Dynamics")

    # ---- helpers ----
    def clock(self):
        """Seconds since start: simulator time in --sim, wall time on the real stage."""
        return self.sim.now() if self.sim is not None else time.time() - self.t_wall0

    def log(self, text):
        self.log_box.appendPlainText(text)
        print(text)

    def send(self, cmd):
        cmd = cmd if cmd.endswith("\n") else cmd + "\n"
        try:
            self.link.write(cmd.encode())
        except Exception as e:
            self.log(f"[serial error] {e}")
        self.log(f">> {cmd.strip()}")

    def add_event(self, name):
        self.events.append((self.clock(), name))

    def sample_geometry(self):
        """A_cs (m^2) and L0 (m): the only sample inputs."""
        return self.area_input.value() * 1e-6, self.length_input.value() * 1e-3

    def last_d(self):
        d = self.data.last("d")
        return float(d) if np.isfinite(d) else 0.0

    def last_F(self):
        F = self.data.last("F")
        return float(F) if np.isfinite(F) else 0.0

    def check_protocol(self):
        if self.protocol == "v2":
            return True
        QtWidgets.QMessageBox.warning(
            self, "Firmware",
            "No protocol-v2 data line (omega,I,d,F,u) has been received yet.\n"
            "Upload firmware/StrainingStage/StrainingStage.ino and check the serial port.")
        return False

    def _update_hint_labels(self, *_):
        A_cs, _L0 = self.sample_geometry()
        L0_mm = self.length_input.value()
        self.disp_ref_strain.setText(f"= {self.disp_ref_input.value() / L0_mm * 100:.4g} % strain")
        self.stress_ref_force.setText(f"= {self.stress_ref_input.value() * 1e6 * A_cs:.4g} N")
        self.rate_strain.setText(f"= {self.rate_input.value() / L0_mm * 100:.4g} %/s strain rate")

    def _mode_changed(self, idx):
        self.stress_ref_input.setEnabled(idx == MODE_FORCE)
        self.disp_ref_input.setEnabled(idx == MODE_DISP)

    def reference_from_inputs(self, mode):
        A_cs, _ = self.sample_geometry()
        if mode == MODE_DISP:
            return self.d0 + self.disp_ref_input.value() * 1e-3
        return self.stress_ref_input.value() * 1e6 * A_cs

    def current_reference(self, t):
        """The controller reference as the GUI expects it (ramps run in the firmware)."""
        if self.ramp_gui is None:
            return self.r_gui
        t0, r0, target, rate = self.ramp_gui
        step = rate * max(t - t0, 0.0)
        return target if abs(target - r0) <= step else r0 + math.copysign(step, target - r0)

    def update_status_display(self):
        if self.running:
            text, color = f"Status: ACTIVE ({MODE_NAMES[self.ctrl_mode]} control)", "green"
        elif self.jogging:
            text, color = "Status: JOG (closed-loop move)", "darkorange"
        else:
            text, color = "Status: OFF", "red"
        self.status_label.setText(text)
        self.status_label.setStyleSheet(f"QLabel {{ color: {color}; font-weight: bold; }}")
        self.start_button.setText("Stop Controller" if (self.running or self.jogging) else "Start Controller")
        regime = self.results.get("regime", "unmeasured")
        self.regime_label.setText(f"Regime: {regime}   |   plant model: "
                                  + ("measured" if self.model.measured else "stiffness-free (nothing measured yet)"))

    # ---- controller ----
    def apply_design(self, mode, r, reason, keep_ramp=False):
        if mode == MODE_FORCE and self.model.softening:
            QtWidgets.QMessageBox.warning(self, "Force control",
                                          "The sample is softening (necking): force control is ill-posed. "
                                          "Use displacement control.")
            return False
        try:
            design = design_lqr(mode, self.model)
        except ValueError as e:
            QtWidgets.QMessageBox.warning(self, "Controller design", str(e))
            return False
        self.design, self.ctrl_mode = design, mode
        self.sent_K = design.K.copy()
        self.send(format_gain_command(r, design.K, mode))
        now = self.clock()
        if keep_ramp and self.ramp_gui is not None:
            self.ramp_gui = (now, r, self.ramp_gui[2], self.ramp_gui[3])
        else:
            self.ramp_gui = None
        self.r_gui = r
        self.pole_hist.append((now, dominant_poles(design.poles_closed, 100.0)))
        self.add_event(f"GAINS {MODE_NAMES[mode]} ({reason})")
        self.update_dynamics_view()
        return True

    def send_parameters(self):
        if not self.check_protocol():
            return
        mode = self.mode_combo.currentIndex()
        self.apply_design(mode, self.reference_from_inputs(mode), "Send parameters")

    def toggle_controller(self):
        if self.running or self.jogging:
            self.stop_controller()
        else:
            self.start_controller()

    def start_controller(self):
        if not self.check_protocol():
            return
        if self.sent_K is None and not self.apply_design(MODE_DISP, self.last_d(), "start: hold position"):
            return
        self.send("START")
        self.running = True
        self.add_event("START")
        self.update_status_display()

    def stop_controller(self):
        self.send("STOP")
        self.running = self.jogging = False
        self.test, self.ramp_gui = None, None
        self.test_label.setText("Test: idle")
        self.add_event("STOP")
        self.update_status_display()

    def update_reference(self):
        if self.mode_combo.currentIndex() != self.ctrl_mode:
            QtWidgets.QMessageBox.information(self, "Update reference",
                                              "The mode was changed: press 'Send Parameters to Arduino' first.")
            return
        r = self.reference_from_inputs(self.ctrl_mode)
        self.send(f"R:{r:.9g}")
        self.r_gui, self.ramp_gui = r, None
        self.add_event("REFERENCE")

    def maybe_reschedule(self):
        """Resend gains designed on the latest measured model (bumpless in the firmware)."""
        if not (self.schedule_pending and self.auto_sched.isChecked() and self.running
                and self.sent_K is not None):
            return
        now = self.clock()
        if now - self.last_schedule < 1.0:
            return
        self.schedule_pending = False
        self.last_schedule = now
        if self.ctrl_mode == MODE_FORCE and self.model.softening:
            self.log("Softening detected: switching from force to displacement control")
            self.apply_design(MODE_DISP, self.last_d(), "softening: displacement hold")
            return
        self.apply_design(self.ctrl_mode, self.current_reference(now),
                          f"gain schedule, {self.model.regime}", keep_ramp=True)

    # ---- test program ----
    def _ensure_displacement_control(self, reason):
        """Displacement control running, starting bumplessly at the present position."""
        if self.running and self.ctrl_mode == MODE_DISP:
            return True
        if not self.apply_design(MODE_DISP, self.last_d(), reason):
            return False
        if not self.running:
            self.send("START")
            self.running = True
            self.add_event("START")
            self.update_status_display()
        return True

    def start_ramp_hold(self):
        if not self.check_protocol():
            return
        target = self.d0 + self.target_input.value() * 1e-3
        rate = self.rate_input.value() * 1e-3
        if not self._ensure_displacement_control("test program"):
            return
        r_start = self.current_reference(self.clock())
        self.send(f"RAMP:{target:.9g},{rate:.9g}")
        self.ramp_gui = (self.clock(), r_start, target, rate)
        self.test = dict(stage="ramp", hold=self.hold_input.value())
        self.add_event(f"RAMP to d = {self.target_input.value():.4g} mm at {self.rate_input.value():.4g} mm/s")
        self.test_label.setText("Test: ramping")

    def start_unload(self):
        if not self.check_protocol():
            return
        rate = self.rate_input.value() * 1e-3
        if not self._ensure_displacement_control("unload"):
            return
        r_start = self.current_reference(self.clock())
        self.send(f"RAMP:{self.d0:.9g},{rate:.9g}")
        self.ramp_gui = (self.clock(), r_start, self.d0, rate)
        self.test = dict(stage="unload")
        self.add_event("UNLOAD")
        self.test_label.setText("Test: unloading to the preload force")

    def on_ramp_done(self):
        if self.ramp_gui is not None:
            self.r_gui, self.ramp_gui = self.ramp_gui[2], None
        if self.test is None:
            return
        if self.test["stage"] == "ramp":
            self.test.update(stage="hold", t_hold=self.clock())
            self.add_event("HOLD start")
            self.test_label.setText(f"Test: holding {self.test['hold']:.0f} s (relaxation recorded)")
        elif self.test["stage"] == "unload":
            self.finish_unload("reference reached zero elongation")

    def finish_unload(self, why):
        _, L0 = self.sample_geometry()
        d = self.last_d()
        self.send(f"R:{d:.9g}")
        self.r_gui, self.ramp_gui = d, None
        self.permanent_set = (d - self.d0) / L0
        self.test = None
        self.add_event(f"UNLOADED ({why})")
        self.test_label.setText(f"Test: unloaded, permanent set = {self.permanent_set * 100:.4g} %")

    def update_test_program(self):
        t = self.test
        if t is None:
            return
        if t["stage"] == "hold" and self.clock() - t["t_hold"] >= t["hold"]:
            self.test = None
            self.add_event("HOLD end")
            self.test_label.setText("Test: hold finished, controller holds position")
        elif t["stage"] == "unload" and len(self.data) and self.last_F() <= self.preload_input.value():
            self.finish_unload("force reached the preload")

    # ---- load cell ----
    def tare(self):
        self.send("TARE")
        self.add_event("TARE")

    def calibrate(self):
        F = self.cal_mass_input.value() * 1e-3 * 9.80665
        answer = QtWidgets.QMessageBox.question(
            self, "Calibrate load cell",
            f"Tare first with the cell unloaded. Then hang {self.cal_mass_input.value():g} g "
            f"(F = {F:.5g} N) on the LSB200 and press Yes.")
        if answer != QtWidgets.QMessageBox.Yes:
            return
        self.send(f"CAL:{F:.6g}")
        self.add_event(f"CAL {F:.5g} N")

    def pre_tension(self):
        if not self.check_protocol():
            return
        self.send(f"PRELOAD:{self.preload_input.value():.6g}")
        self.running = True
        self.test, self.ramp_gui = None, None
        self.add_event("PRETENSION")
        self.update_status_display()

    def set_force_limit(self):
        self.send(f"FLIM:{self.flim_input.value():.4f}")
        self.add_event(f"FLIM {self.flim_input.value():.4g} N")

    def set_strain_origin(self):
        if len(self.data) == 0:
            return
        self.d0, self.t_origin = self.last_d(), float(self.data.last("t"))
        self.r_gui, self.ramp_gui = self.d0, None
        self.disp_ref_input.setValue(0.0)
        self.add_event(f"STRAIN ZERO at d = {self.d0 * 1e3:.5f} mm (preload {self.last_F():.4g} N)")

    # ---- manual control (distance-based) ----
    def jog(self, direction):
        if not self.check_protocol():
            return
        if self.running and self.ctrl_mode == MODE_FORCE:
            QtWidgets.QMessageBox.information(self, "Manual control", "Stop force control before jogging.")
            return
        dist = self.jog_distance_input.value() * 1e-3
        speed = self.jog_speed_input.value() * 1e-3
        self.send(f"JOG:{direction},{dist:.9g},{speed:.9g}")
        if self.running:
            r0 = self.current_reference(self.clock())
            target = r0 + direction * dist if dist > 0 else direction * SimulatedArduino.D_LIMIT
            self.ramp_gui = (self.clock(), r0, target, speed)
        else:
            self.jogging = True
        self.add_event(f"JOG {'+' if direction > 0 else '-'}{dist * 1e3:.4g} mm at {speed * 1e3:.4g} mm/s")
        self.update_status_display()

    def stop_manual_control(self):
        self.stop_controller()

    # ---- data ----
    def reset_data(self):
        self.data.clear()
        self.events.clear()
        self.results, self.permanent_set = {}, None
        self.d0, self.t_origin = 0.0, None
        self.model, self.design, self.sent_K = PlantModel(), stiffness_free_design(), None
        self.model_hist.clear()
        self.pole_hist.clear()
        self.last_loading = ("elastic", None)
        self.overlay, self.fracture_handled, self.schedule_pending = None, False, False
        self.mode_combo.setCurrentIndex(MODE_DISP)
        self.mode_combo.model().item(MODE_FORCE).setEnabled(False)
        self.send("RESET")
        self.running = self.jogging = False
        self.test, self.r_gui, self.ramp_gui = None, 0.0, None
        self.test_label.setText("Test: idle")
        for c in (self.c_d, self.c_r, self.c_F, self.c_F_model, self.c_F_drive, self.c_w, self.c_u, self.c_I,
                  self.c_ss, self.c_ss_true, self.c_ss_fit, self.c_ss_off, self.c_fd, self.c_fd_fit,
                  self.c_et, self.c_hold, self.c_hold_fit, self.c_ki, self.c_kr):
            c.setData([], [])
        self.c_ss_pts.setData([])
        self.c_pole_hist.setData([])
        self.c_pole_now.setData([])
        self.results_label.setText("Not measured yet.")
        self.update_status_display()
        self.update_dynamics_view()

    def save_data_to_file(self, path=None):
        interactive = path is None
        if interactive:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Data", "",
                                                            "CSV Files (*.csv);;Text Files (*.txt)")
            if not path:
                return
        try:
            a = self.data.arrays()
            A_cs, L0 = self.sample_geometry()
            eps, sig, eps_t, sig_t = strain_stress(a["d"], a["F"], A_cs, L0, self.d0)
            with open(path, "w") as f:
                f.write("# Straining stage measurement, load cell FUTEK LSB200 + HX711\n")
                f.write("# Sample Parameters:\n")
                f.write(f"# Cross-sectional Area (A_cs): {A_cs:.6e} m^2\n")
                f.write(f"# Length (L0): {L0:.6e} m\n")
                f.write(f"# Strain origin d0 (preload): {self.d0:.9e} m\n")
                f.write("# Calculated properties:\n")
                for line in self.results_text_lines():
                    f.write(f"#   {line}\n")
                f.write("\n# Measurement Events:\n")
                for event_time, name in self.events:
                    f.write(f"# {name} at t={event_time:.6f}s\n")
                f.write("\nTime(s),Displacement(m),Elongation(m),Strain(%),Force(N),Stress(Pa),"
                        "TrueStrain(-),TrueStress(Pa),Omega(rad/s),Current_est(A),Voltage(V),"
                        "Reference(m or N),Mode(0=disp,1=force,-1=off)\n")
                M = np.column_stack([a["t"], a["d"], a["d"] - self.d0, eps * 100, a["F"], sig, eps_t, sig_t,
                                     a["omega"], a["I"], a["u"], a["r"], a["mode"]])
                np.savetxt(f, M, delimiter=",", fmt="%.9g")
            self.log(f"Data saved to: {path}")
            if interactive:
                QtWidgets.QMessageBox.information(self, "Success", f"Data saved to:\n{path}")
        except Exception as e:
            self.log(f"Error saving file: {e}")
            if interactive:
                QtWidgets.QMessageBox.critical(self, "Error", f"Failed to save data:\n{e}")

    # ---- serial ----
    def on_tick(self):
        self.read_serial()
        self.update_test_program()
        self.update_signal_plots()

    def read_serial(self):
        try:
            if self.sim is not None:              # simulated stage: complete lines with their own time
                count = 0
                while self.link.in_waiting and count < 20000:
                    line = self.link.readline().decode(errors="replace").strip()
                    count += 1
                    if line:
                        self.handle_line(line, self.link.last_line_time)
                return
            # real port: read what is there without blocking, keep a partial line for later
            n = self.link.in_waiting
            if n:
                self._rx += self.link.read(n)
            if b"\n" not in self._rx:
                return
            *lines, self._rx = self._rx.split(b"\n")
            lines = [ln.decode(errors="replace").strip() for ln in lines]
            lines = [ln for ln in lines if ln]
            # data lines arrive at the HX711 rate: spread their time stamps over the batch
            now = self.clock()
            n_data = sum(1 for ln in lines if ln[0] in "-0123456789.")
            k = 0
            for ln in lines:
                if ln[0] in "-0123456789.":
                    t = now - (n_data - 1 - k) / HX711_RATE_HZ
                    t = max(t, self._t_last + 1e-6)
                    self._t_last = t
                    k += 1
                else:
                    t = now
                self.handle_line(ln, t)
        except Exception as e:
            self.log(f"[serial error] {e}")

    def handle_line(self, line, t):
        parts = line.split(",")
        vals = None
        if len(parts) in (2, 5):
            try:
                vals = [float(v) for v in parts]
            except ValueError:
                vals = None
        if vals is None:
            self.handle_message(line)
            return
        if len(vals) == 5:
            omega, I, d, F, u = vals
            if self.protocol != "v2":
                self.protocol = "v2"
                self.fw_label.setText("Firmware protocol: v2 (omega,I,d,F,u)")
        else:
            omega, d = vals
            I = F = u = float("nan")
            if self.protocol is None:
                self.protocol = "legacy"
                self.fw_label.setText("Firmware protocol: LEGACY (omega,d), no force: update the sketch")
        self.data.append((t, omega, I, d, F, u, self.current_reference(t),
                          self.ctrl_mode if self.running else -1))

    def handle_message(self, line):
        self.log(f"[Arduino] {line}")
        if "Pre-tensioning complete" in line:
            self.running = False
            self.set_strain_origin()
            self.add_event("PRETENSION complete")
            self.update_status_display()
        elif line.startswith("RAMP done"):
            self.on_ramp_done()
        elif line.startswith("JOG done"):
            self.jogging = False
            if self.running and self.ramp_gui is not None:
                self.r_gui, self.ramp_gui = self.ramp_gui[2], None
            self.update_status_display()
        elif "FORCE LIMIT" in line or "DISPLACEMENT LIMIT" in line:
            self.running = self.jogging = False
            self.test, self.ramp_gui = None, None
            self.add_event(line)
            self.test_label.setText("Test: aborted by the " + ("force" if "FORCE" in line else "displacement")
                                    + " limit")
            self.update_status_display()
        elif "protocol v2" in line:
            self.protocol = "v2"
            self.fw_label.setText("Firmware protocol: v2 (omega,I,d,F,u)")

    # ---- analysis (every 0.5 s, everything calculated from d, F, A_cs, L0) ----
    def fit_window(self):
        lo, hi = self.fit_start_input.value(), self.fit_end_input.value()
        if hi <= lo:
            return None
        L0_mm = self.length_input.value()
        return (lo / L0_mm, hi / L0_mm)

    def live_tangent_stiffness(self, t, d, F, load):
        sel = load & (t >= t[-1] - 3.0)
        if sel.sum() < 10 or np.ptp(d[sel]) < 5 * self.D_TOL:
            return None
        return float(np.polyfit(d[sel], F[sel], 1)[0])

    def _loaded_before(self, t, d, i0, quantum):
        """True when the displacement rose over the 2 s before index i0 (loading, not unloading)."""
        j = int(np.searchsorted(t, t[i0] - 2.0))
        return i0 > 0 and d[i0] - d[min(j, i0 - 1)] > max(2e-6, 5.0 * quantum)

    def analyse_holds(self, t, d, F, mode, labels, rate, rate_tol, quantum, A_cs, L0):
        loaded = F > 0.5 * self.preload_input.value()
        # stress relaxation: a displacement hold that follows loading (a hold after unloading
        # shows force recovery instead, which is not fitted as relaxation)
        runs = [("relaxation", (i0, i1)) for i0, i1 in
                detect_holds(t, rate, rate_tol, self.HOLD_MIN_S, mask=loaded & (mode != MODE_FORCE))
                if self._loaded_before(t, d, i0, quantum)]
        F_rate = windowed_rate(t, F, self.RATE_WINDOW)
        runs += [("creep", r) for r in detect_holds(t, F_rate, self.CREEP_RATE_TOL, self.HOLD_MIN_S,
                                                    mask=loaded & (mode == MODE_FORCE))]
        if not runs:
            return self.results.get("relax")
        kind, (i0, i1) = max(runs, key=lambda c: c[1][1])
        if kind == "relaxation":
            fit = fit_hold(t[i0:i1], F[i0:i1])
            if fit is None:
                return self.results.get("relax")
            # start of the loading ramp before the hold (skip the last second of settling)
            j = int(np.searchsorted(t, t[i0] - 1.0))
            while j > 0 and labels[j - 1] in (LOAD, RELOAD):
                j -= 1
            t_load = t[i0] - t[j] if j < i0 else 0.0
            v_load = (d[i0] - d[j]) / t_load if t_load > 0 else None
            res = relaxation_results(fit, float(np.mean(d[i0:i1])), self.d0, A_cs, L0, v_load, t_load,
                                     t_ref=t[i0] - self.RATE_WINDOW / 2)
            y_seg = F[i0:i1]
        else:
            fit = fit_hold(t[i0:i1], d[i0:i1])
            if fit is None:
                return self.results.get("relax")
            res = dict(d_start=fit["y0"], d_inf=fit["y_inf"], dd=fit["y_inf"] - fit["y0"], tau=fit["tau"],
                       r2=fit["r2"], significant=fit["significant"], F_hold=float(np.mean(F[i0:i1])),
                       duration=fit["duration"], t_start=fit["t_start"])
            y_seg = d[i0:i1]
        res.update(kind=kind, fit=fit, t_seg=t[i0:i1], y_seg=y_seg)
        return res

    def run_analysis(self):
        if len(self.data) < 20 or self.protocol != "v2":
            return
        a = self.data.arrays()
        A_cs, L0 = self.sample_geometry()
        sel = np.isfinite(a["F"])
        if self.t_origin is not None:
            sel &= a["t"] >= self.t_origin
        if sel.sum() < 20:
            return
        t, d, F, mode = a["t"][sel], a["d"][sel], a["F"][sel], a["mode"][sel]
        eps, sig, eps_t, sig_t = strain_stress(d, F, A_cs, L0, self.d0)
        quantum = displacement_quantum(d)                 # about one encoder count
        rate_tol = max(self.RATE_TOL, 4.0 * quantum / self.RATE_WINDOW)
        labels, rate = segment_branches(t, d, rate_tol, self.D_TOL, self.RATE_WINDOW)
        load = labels == LOAD
        res = {}
        fit = fit_elastic(eps, sig, d, F, load, self.fit_window())
        yld = offset_yield(eps, sig, load, fit, self.offset_input.value() * 1e-2) if fit else None
        ult = ultimate(eps, sig)
        res["fit"] = fit
        res["yield"] = yld
        res["prop"] = proportional_limit(eps, sig, load, fit) if fit else None
        res["ult"] = ult
        res["tangent"] = tangent_modulus(eps, sig, load)
        res["hard"] = hardening(eps, sig, load, fit, yld, ult)
        res["hollomon"] = hollomon_fit(eps, sig, load, yld, ult)
        res["energy"] = energies(eps, sig, load, fit, yld, ult)
        res["relax"] = self.analyse_holds(t, d, F, mode, labels, rate, rate_tol, quantum, A_cs, L0)
        res["k_tan_now"] = self.live_tangent_stiffness(t, d, F, load)
        res["fractured"] = bool(ult and ult["broken"])
        res["regime"] = classify_regime(t, d, F, res, eps[-1], 0.5 * self.preload_input.value(),
                                        rate_tol, self.D_TOL, window=self.RATE_WINDOW)
        res["eps_now"], res["sig_now"], res["F_now"] = float(eps[-1]), float(sig[-1]), float(F[-1])
        self.results = res
        self.update_model(res)
        self.maybe_reschedule()
        self.results_label.setText("\n".join(self.results_text_lines(res)))
        self.update_material_plots(d, F, eps, sig, eps_t, sig_t, res)
        self.update_overlay(a, res, A_cs, L0)
        self.update_dynamics_view()
        self.update_status_display()
        if res["fractured"] and not self.fracture_handled:
            self.fracture_handled = True
            self.add_event("FRACTURE detected")
            self.log("Fracture detected: motor stopped")
            if self.running or self.jogging:
                self.stop_controller()
            self.test_label.setText("Test: sample fractured")

    def update_model(self, res):
        """Measured plant model: k_e from the F-d fit, k_tan on the plastic envelope,
        k_v and tau from the last significant relaxation fit."""
        fit = res.get("fit")
        k_e = fit["k"] if fit and fit["k"] > 0 else None
        regime = res["regime"]
        if regime in ("elastic", "plastic", "softening"):
            self.last_loading = (regime, res.get("k_tan_now"))
        reg, k_tan = self.last_loading
        relax = res.get("relax")
        k_v, tau = 0.0, None
        if relax and relax.get("kind") == "relaxation" and relax.get("significant") and relax.get("k_v"):
            k_v, tau = max(relax["k_v"], 0.0), relax["tau"]
        new = regime_parameters(reg, k_e, k_tan, k_v, tau)
        old = self.model
        changed = (new.measured != old.measured or new.softening != old.softening
                   or _rel_change(new.k_i, old.k_i) > 0.25 or _rel_change(new.k_r, old.k_r) > 0.25
                   or _rel_change(new.tau, old.tau) > 0.25)
        new.regime = regime
        if not changed:
            self.model.regime = regime
            return
        self.model = new
        self.model_hist.append((self.clock(), new.k_i, new.k_r, new.tau, regime))
        if new.measured:
            self.mode_combo.model().item(MODE_FORCE).setEnabled(True)
        mode = self.ctrl_mode if new.measured else MODE_DISP
        self.design = design_lqr(mode, new)
        self.schedule_pending = True

    def results_text_lines(self, res=None):
        res = self.results if res is None else res
        nm = "not measured yet"
        lines = []
        fit = res.get("fit")
        if fit:
            lines.append(f"E (Young's modulus) : {fit['E'] / 1e9:.4g} GPa  (R^2 = {fit['r2']:.4f}, {fit['method']})")
            lines.append(f"k (stiffness, F-d)  : {fit['k']:.4g} N/m = {fit['k'] * 1e-3:.4g} N/mm")
        else:
            lines += [f"E (Young's modulus) : {nm}", f"k (stiffness, F-d)  : {nm}"]
        yld = res.get("yield")
        off = self.offset_input.value()
        lines.append(f"Yield ({off:g} % offset) : {yld['sig'] / 1e6:.4g} MPa at {yld['eps'] * 100:.4g} % strain"
                     if yld else f"Yield ({off:g} % offset) : not reached")
        prop = res.get("prop")
        if prop:
            lines.append(f"Proportional limit  : {prop['sig'] / 1e6:.4g} MPa at {prop['eps'] * 100:.4g} %")
        ult = res.get("ult")
        if ult:
            lines.append(f"UTS                 : {ult['uts'] / 1e6:.4g} MPa at {ult['eps_u'] * 100:.4g} % (uniform elongation)")
            if ult.get("broken"):
                lines.append(f"Fracture            : at {ult['eps_break'] * 100:.4g} % strain")
        hard = res.get("hard")
        if hard:
            lines.append(f"Plastic regime      : E_t = {hard['E_t'] / 1e6:.4g} MPa, H = d(sigma)/d(eps_p) = {hard['H'] / 1e6:.4g} MPa")
        hol = res.get("hollomon")
        if hol:
            lines.append(f"Hollomon            : K = {hol['K'] / 1e6:.4g} MPa, n = {hol['n']:.3f} (R^2 = {hol['r2']:.3f})")
        if self.permanent_set is not None:
            lines.append(f"Permanent set       : {self.permanent_set * 100:.4g} % (after unloading)")
        en = res.get("energy") or {}
        if "resilience" in en:
            lines.append(f"Resilience          : {en['resilience'] / 1e3:.4g} kJ/m^3")
        if "toughness" in en:
            lines.append(f"Toughness           : {en['toughness'] / 1e6:.4g} MJ/m^3 (area under loading curve)")
        rx = res.get("relax")
        if rx and rx.get("kind") == "relaxation":
            lines.append(f"Relaxation (hold at {rx['eps_hold'] * 100:.4g} %, {rx['duration']:.1f} s):")
            if rx["significant"]:
                lines.append(f"  F0 = {rx['F0']:.5g} N, F_inf = {rx['F_inf']:.5g} N, relaxation force dF = {rx['dF']:.4g} N ({rx['pct']:.3g} %)")
                extra = f", E_v = {rx['E_v'] / 1e9:.4g} GPa" if "E_v" in rx else ""
                lines.append(f"  tau = {rx['tau']:.4g} s{extra}")
                if "E0" in rx:
                    lines.append(f"  E0 = {rx['E0'] / 1e9:.4g} GPa, E_inf = {rx['E_inf'] / 1e9:.4g} GPa (secant at hold)")
            else:
                lines.append(f"  no significant relaxation (dF = {rx['dF']:.2g} N)")
        elif rx and rx.get("kind") == "creep":
            lines.append(f"Creep (hold at F = {rx['F_hold']:.4g} N, {rx['duration']:.1f} s): "
                         + (f"dd = {rx['dd'] * 1e3:.4g} mm, tau = {rx['tau']:.4g} s" if rx["significant"]
                            else "no significant creep"))
        if "F_now" in res:
            lines.append(f"Now                 : F = {res['F_now']:.5g} N, sigma = {res['sig_now'] / 1e6:.4g} MPa, "
                         f"eps = {res['eps_now'] * 100:.4g} %, regime = {res.get('regime', '-')}")
        return lines

    # ---- plots ----
    def update_signal_plots(self):
        n = len(self.data)
        if n == 0:
            return
        a = self.data.arrays()
        s = _decimate(n, self.PLOT_POINTS)
        t = a["t"][s]
        mode = a["mode"][s]
        self.c_d.setData(t, (a["d"][s] - self.d0) * 1e3)
        self.c_r.setData(t, np.where(mode == MODE_DISP, (a["r"][s] - self.d0) * 1e3, np.nan))
        self.c_F.setData(t, a["F"][s])
        if self.overlay is not None:
            self.c_F_model.setData(*self.overlay)
        if self.drive_check.isChecked():
            _, _, tau_m = motor_quantities(a["omega"][s], a["I"][s], a["u"][s])
            self.c_F_drive.setData(t, sample_force(a["omega"][s], tau_m))
        else:
            self.c_F_drive.setData([], [])
        self.c_w.setData(t, a["omega"][s])
        self.c_u.setData(t, a["u"][s])
        self.c_I.setData(t, a["I"][s] * 1e3)
        A_cs, _ = self.sample_geometry()
        F = self.last_F()
        self.force_label.setText(f"F = {F:.5f} N   (sigma = {F / A_cs / 1e6:.4g} MPa)")

    def update_material_plots(self, d, F, eps, sig, eps_t, sig_t, res):
        _, L0 = self.sample_geometry()
        s = _decimate(len(eps), self.PLOT_POINTS)
        self.c_ss.setData(eps[s] * 100, sig[s] / 1e6)
        if self.true_curve_check.isChecked():
            self.c_ss_true.setData(eps_t[s] * 100, sig_t[s] / 1e6)
        else:
            self.c_ss_true.setData([], [])
        self.c_fd.setData((d[s] - self.d0) * 1e3, F[s])
        fit, yld, ult = res.get("fit"), res.get("yield"), res.get("ult")
        if fit:
            e_end = yld["eps"] * 1.15 if yld else min(float(np.nanmax(eps)), fit["eps_hi"] * 2.5)
            e_line = np.array([max(-fit["sigma0"] / fit["E"], 0.0), e_end])
            self.c_ss_fit.setData(e_line * 100, (fit["E"] * e_line + fit["sigma0"]) / 1e6)
            off = self.offset_input.value() * 1e-2
            e_off = np.array([off, e_end + off])
            self.c_ss_off.setData(e_off * 100, (fit["E"] * (e_off - off) + fit["sigma0"]) / 1e6)
            d_line = self.d0 + e_line * L0
            self.c_fd_fit.setData((d_line - self.d0) * 1e3, fit["k"] * d_line + fit["F0"])
        else:
            for c in (self.c_ss_fit, self.c_ss_off, self.c_fd_fit):
                c.setData([], [])
        spots = []
        if yld:
            spots.append(dict(pos=(yld["eps"] * 100, yld["sig"] / 1e6), symbol='star', size=16, brush='y'))
        if ult:
            spots.append(dict(pos=(ult["eps_u"] * 100, ult["uts"] / 1e6), symbol='d', size=13, brush='c'))
            if ult.get("broken"):
                spots.append(dict(pos=(ult["eps_break"] * 100, ult["sig_break"] / 1e6), symbol='x', size=14, brush='r'))
        self.c_ss_pts.setData(spots)
        tm = res.get("tangent")
        if tm is not None:
            self.c_et.setData(tm[0] * 100, tm[1] / 1e9)
        rx = res.get("relax")
        if rx and "t_seg" in rx:
            fitp = rx["fit"]
            tt = rx["t_seg"] - rx["t_seg"][0]
            scale = 1.0 if rx["kind"] == "relaxation" else 1e3
            offset = 0.0 if rx["kind"] == "relaxation" else self.d0
            ss = _decimate(len(tt), 2000)
            self.c_hold.setData(tt[ss], (rx["y_seg"][ss] - offset) * scale)
            self.c_hold_fit.setData(tt[ss], (_exp_model(tt[ss], fitp["y_inf"], fitp["delta"], fitp["tau"]) - offset) * scale)
            self.p_hold.setLabel('left', 'F (N)' if rx["kind"] == "relaxation" else 'd (mm)')
            self.p_hold.setTitle("Last hold: stress relaxation (F at constant d)" if rx["kind"] == "relaxation"
                                 else "Last hold: creep (d at constant F)")

    def update_overlay(self, a, res, A_cs, L0):
        """Force predicted by the identified model, driven by the measured d (every 2 s)."""
        if time.time() - self.last_overlay_wall < 2.0:
            return
        self.last_overlay_wall = time.time()
        P = identified_material(res, A_cs, L0, self.d0)
        if P is None:
            self.overlay = None
            self.c_F_model.setData([], [])
            return
        sel = np.isfinite(a["F"])
        if self.t_origin is not None:
            sel &= a["t"] >= self.t_origin
        t, d = a["t"][sel], a["d"][sel]
        s = _decimate(len(t), 3000)
        self.overlay = (t[s], simulate_material_force(t[s], d[s], P))

    def update_dynamics_view(self):
        A_cs, L0 = self.sample_geometry()
        text = describe_design(self.design, A_cs, L0, self.sent_K)
        if text != self._dyn_text_cache:
            self._dyn_text_cache = text
            self.dyn_text.setPlainText(text)
        if self.model_hist:
            th = [h[0] for h in self.model_hist]
            self.c_ki.setData(th, [h[1] for h in self.model_hist])
            self.c_kr.setData(th, [h[2] for h in self.model_hist])
        hist = [dict(pos=(z.real, z.imag)) for _, poles in self.pole_hist[:-1] for z in poles]
        self.c_pole_hist.setData(hist)
        now = dominant_poles(self.design.poles_closed, 100.0)
        self.c_pole_now.setData([dict(pos=(z.real, z.imag)) for z in now])


# -----------------------
# Demo (--demo): scripted virtual test through the real button handlers
# -----------------------
class DemoScript:
    """jog -> tare -> preload -> ramp to 0.5 mm + 30 s hold -> unload -> ramp to fracture.
    Prints the identified properties next to the simulator ground truth."""

    def __init__(self, gui, shot_dir=None):
        self.gui, self.dir = gui, shot_dir
        self.steps = [self._jog, self._wait_jog, self._tare, self._preload, self._wait_preload,
                      self._ramp_hold, self._wait_test, self._unload, self._wait_test,
                      self._fracture, self._wait_fracture, self._report]
        self.i, self.t_mark = 0, 0.0
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self._tick)

    def start(self):
        self.t_mark = self.gui.clock()
        self.timer.start(100)

    def _elapsed(self):
        return self.gui.clock() - self.t_mark

    def _tick(self):
        if self.i < len(self.steps) and self.steps[self.i]():
            self.i += 1
            self.t_mark = self.gui.clock()

    def _jog(self):
        if self._elapsed() < 1.0:
            return False
        g = self.gui
        g.jog_distance_input.setValue(0.005)
        g.jog_speed_input.setValue(0.05)
        g.jog(+1)
        return True

    def _wait_jog(self):
        return (not self.gui.jogging and self._elapsed() > 0.5) or self._elapsed() > 10.0

    def _tare(self):
        if self._elapsed() < 0.5:
            return False
        self.gui.tare()
        return True

    def _preload(self):
        if self._elapsed() < 1.0:
            return False
        self.gui.pre_tension()
        return True

    def _wait_preload(self):
        return self.gui.t_origin is not None or self._elapsed() > 60.0

    def _ramp_hold(self):
        if self._elapsed() < 1.0:
            return False
        g = self.gui
        g.target_input.setValue(0.5)
        g.rate_input.setValue(0.01)
        g.hold_input.setValue(30.0)
        g.start_ramp_hold()
        return True

    def _wait_test(self):
        return self.gui.test is None and self._elapsed() > 1.0

    def _unload(self):
        self.gui.rate_input.setValue(0.02)
        self.gui.start_unload()
        return True

    def _fracture(self):
        if self._elapsed() < 2.0:
            return False
        g = self.gui
        g.target_input.setValue(2.0)
        g.rate_input.setValue(0.02)
        g.hold_input.setValue(0.0)
        g.start_ramp_hold()
        return True

    def _wait_fracture(self):
        return (self.gui.results.get("fractured", False)
                or (self.gui.test is None and self._elapsed() > 2.0) or self._elapsed() > 300.0)

    def _report(self):
        g = self.gui
        g.run_analysis()
        P = g.sim.P if g.sim is not None else MaterialParams()
        strain_rate = 0.01e-3 / P.L0
        dF_expected = P.E_v * P.tau_v * strain_rate * P.A_cs * (1.0 - math.exp(-50.0 / P.tau_v))
        print("\n=== Demo: properties calculated from d, F, A_cs and L0 ===")
        for line in g.results_text_lines():
            print("  " + line)
        print("=== Simulator ground truth (unknown to the analysis) ===")
        print(f"  backbone E = {P.E / 1e9:.3g} GPa (k = {P.E * P.A_cs / P.L0:.4g} N/m), Maxwell arm E_v = "
              f"{P.E_v / 1e9:.3g} GPa, tau_v = {P.tau_v:.3g} s (instantaneous E = {(P.E + P.E_v) / 1e9:.3g} GPa)")
        print(f"  sigma_y0 = {P.sigma_y0 / 1e6:.3g} MPa, Voce Q = {P.Q_voce / 1e6:.3g} MPa, "
              f"H = {P.H / 1e6:.3g} MPa, necking from {P.eps_neck * 100:.3g} %, fracture at {P.eps_fracture * 100:.3g} %")
        print(f"  expected relaxation after the 0.5 mm ramp at 0.01 mm/s: dF = {dF_expected:.4g} N, tau = {P.tau_v:.3g} s")
        print("  (a ramp measures an apparent modulus between E and E + E_v, because the material is viscoelastic)")
        if self.dir:
            os.makedirs(self.dir, exist_ok=True)
            for i, name in enumerate(("signals", "material", "dynamics")):
                g.tabs.setCurrentIndex(i)
                QtWidgets.QApplication.processEvents()
                QtWidgets.QApplication.processEvents()
                g.grab().save(os.path.join(self.dir, f"demo_{name}.png"))
            g.save_data_to_file(os.path.join(self.dir, "demo_data.csv"))
            print(f"Screenshots and CSV saved in {self.dir}")
        QtCore.QTimer.singleShot(300, QtWidgets.QApplication.quit)
        return True


# -----------------------
# Built-in self-test (python ControlInterface4.py --selftest)
# -----------------------
def _elastic_params(**kw):
    """Ground truth without viscosity, hardening saturation or necking (for exact checks)."""
    base = dict(E=1.0e9, sigma_y0=1e15, Q_voce=0.0, beta_voce=1.0, H=0.0, E_v=0.0, tau_v=1.0,
                eps_neck=10.0, eps_fracture=11.0, slack=0.0)
    base.update(kw)
    return MaterialParams(**base)


def _data_lines(lines):
    rows = []
    for ln in lines:
        parts = ln.split(",")
        if len(parts) == 5:
            try:
                rows.append([float(v) for v in parts])
            except ValueError:
                pass
    return np.array(rows)


class TestDynamics(unittest.TestCase):
    def test_reachability_determinant(self):
        for ki, kr, tau in [(1500.0, 1000.0, 5.0), (1000.0, 50.0, 2.0)]:
            A, B, _, _ = plant_matrices(ki, kr, tau)
            W = np.hstack([np.linalg.matrix_power(A, i) @ B for i in range(4)])
            analytic = reachability(A, B, PlantModel(ki, kr, tau, "x", True))["det"]
            self.assertAlmostEqual(np.linalg.det(W) / analytic, 1.0, places=6)

    def test_elastic_mode_uncontrollable_but_stable(self):
        m = PlantModel(1000.0, 1000.0, 2.0, "elastic", True)
        A, B, _, _ = plant_matrices(m.k_i, m.k_r, m.tau)
        r = reachability(A, B, m)
        self.assertEqual(r["rank"], 3)
        self.assertAlmostEqual(r["uncontrollable"][0].real, -1.0 / m.tau, places=6)
        self.assertTrue(r["stabilizable"])
        m2 = PlantModel(1500.0, 1000.0, 2.0, "viscoelastic", True)
        self.assertEqual(reachability(*plant_matrices(m2.k_i, m2.k_r, m2.tau)[:2], m2)["rank"], 4)

    def test_lqr_stable_in_all_regimes(self):
        models = [PlantModel(1e3, 1e3, 1.0, "elastic", True), PlantModel(1.5e3, 1e3, 5.0, "viscoelastic", True),
                  PlantModel(1e3, 50.0, 2.0, "plastic", True), PlantModel(1e3, 10.0, 2.0, "plastic", True)]
        for m in models:
            for mode in (MODE_DISP, MODE_FORCE):
                d = design_lqr(mode, m)
                self.assertLess(d.poles_closed.real.max(), 0.0)
                if mode == MODE_FORCE:
                    slow = dominant_poles(d.poles_closed, 100.0)
                    self.assertTrue(np.any(np.abs(slow.real + 12.0) < 2.0), msg=str(slow))

    def test_stiffness_free_design_is_robust(self):
        K0 = stiffness_free_design().K
        for ki, kr, tau in [(10.0, 10.0, 1.0), (1e3, 1e3, 1.0), (1.5e3, 1e3, 5.0), (1e3, 50.0, 2.0), (1e6, 1e6, 1.0)]:
            A, B, _, _ = plant_matrices(ki, kr, tau)
            Aa = np.block([[A, np.zeros((4, 1))], [-C_TRACK[MODE_DISP], np.zeros((1, 1))]])
            Ba = np.vstack((B, [[0.0]]))
            self.assertLess(np.linalg.eigvals(Aa - Ba @ K0[None, :]).real.max(), 0.0)

    def test_force_mode_needs_measured_stiffness(self):
        with self.assertRaises(ValueError):
            design_lqr(MODE_FORCE, PlantModel())

    def test_observable_from_encoder_and_load_cell(self):
        for m in (PlantModel(1e3, 1e3, 1.0, "elastic", True), PlantModel(1e3, 50.0, 2.0, "plastic", True)):
            A = plant_matrices(m.k_i, m.k_r, m.tau)[0]
            self.assertEqual(observability(A, C_MEASURED, m)["rank"], 4)


class TestMaterialAnalysis(unittest.TestCase):
    A_CS, L0 = 0.01e-6, 10e-3

    def test_elastic_yield_hardening_stiffness(self):
        E, sy, H = 2e9, 40e6, 200e6
        P = _elastic_params(E=E, sigma_y0=sy, H=H, A_cs=self.A_CS, L0=self.L0)
        t = np.linspace(0.0, 60.0, 3000)
        d = np.linspace(0.0, 0.06 * self.L0, 3000)
        F = simulate_material_force(t, d, P)
        eps, sig, _, _ = strain_stress(d, F, self.A_CS, self.L0)
        load = segment_branches(t, d, 5e-7, 1e-6)[0] == LOAD
        fit = fit_elastic(eps, sig, d, F, load)
        self.assertAlmostEqual(fit["E"] / E, 1.0, delta=0.02)
        self.assertAlmostEqual(fit["k"] / (E * self.A_CS / self.L0), 1.0, delta=0.02)
        self.assertAlmostEqual(fit["k"] * self.L0 / self.A_CS / fit["E"], 1.0, delta=1e-6)
        yld = offset_yield(eps, sig, load, fit, 0.002)
        self.assertAlmostEqual(yld["sig"] / (sy + H * 0.002), 1.0, delta=0.02)
        hard = hardening(eps, sig, load, fit, yld, ultimate(eps, sig))
        self.assertAlmostEqual(hard["H"] / H, 1.0, delta=0.02)
        self.assertAlmostEqual(hard["E_t"] / (E * H / (E + H)), 1.0, delta=0.02)

    def test_relaxation_fit(self):
        rng = np.random.default_rng(3)
        t = np.linspace(0.0, 40.0, 3200)
        F = 0.30 + 0.05 * np.exp(-t / 6.0) + rng.normal(0.0, 1e-4, t.size)
        fit = fit_hold(t, F)
        self.assertAlmostEqual(fit["tau"] / 6.0, 1.0, delta=0.02)
        res = relaxation_results(fit, d_hold=5e-4, d0=0.0, A_cs=self.A_CS, L0=self.L0)
        self.assertAlmostEqual(res["dF"] / 0.05, 1.0, delta=0.02)
        self.assertTrue(res["significant"])

    def test_regime_parameters(self):
        self.assertFalse(regime_parameters("elastic").measured)
        m = regime_parameters("plastic", k_e=1000.0, k_tan=50.0, k_v=300.0, tau_v=8.0)
        self.assertEqual((m.k_i, m.k_r, m.tau), (350.0, 50.0, 8.0))


class TestClosedLoop(unittest.TestCase):
    """Simulated stage (same protocol and control law as the sketch), in virtual time."""

    def _sim(self, params=None):
        sim = SimulatedArduino(params=params or MaterialParams(), realtime=False)
        sim.advance(0.2)
        sim.drain()
        return sim

    @staticmethod
    def _cmd(sim, text, duration):
        sim.write((text + "\n").encode())
        sim.advance(duration)
        return sim.drain()

    def test_displacement_step(self):
        sim = self._sim(_elastic_params(slack=1e-3))
        self._cmd(sim, "START", 0.5)
        d = _data_lines(self._cmd(sim, "R:0.0001", 10.0))[-40:, 2]
        self.assertLess(abs(d.mean() - 1e-4), 3e-7)

    def test_force_step_with_measured_stiffness(self):
        sim = self._sim(_elastic_params(slack=0.0))
        k = sim.P.E * sim.P.A_cs / sim.P.L0
        self._cmd(sim, "TARE", 0.5)
        design = design_lqr(MODE_FORCE, PlantModel(k, k, 1.0, "elastic", True))
        self._cmd(sim, format_gain_command(0.2, design.K, MODE_FORCE).strip(), 0.1)
        F = _data_lines(self._cmd(sim, "START", 15.0))[-80:, 3]
        self.assertLess(abs(F.mean() - 0.2), 2e-3)

    def test_preload_stops_at_force(self):
        sim = self._sim()
        self._cmd(sim, "TARE", 0.5)
        lines = self._cmd(sim, "PRELOAD:0.01", 20.0)
        self.assertIn("Pre-tensioning complete, controller stopped", lines)
        self.assertEqual(sim.state, "idle")
        self.assertTrue(0.005 <= sim.F_true < 0.015, msg=str(sim.F_true))

    def test_force_limit_cuts_motor(self):
        sim = self._sim()
        self._cmd(sim, "TARE", 0.5)
        self._cmd(sim, "FLIM:0.05", 0.1)
        lines = self._cmd(sim, "JOG:1,0.001,0.0001", 15.0)
        self.assertTrue(any(ln.startswith("FORCE LIMIT") for ln in lines))
        self.assertEqual(sim.state, "idle")
        self.assertLess(sim.F_true, 0.06)

    def test_jog_distance(self):
        sim = self._sim(_elastic_params(slack=1e-2))
        d_start = sim._d_meas()
        lines = self._cmd(sim, "JOG:1,0.0001,0.00005", 8.0)
        self.assertIn("JOG done", lines)
        self.assertLess(abs(sim._d_meas() - d_start - 1e-4), 3e-7)
        self.assertEqual(sim.state, "idle")

    def test_ramp_rate(self):
        sim = self._sim(_elastic_params(slack=1e-2))
        self._cmd(sim, "START", 0.5)
        lines = self._cmd(sim, "RAMP:0.0002,0.00001", 25.0)
        self.assertIn("RAMP done", lines)
        data = _data_lines(lines)
        t = np.arange(len(data)) / HX711_RATE_HZ
        mid = (t > 5.0) & (t < 15.0)
        rate = np.polyfit(t[mid], data[mid, 2], 1)[0]
        self.assertAlmostEqual(rate / 1e-5, 1.0, delta=0.02)
        self.assertLess(abs(data[-1, 2] - 2e-4), 3e-7)


def run_selftest():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    for case in (TestDynamics, TestMaterialAnalysis, TestClosedLoop):
        suite.addTests(loader.loadTestsFromTestCase(case))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


# -----------------------
# Start app
# -----------------------
def open_link(args):
    """Serial port to the Arduino, or the simulated stage."""
    if args.sim or args.demo:
        sim = SimulatedArduino(speed=args.sim_speed)
        return sim, sim, f"simulated stage (x{args.sim_speed:g} speed)"
    try:
        import serial
        ser = serial.Serial(args.port, args.baud, timeout=1)
        time.sleep(2)                 # the Uno resets when the port opens
        return ser, None, f"{args.port} @ {args.baud} baud"
    except Exception as e:
        answer = QtWidgets.QMessageBox.question(
            None, "Serial port", f"Could not open {args.port}:\n{e}\n\nStart with the simulated stage instead?")
        if answer == QtWidgets.QMessageBox.Yes:
            sim = SimulatedArduino(speed=args.sim_speed)
            return sim, sim, "simulated stage (serial port unavailable)"
        return None, None, ""


def main(argv=None):
    ap = argparse.ArgumentParser(description="Straining-stage control GUI with load-cell force feedback")
    ap.add_argument("--port", default=SERIAL_PORT, help=f"serial port (default {SERIAL_PORT})")
    ap.add_argument("--baud", type=int, default=BAUD_RATE)
    ap.add_argument("--sim", action="store_true", help="use the simulated stage (no hardware)")
    ap.add_argument("--sim-speed", type=float, default=1.0, help="simulation speed factor")
    ap.add_argument("--demo", action="store_true", help="scripted test in the simulator")
    ap.add_argument("--screenshot", default=None, help="folder for demo screenshots and CSV")
    ap.add_argument("--selftest", action="store_true", help="run the built-in checks and exit")
    args = ap.parse_args(argv)
    if args.selftest:
        return run_selftest()
    if not QT_AVAILABLE:
        print("PyQt5 and pyqtgraph are needed: pip install PyQt5 pyqtgraph pyserial")
        return 1
    if args.demo and args.sim_speed == 1.0:
        args.sim_speed = 20.0
    app = QtWidgets.QApplication(sys.argv)
    link, sim, name = open_link(args)
    if link is None:
        return 1
    gui = StageGUI(link, sim, name)
    gui.show()
    if args.demo:
        gui.demo = DemoScript(gui, args.screenshot)
        QtCore.QTimer.singleShot(500, gui.demo.start)
    code = app.exec_()
    try:
        link.close()
    except Exception:
        pass
    return code


if __name__ == "__main__":
    sys.exit(main())

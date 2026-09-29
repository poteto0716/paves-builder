"""Fits and aggregation for the thermomechanical workflow."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np


def smooth_hyperbola(temperature, rho0, tg, slope_glass, slope_change, width):
    """Continuous two-asymptote density model.

    The low/high-temperature asymptotic slopes are ``slope_glass`` and
    ``slope_glass + slope_change``.  ``tg`` is the centre/intersection of the
    rounded transition and ``width`` controls only its smoothness.
    """
    x = np.asarray(temperature, float) - tg
    return rho0 + slope_glass * x + 0.5 * slope_change * (x + np.sqrt(x * x + width * width))


def bin_curve(temperature, values, width_K=10.0):
    """Average a cooling trace in target-temperature bins."""
    t = np.asarray(temperature, float)
    y = np.asarray(values, float)
    good = np.isfinite(t) & np.isfinite(y)
    t, y = t[good], y[good]
    if len(t) < 5:
        raise ValueError('at least five finite cooling samples are required')
    origin = math.floor(float(t.min()) / width_K) * width_K
    index = np.floor((t - origin) / width_K + 1e-12).astype(int)
    out = []
    for i in np.unique(index):
        take = index == i
        out.append((float(np.mean(t[take])), float(np.mean(y[take])), int(take.sum())))
    return np.asarray(out, float)


def fit_tg(temperature, density, bin_width_K=10.0):
    """Fit density versus temperature and return Tg and linear CTEs.

    Density normally decreases with temperature and its high-temperature
    slope is more negative.  The linear CTE follows directly from the fitted
    density slope: alpha_L = -(1/(3 rho)) d(rho)/dT.
    """
    from scipy.optimize import least_squares

    binned = bin_curve(temperature, density, bin_width_K)
    t, rho = binned[:, 0], binned[:, 1]
    if len(t) < 8:
        raise ValueError('at least eight temperature bins are required for the Tg fit')
    order = np.argsort(t)
    t, rho = t[order], rho[order]
    span = float(np.ptp(t))
    edge = max(3, len(t) // 5)
    low_fit = np.polyfit(t[:edge], rho[:edge], 1)
    high_fit = np.polyfit(t[-edge:], rho[-edge:], 1)
    m_low = float(low_fit[0])
    m_high = float(high_fit[0])
    # The high-temperature density slope is expected to be more negative.
    dm = min(-1e-10, m_high - m_low)
    tg0 = float(np.median(t))
    rho0 = float(np.interp(tg0, t, rho))
    scale = max(float(np.std(rho)), 1e-8)
    slope_limit = max(0.05, 20 * scale / span)
    lower = [float(rho.min() - 5 * scale), float(t.min() + 0.03 * span),
             -slope_limit, -2 * slope_limit, 0.25]
    upper = [float(rho.max() + 5 * scale), float(t.max() - 0.03 * span),
             slope_limit, 0.0, span]
    x0 = [rho0, tg0, np.clip(m_low, lower[2], upper[2]),
          np.clip(dm, lower[3], -1e-10), min(20.0, span / 3)]
    result = least_squares(
        lambda p: smooth_hyperbola(t, *p) - rho,
        x0,
        bounds=(lower, upper),
        max_nfev=100000,
    )
    if not result.success:
        raise RuntimeError(f'hyperbolic Tg fit failed: {result.message}')
    rho0, tg, mg, dmh, width = map(float, result.x)
    pred = smooth_hyperbola(t, *result.x)
    ss_res = float(np.sum((rho - pred) ** 2))
    ss_tot = float(np.sum((rho - np.mean(rho)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot else None
    rho_tg = float(smooth_hyperbola(np.asarray([tg]), *result.x)[0])
    alpha_glass = -mg / (3.0 * rho_tg)
    alpha_rubber = -(mg + dmh) / (3.0 * rho_tg)
    return {
        'tg_K': tg,
        'density_at_tg_g_cm3': rho_tg,
        'density_slope_glass_g_cm3_K': mg,
        'density_slope_rubber_g_cm3_K': mg + dmh,
        'alpha_linear_glass_1_K': alpha_glass,
        'alpha_linear_rubber_1_K': alpha_rubber,
        'transition_width_K': width,
        'r_squared': r2,
        'rmse_g_cm3': math.sqrt(ss_res / len(rho)),
        'bins': [
            {'temperature_K': float(tt), 'density_g_cm3': float(rr),
             'samples': int(nn), 'fit_density_g_cm3': float(pp)}
            for (tt, rr, nn), pp in zip(binned[order], pred)
        ],
    }


def fit_modulus(strain, stress_MPa, min_strain=0.0, max_strain=0.02):
    """OLS fit of tensile stress versus engineering strain."""
    e = np.asarray(strain, float)
    s = np.asarray(stress_MPa, float)
    use = np.isfinite(e) & np.isfinite(s) & (e >= min_strain) & (e <= max_strain)
    if use.sum() < 3:
        raise ValueError('at least three stress samples are required in the elastic fit range')
    slope, intercept = np.polyfit(e[use], s[use], 1)
    pred = slope * e[use] + intercept
    ss_res = float(np.sum((s[use] - pred) ** 2))
    ss_tot = float(np.sum((s[use] - np.mean(s[use])) ** 2))
    return {
        'tensile_modulus_MPa': float(slope),
        'tensile_modulus_GPa': float(slope / 1000.0),
        'intercept_MPa': float(intercept),
        'r_squared': 1.0 - ss_res / ss_tot if ss_tot else None,
        'rmse_MPa': math.sqrt(ss_res / int(use.sum())),
        'samples': int(use.sum()),
        'fit_min_strain': float(min_strain),
        'fit_max_strain': float(max_strain),
    }


def read_csv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh))


def statistics(values):
    x = np.asarray(values, float)
    finite = x[np.isfinite(x)]
    if not len(finite):
        return {'n': 0, 'mean': None, 'std': None, 'sem': None}
    std = float(np.std(finite, ddof=1)) if len(finite) > 1 else 0.0
    return {'n': int(len(finite)), 'mean': float(np.mean(finite)), 'std': std,
            'sem': std / math.sqrt(len(finite))}


def write_rows(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

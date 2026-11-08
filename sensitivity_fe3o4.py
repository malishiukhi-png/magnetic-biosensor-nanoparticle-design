"""Sensitivity of the Fe3O4 usable window to the effective anisotropy constant and to the chaining criterion,
and a simple estimate for multicore iron oxide particles -> sensitivity_fe3o4.json"""
import json, numpy as np
from model3 import *

T0 = 298.15; DELTA = 12e-9; LAYER = 5e-9


def lam(d, dh, Ms):
    m = Ms * np.pi * d**3 / 6
    return mu0 * m**2 / (4 * np.pi * dh**3 * kB * T0)


def window(K, lam_max, mat='Fe3O4', delta=DELTA, Ms=None):
    Ms = MAT[mat]['Ms'] if Ms is None else Ms
    d = np.linspace(4e-9, 40e-9, 7201)
    dh = d + 2 * delta
    tN, tB, _ = taus(d, dh, mat, T0, K=K)
    ok = (tN / (tN + tB) > 0.9) & (lam(d, dh, Ms) < lam_max)
    if not ok.any():
        return None
    return [round(float(d[ok][0] * 1e9), 1), round(float(d[ok][-1] * 1e9), 1)]   # 0.1 nm, as in analysis.py


def mass_frac(lo, hi, s=0.1):
    best = (0.0, None)
    for d0 in np.linspace(lo * 0.7, hi * 1.3, 600):
        x = np.linspace(np.log(d0) - 5 * s, np.log(d0) + 5 * s, 4001); d = np.exp(x)
        w = np.exp(-(x - np.log(d0))**2 / (2 * s * s)) * d**3; w /= w.sum()
        fr = w[(d >= lo) & (d <= hi)].sum()
        if fr > best[0]:
            best = (fr, d0)
    return round(100 * best[0], 1), round(best[1], 2)


out = {'grid': []}
for K in (8e3, 13e3, 20e3, 30e3):
    for lm in (0.5, 1.0, 2.0):
        w = window(K, lm)
        row = dict(K_kJ=K / 1e3, lam_max=lm, window=w)
        if w:
            row['width_nm'] = round(w[1] - w[0], 2)
            row['mass_frac_sigma0.1'], row['best_median'] = mass_frac(*w)
            f0, _ = peak(w[0] * 1e-9, DELTA, 'Fe3O4', K=K)
            f1, _ = peak(w[0] * 1e-9, DELTA, 'Fe3O4', K=K, cover=1.0, layer=LAYER)
            row['shift_full_at_lower_edge'] = round(100 * (1 - f1 / f0), 1)
        out['grid'].append(row)
# Thicker shells widen the Fe3O4 window (K = 13 kJ m-3, lambda < 1)
out['shell'] = {str(s): window(13e3, 1.0, delta=s * 1e-9) for s in (12, 15, 20, 25)}
# CoFe2O4 reference for the same criteria
out['CoFe2O4_lam'] = {str(lm): window(180e3, lm, mat='CoFe2O4') for lm in (0.5, 1.0, 2.0)}

# ---- multicore particles: blocked collective moment, Brownian relaxation of the whole particle -----------
# A multicore particle of hydrodynamic diameter Dh contains magnetic material with volume fraction phi_m of the
# magnetic cluster (diameter Dc = Dh - 2*shell); only a fraction q of the core moments is aligned (net moment).
mc = []
for Dh in (50e-9, 80e-9, 100e-9, 150e-9):
    for q in (0.2, 0.5):
        Dc = Dh - 2 * 10e-9; phi_m = 0.5
        m = q * phi_m * MAT['Fe3O4']['Ms'] * np.pi * Dc**3 / 6
        lam_mc = mu0 * m**2 / (4 * np.pi * Dh**3 * kB * T0)
        shift = 100 * (1 - (Dh / (Dh + 2 * LAYER))**3)
        tB = 3 * eta_water(T0) * np.pi * Dh**3 / 6 / (kB * T0)
        mc.append(dict(Dh_nm=Dh * 1e9, q=q, lam=round(float(lam_mc), 2), shift_full_pct=round(float(shift), 1),
                       fB_kHz=round(float(1 / (2 * np.pi * tB) / 1e3), 3)))
out['multicore'] = mc
json.dump(out, open('sensitivity_fe3o4.json', 'w'), indent=1)
print(json.dumps(out, indent=1))

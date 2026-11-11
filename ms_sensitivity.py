"""Revision: sensitivity of the usable windows to the saturation magnetization of the cores (Ms enters only the
dipolar coupling constant) -> ms_sensitivity.json"""
import json, sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model3 import MAT, taus, mu0, kB

T0 = 298.15; DELTA = 12e-9


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
    return [round(float(d[ok][0] * 1e9), 1), round(float(d[ok][-1] * 1e9), 1)]


out = {}
for Ms in (350e3, 390e3, 446e3):
    out['Fe3O4_Ms%d' % (Ms / 1e3)] = {str(K / 1e3): window(K, 1.0, Ms=Ms) for K in (10e3, 13e3, 20e3, 30e3)}
for Ms in (300e3, 425e3):
    out['CoFe2O4_Ms%d' % (Ms / 1e3)] = window(180e3, 1.0, mat='CoFe2O4', Ms=Ms)
print(json.dumps(out, indent=1))
json.dump(out, open('ms_sensitivity.json', 'w'), indent=1)


def mass_frac(lo, hi, s=0.1):
    best = (0.0, None)
    for d0 in np.linspace(lo * 0.7, hi * 1.3, 600):
        x = np.linspace(np.log(d0) - 5 * s, np.log(d0) + 5 * s, 4001); d = np.exp(x)
        w = np.exp(-(x - np.log(d0))**2 / (2 * s * s)) * d**3; w /= w.sum()
        fr = w[(d >= lo) & (d <= hi)].sum()
        if fr > best[0]:
            best = (fr, d0)
    return round(100 * best[0], 1), round(best[1], 2)


out['mass_frac_sigma0.1'] = {}
for key in ('Fe3O4_Ms350', 'Fe3O4_Ms390', 'Fe3O4_Ms446'):
    for K, w in out[key].items():
        if w:
            out['mass_frac_sigma0.1']['%s_K%s' % (key, K)] = mass_frac(*w)
# lowest anisotropy with a usable window, per Ms
out['K_min'] = {}
for Ms in (350e3, 390e3, 446e3):
    for K in np.arange(6e3, 14e3, 0.25e3):
        if window(K, 1.0, Ms=Ms):
            out['K_min']['%d' % (Ms / 1e3)] = float(K / 1e3); break
print(json.dumps({k: out[k] for k in ('mass_frac_sigma0.1', 'K_min')}, indent=1))
json.dump(out, open('ms_sensitivity.json', 'w'), indent=1)

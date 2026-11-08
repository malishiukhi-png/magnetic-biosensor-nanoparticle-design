"""Design analysis for Brownian-relaxation biosensing -> results.json"""
import json, numpy as np
from model3 import *
R = {}; T0 = 298.15; DELTA = 12e-9; LAYER = 5e-9
def lam(d, dh, mat, T=T0):
    """Dipolar coupling constant at contact: mu0 m^2 / (4 pi dh^3 kB T)."""
    m = MAT[mat]['Ms'] * np.pi * d ** 3 / 6
    return mu0 * m ** 2 / (4 * np.pi * dh ** 3 * kB * T)
d = np.linspace(5e-9, 30e-9, 501)
win = {}
for mat in MAT:
    dh = d + 2 * DELTA; tN, tB, t = taus(d, dh, mat, T0)
    b = tN / (tN + tB); L = lam(d, dh, mat)
    dmin = float(d[np.argmax(b > 0.9)] * 1e9); ok = (b > 0.9) & (L < 1)
    dmax_l = float(d[np.where(L < 1)[0][-1]] * 1e9) if np.any(L < 1) else None
    win[mat] = dict(d_min_brownian=round(dmin, 2), d_max_lambda1=round(dmax_l, 2) if dmax_l else None,
                    window=[round(float(d[ok][0] * 1e9), 2), round(float(d[ok][-1] * 1e9), 2)] if ok.any() else None)
    R.setdefault('curves', {})[mat] = dict(d=(d * 1e9).round(3).tolist(), b=b.round(5).tolist(), lam=L.round(5).tolist())
R['window'] = win
# sensitivity at full coverage (monodisperse) and amplitude per mass
sens = {}
for mat in MAT:
    rows = []
    for dc in np.arange(6, 30.01, 0.5) * 1e-9:
        f0, a0 = peak(dc, DELTA, mat); f1, _ = peak(dc, DELTA, mat, cover=1.0, layer=LAYER)
        rows.append([round(dc * 1e9, 2), f0, 100 * (1 - f1 / f0), a0])
    sens[mat] = rows
R['sens'] = sens
# design points: smallest core inside the window
design = {'Fe3O4': win['Fe3O4']['window'][0] if win['Fe3O4']['window'] else 20.0, 'CoFe2O4': win['CoFe2O4']['window'][0]}
R['design'] = {}
for mat, dc in design.items():
    f0, a0 = peak(dc * 1e-9, DELTA, mat); f1, _ = peak(dc * 1e-9, DELTA, mat, cover=1.0, layer=LAYER)
    R['design'][mat] = dict(d=dc, dh=round(dc + 2 * DELTA * 1e9, 1), f0=round(f0, 1), shift_full=round(100 * (1 - f1 / f0), 2), amp=a0,
                            lam=round(float(lam(dc * 1e-9, dc * 1e-9 + 2 * DELTA, mat)), 3))
ref = R['design']['Fe3O4']['amp']
for mat in R['design']: R['design'][mat]['amp_rel'] = round(R['design'][mat]['amp'] / ref, 3)
# polydispersity
poly = {}
for mat, dc in design.items():
    rows = []
    for s in (0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3):
        f0, a0 = peak(dc * 1e-9, DELTA, mat, sigma=s); f1, _ = peak(dc * 1e-9, DELTA, mat, sigma=s, cover=1.0, layer=LAYER)
        fr = np.logspace(-1, 9, 2001); sp = spectrum(fr, dc * 1e-9, DELTA, mat, sigma=s)
        hm = fr[sp >= sp.max() / 2]; width = float(np.log10(hm[-1] / hm[0]))
        rows.append(dict(sigma=s, shift=round(100 * (1 - f1 / f0), 2), f0=round(f0, 1), width_decades=round(width, 2)))
    poly[mat] = rows
R['poly'] = poly
# example spectra for figure
fr = np.logspace(1, 7, 600)
SPD = {'Fe3O4': design['Fe3O4'], 'CoFe2O4': 16.0}
R['spectra'] = {mat: dict(d0=SPD[mat], f=fr.tolist(), free=spectrum(fr, SPD[mat] * 1e-9, DELTA, mat, sigma=0.1).tolist(),
                          bound=spectrum(fr, SPD[mat] * 1e-9, DELTA, mat, sigma=0.1, cover=1.0, layer=LAYER).tolist()) for mat in MAT}
# practical CoFe2O4 design for a polydisperse sample (sigma 0.1): median 16 nm keeps > 96 % of the sample in the window
# for every effective anisotropy from 40 to 250 kJ/m3 (see noise_lod.py)
DP = 16.0
dc = DP * 1e-9
f25, _ = peak(dc, DELTA, 'CoFe2O4', sigma=0.1, T=298.15); f26, _ = peak(dc, DELTA, 'CoFe2O4', sigma=0.1, T=299.15)
R['temp_per_K'] = round(100 * (f26 / f25 - 1), 2)
f37w, _ = peak(dc, DELTA, 'CoFe2O4', sigma=0.1, T=310.15)
R['plasma'] = {}
for eta in (1.1e-3, 1.2e-3, 1.3e-3):
    fp, _ = peak(dc, DELTA, 'CoFe2O4', sigma=0.1, T=310.15, eta=eta); R['plasma'][str(eta)] = round(100 * (1 - fp / f37w), 1)
R['f37_water'] = round(f37w, 1); R['f25_water'] = round(f25, 1)
# binding layer thickness sweep
R['layer'] = {}
for L in (2e-9, 3e-9, 5e-9, 7.5e-9, 10e-9, 15e-9):
    f1, _ = peak(dc, DELTA, 'CoFe2O4', sigma=0.1, cover=1.0, layer=L); R['layer'][str(round(L * 1e9, 1))] = round(100 * (1 - f25 / f25 * f1 / f25), 2)
# dose-response and LOD (Langmuir), CoFe2O4 design, sigma 0.1
th = np.linspace(0, 1, 401); sh = []
for c in th:
    f1, _ = peak(dc, DELTA, 'CoFe2O4', sigma=0.1, cover=c, layer=LAYER); sh.append(100 * (1 - f1 / f25))
sh = np.array(sh); R['shift_vs_cover'] = dict(theta=th.round(4).tolist(), shift=sh.round(4).tolist())
lod = {}
for noise in (0.25, 0.5, 1.0, 2.0):
    tmin = float(np.interp(3 * noise, sh, th))
    lod[str(noise)] = dict(theta_min=round(tmin, 4), lod_over_Kd=round(tmin / (1 - tmin), 4),
                           lod_pM_Kd1nM=round(1000 * tmin / (1 - tmin), 1), temp_equiv_K=round(3 * noise / abs(R['temp_per_K']), 2))
R['lod'] = lod
# effect of the anisotropy uncertainty of CoFe2O4 on the Brownian window
R['K_sens'] = {}
for K in (40e3, 45e3, 50e3, 100e3, 180e3, 250e3):
    dd = np.linspace(5e-9, 20e-9, 1501); tN, tB, _ = taus(dd, dd + 2 * DELTA, 'CoFe2O4', T0, K=K)
    R['K_sens'][str(int(K / 1e3))] = round(float(dd[np.argmax(tN / (tN + tB) > 0.9)] * 1e9), 2)
f0p, a0p = peak(dc, DELTA, 'CoFe2O4', sigma=0.1); f1p, _ = peak(dc, DELTA, 'CoFe2O4', sigma=0.1, cover=1.0, layer=LAYER)
f0m, a0m = peak(design['Fe3O4'] * 1e-9, DELTA, 'Fe3O4', sigma=0.1)
R['practical'] = dict(d0=DP, sigma=0.1, f0=round(f0p, 1), shift_full=round(100 * (1 - f1p / f0p), 2), amp_rel_to_Fe3O4=round(a0p / a0m, 3))
json.dump(R, open('results.json', 'w'))
print(json.dumps({k: v for k, v in R.items() if k not in ('curves', 'sens', 'spectra', 'shift_vs_cover')}, indent=1))

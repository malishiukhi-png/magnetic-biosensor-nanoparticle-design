"""Values quoted in the text that are not stored by analysis.py / shell.py / window_frac.py -> extra3.json"""
import json, numpy as np
from model3 import *
T0 = 298.15; DELTA = 12e-9
R = json.load(open('results.json'))
def lam(d, dh, mat): m = MAT[mat]['Ms'] * np.pi * d ** 3 / 6; return mu0 * m ** 2 / (4 * np.pi * dh ** 3 * kB * T0)
out = {}
# 1) thinnest shell giving a usable Fe3O4 window (b > 0.9 and lambda < 1), 0.1 nm steps
d = np.linspace(5e-9, 35e-9, 3001)
for delta in np.arange(10.0, 13.01, 0.1) * 1e-9:
    dh = d + 2 * delta; tN, tB, _ = taus(d, dh, 'Fe3O4', T0)
    if ((tN / (tN + tB) > 0.9) & (lam(d, dh, 'Fe3O4') < 1)).any():
        out['Fe3O4_min_shell_nm'] = round(delta * 1e9, 1); break
# 2) chi'' amplitude per mass of the practical CoFe2O4 design (12 nm, sigma 0.1) relative to monodisperse 19.6 nm Fe3O4
_, a0p = peak(12e-9, DELTA, 'CoFe2O4', sigma=0.1)
out['amp_rel_practical_to_mono_Fe3O4'] = round(a0p / R['design']['Fe3O4']['amp'], 3)
# 3) mass fraction of the practical design (median 12 nm, sigma 0.1) inside the CoFe2O4 window
lo, hi = R['window']['CoFe2O4']['window']; s = 0.1; d0 = 12.0
x = np.linspace(np.log(d0) - 5 * s, np.log(d0) + 5 * s, 4001); dd = np.exp(x)
w = np.exp(-(x - np.log(d0)) ** 2 / (2 * s * s)) * dd ** 3; w /= w.sum()
out['practical_mass_in_window_pct'] = round(100 * float(w[(dd >= lo) & (dd <= hi)].sum()), 1)
json.dump(out, open('extra3.json', 'w'), indent=1); print(out)

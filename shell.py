"""Shell-thickness trade-off: Brownian window (b>0.9, lambda<1) and full-coverage shift vs shell thickness."""
import json, numpy as np
from model3 import *
T0 = 298.15
def lam(d, dh, mat): m = MAT[mat]['Ms'] * np.pi * d ** 3 / 6; return mu0 * m ** 2 / (4 * np.pi * dh ** 3 * kB * T0)
out = {}
d = np.linspace(5e-9, 35e-9, 601)  # 0.05 nm steps, the same grid as analysis.py
for mat in MAT:
    rows = []
    for delta in np.arange(5, 30.1, 1.0) * 1e-9:
        dh = d + 2 * delta; tN, tB, _ = taus(d, dh, mat, T0); b = tN / (tN + tB); L = lam(d, dh, mat)
        ok = (b > 0.9) & (L < 1)
        if ok.any():
            dc = float(d[ok][0]); f0, _ = peak(dc, delta, mat); f1, _ = peak(dc, delta, mat, cover=1.0, layer=5e-9)
            rows.append(dict(delta=round(delta * 1e9, 1), d=round(dc * 1e9, 2), dh=round((dc + 2 * delta) * 1e9, 1), shift=round(100 * (1 - f1 / f0), 2)))
        else:
            rows.append(dict(delta=round(delta * 1e9, 1), d=None, dh=None, shift=None))
    out[mat] = rows
json.dump(out, open('shell.json', 'w'), indent=1)
for mat in out:
    print(mat, [(r['delta'], r['d'], r['shift']) for r in out[mat] if r['delta'] in (5, 8, 10, 12, 15, 20, 25, 30)])

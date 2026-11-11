"""Mass fraction of a log-normal sample inside the usable window (Brownian fraction > 0.9 and lambda < 1), optimized over the median."""
import json, numpy as np
R = json.load(open('results.json'))
W = {m: R['window'][m]['window'] for m in R['window']}
out = {}
for mat, (lo, hi) in W.items():
    rows = []
    for s in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3):
        best = (0, None)
        for d0 in np.linspace(lo * 0.7, hi * 1.3, 600):
            x = np.linspace(np.log(d0) - 5 * s, np.log(d0) + 5 * s, 4001); d = np.exp(x)
            w = np.exp(-(x - np.log(d0)) ** 2 / (2 * s * s)) * d ** 3; w /= w.sum()
            fr = w[(d >= lo) & (d <= hi)].sum()
            if fr > best[0]: best = (fr, d0)
        rows.append(dict(sigma=s, frac=round(100 * best[0], 1), d0=round(best[1], 2)))
    out[mat] = rows
json.dump(out, open('window_frac.json', 'w'), indent=1); print(json.dumps(out, indent=1))

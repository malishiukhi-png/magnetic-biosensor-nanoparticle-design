"""Limit of detection with a read-out noise that depends on the signal amplitude and the particle concentration,
and sensitivity of the CoFe2O4 window to the cubic-anisotropy barrier (K1 V/4).  -> noise_lod.json
Same relaxation, binding and bridging model as binding.py; here the design (material, median size, width, K) is a parameter.

Noise model (Section 2.7 of the manuscript)
    eps_peak(cm)  = sqrt( eps0**2 + ( kappa * n_chi / chi2max(cm) )**2 )            relative error of the fitted peak frequency
    eps_flank(cm) = sqrt( eps0**2 + ( n_chi / ( sqrt(M) * chi2(3f0, cm) ) )**2 )    relative error of chi'' at 3 f0
n_chi    noise-equivalent volume susceptibility per spectral point (SI)
kappa    = sqrt(3/(2*rho_u)); rho_u spectral points per e-fold of frequency (Fisher information for the position of a Debye peak)
M        number of spectral points of the sweep; the single-frequency read-out uses the same total acquisition time
eps0     multiplicative floor (drift, repeatability), independent of the concentration
"""
import json, numpy as np
from model3 import *

T0 = 298.15; DELTA = 12e-9; LAYER = 5e-9; KD = 1e-9; NS = 20; DF = 1.8; NA = 6.02214076e23
ETA = eta_water(T0); KS = 8 * kB * T0 / (3 * ETA)
RHO_U = 10.0; M = 46; KAPPA = np.sqrt(3 / (2 * RHO_U)); EPS0 = 0.4e-2
FG = np.logspace(1, 7, 2001)
LGRID = np.logspace(-13, -5, 97)
K13 = 13e3


class Design:
    def __init__(self, mat, d0, sig, K=None):
        self.mat, self.d0, self.sig = mat, d0, sig
        self.K = MAT[mat]['K'] if K is None else K
        m = MAT[mat]
        if sig > 0:
            x = np.linspace(np.log(d0) - 4 * sig, np.log(d0) + 4 * sig, 201); self.dd = np.exp(x)
            g = np.exp(-(x - np.log(d0))**2 / (2 * sig**2))
        else:
            self.dd = np.array([d0]); g = np.array([1.0])
        w = g * self.dd**3; self.wmass = w / w.sum()
        self.Vn = np.sum(g * np.pi * self.dd**3 / 6) / g.sum()
        self.chi0 = mu0 * m['Ms']**2 * (np.pi * self.dd**3 / 6) / (3 * kB * T0) / m['rho']   # m3/kg, per kg of magnetic material
        self.rho = m['rho']
        self.free = np.zeros(12); self.free[0] = 1.0
        s = self.spec(FG, 0.0, self.free); i = int(np.argmax(s))
        self.F0 = self._refine(s); self.A = float(s.max())                                   # chi''max per unit mass (m3/kg)
        self.FR = 3 * self.F0; self.X0 = self.spec(np.array([self.FR]), 0.0, self.free)[0]

    def _refine(self, s):
        i = int(np.argmax(s)); i = min(max(i, 1), len(FG) - 2)
        y0, y1, y2 = np.log(s[i - 1:i + 2]); den = y0 - 2 * y1 + y2
        off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        return float(np.exp(np.log(FG[i]) + off * (np.log(FG[i + 1]) - np.log(FG[i]))))

    def spec(self, f, theta, pk):
        f = np.atleast_1d(f); dh = self.dd + 2 * DELTA + 2 * LAYER * theta
        tN, tB1, _ = taus(self.dd, dh, self.mat, T0, K=self.K)
        out = np.zeros(len(f))
        for k, p in enumerate(pk, start=1):
            if p < 1e-12: continue
            tB = tB1 * k**(3 / DF); t = tN * tB / (tN + tB); wt = np.outer(2 * np.pi * f, t)
            out += p * (self.wmass * self.chi0 * wt / (1 + wt**2)).sum(axis=1)
        return out

    def Pt(self, cm):               # cm in ug/mL -> particle molarity
        return cm * 1e-3 / (self.rho * self.Vn) / NA / 1e3

    def signals(self, theta, pk):
        fp = self._refine(self.spec(FG, theta, pk)); xs = self.spec(np.array([self.FR]), theta, pk)[0]
        return 100 * (1 - fp / self.F0), 100 * (1 - xs / self.X0)

    def eps(self, cm, n):
        chi2 = self.A * cm * 1e-3                                            # cm [ug/mL] = 1e-3 kg/m3
        e_p = np.sqrt(EPS0**2 + (KAPPA * n / chi2)**2)
        e_f = np.sqrt(EPS0**2 + (n / (np.sqrt(M) * self.X0 * cm * 1e-3))**2)
        return 100 * e_p, 100 * e_f


def occupancy(Lt, Pt, Ns=NS, Kd=KD):
    St = Ns * Pt; b = Lt + St + Kd
    return (b - np.sqrt(np.maximum(b * b - 4 * Lt * St, 0))) / 2 / St


def cluster(theta, Pt, t_inc, alpha, Ns=NS, kmax=12):
    n = Pt * 1e3 * NA; eps = alpha * 2 * theta * (1 - theta); tau = eps * KS * n * t_inc / 2; b = Ns * theta
    if b < 1: tau = min(tau, b / (1 - b))
    k = np.arange(1, kmax + 1); pk = k * tau**(k - 1) / (1 + tau)**(k + 1); pk[-1] += max(0.0, 1 - pk.sum())
    return pk


def first_cross(sig, thr, L):
    # lowest concentration at which the signal reaches the threshold (rising branch only)
    if sig.max() < thr: return None
    i = int(np.argmax(sig >= thr))
    if i == 0: return float(L[0])
    return float(np.exp(np.interp(thr, [sig[i - 1], sig[i]], np.log([L[i - 1], L[i]]))))


def curves(D, cm, alpha=None):
    Pt = D.Pt(cm); th = occupancy(LGRID, Pt); sf, sx = [], []
    for t in th:
        pk = cluster(t, Pt, 600.0, alpha) if alpha else D.free
        a, b = D.signals(t, pk); sf.append(a); sx.append(b)
    return np.array(sf), np.array(sx)


CMS = np.array([0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000.0])
NCHI = (1e-7, 1e-6, 1e-5)


def lod_table(D, alpha):
    out = {}
    for cm in CMS:
        sf, sx = curves(D, cm, alpha)
        for n in NCHI:
            ep, ef = D.eps(cm, n)
            out.setdefault(str(n), []).append(dict(cm=float(cm), eps_peak=round(float(ep), 3), eps_flank=round(float(ef), 3),
                                                   lod_peak=first_cross(sf, 3 * ep, LGRID), lod_flank=first_cross(sx, 3 * ef, LGRID)))
    return out


def best(rows, key):
    v = [(r[key], r['cm']) for r in rows if r[key] is not None]
    return min(v) if v else (None, None)


if __name__ == '__main__':
    R = {}
    # ---- 1) CoFe2O4 Brownian/chaining window versus effective anisotropy ---------------------------------------
    def lam(d, dh, Ms): m = Ms * np.pi * d**3 / 6; return mu0 * m**2 / (4 * np.pi * dh**3 * kB * T0)
    def window(K, lam_max=1.0, delta=DELTA):
        d = np.linspace(4e-9, 40e-9, 7201); dh = d + 2 * delta; tN, tB, _ = taus(d, dh, 'CoFe2O4', T0, K=K)
        ok = (tN / (tN + tB) > 0.9) & (lam(d, dh, MAT['CoFe2O4']['Ms']) < lam_max)
        return [round(float(d[ok][0] * 1e9), 1), round(float(d[ok][-1] * 1e9), 1)] if ok.any() else None
    def mass_in(d0, s, lo, hi):
        x = np.linspace(np.log(d0) - 5 * s, np.log(d0) + 5 * s, 4001); dd = np.exp(x); w = np.exp(-(x - np.log(d0))**2 / (2 * s * s)) * dd**3; w /= w.sum()
        return float(w[(dd >= lo) & (dd <= hi)].sum() * 100)
    R['K_window'] = {}
    for K in (40e3, 45e3, 50e3, 100e3, 180e3, 250e3):
        w = window(K); row = dict(window=w, width=round(w[1] - w[0], 1))
        row['mass_in'] = {str(d0): round(mass_in(d0, 0.1, *w), 1) for d0 in (12, 14, 15, 16, 17, 18)}
        grid = np.linspace(w[0] * 0.8, w[1] * 1.1, 400); m = [mass_in(x, 0.1, *w) for x in grid]
        row['best_median'] = round(float(grid[int(np.argmax(m))]), 1); row['best_mass'] = round(float(max(m)), 1)
        for lm in (0.5, 2.0): row['window_lam%g' % lm] = window(K, lm)
        R['K_window'][str(int(K / 1e3))] = row
    # peak shift and response fraction of candidate medians for K_eff = 45 and 180 kJ/m3
    R['candidates'] = {}
    for K in (45e3, 180e3):
        for d0 in (12e-9, 14e-9, 16e-9, 18e-9):
            f0, a0 = peak(d0, DELTA, 'CoFe2O4', sigma=0.1, K=K); f1, _ = peak(d0, DELTA, 'CoFe2O4', sigma=0.1, K=K, cover=1.0, layer=LAYER)
            lmd = float(lam(d0, d0 + 2 * DELTA, MAT['CoFe2O4']['Ms']))
            R['candidates']['%d_%d' % (K / 1e3, d0 * 1e9)] = dict(f0=round(f0, 1), shift=round(100 * (1 - f1 / f0), 1), lam=round(lmd, 3), A=a0)
    # ---- 2) LOD with signal- and concentration-dependent noise --------------------------------------------------
    designs = {'CFO12_K180': ('CoFe2O4', 12e-9, 0.1, 180e3), 'CFO16_K180': ('CoFe2O4', 16e-9, 0.1, 180e3),
               'CFO12_K45': ('CoFe2O4', 12e-9, 0.1, 45e3), 'CFO16_K45': ('CoFe2O4', 16e-9, 0.1, 45e3), 'CFO18_K45': ('CoFe2O4', 18e-9, 0.1, 45e3),
               'FO18_K20': ('Fe3O4', 18.5e-9, 0.1, 20e3), 'FO20_K13': ('Fe3O4', 19.8e-9, 0.1, 13e3)}
    R['lod'] = {}
    for name, (mat, d0, s, K) in designs.items():
        D = Design(mat, d0, s, K)
        e = dict(F0=round(D.F0, 1), A=D.A, Pt_nM_per_ugml=round(D.Pt(1.0) * 1e9, 4), sites_nM_per_ugml=round(NS * D.Pt(1.0) * 1e9, 3))
        e['monolayer'] = lod_table(D, None); e['bridging_1e-3'] = lod_table(D, 1e-3)
        for mode in ('monolayer', 'bridging_1e-3'):
            e['best_' + mode] = {}
            for n in NCHI:
                rows = e[mode][str(n)]
                e['best_' + mode][str(n)] = dict(peak=best(rows, 'lod_peak'), flank=best(rows, 'lod_flank'))
        R['lod'][name] = e
        print(name, json.dumps({m: e['best_' + m] for m in ('monolayer', 'bridging_1e-3')}), flush=True)
    R['params'] = dict(kappa=KAPPA, M=M, eps0=EPS0, rho_u=RHO_U)
    json.dump(R, open('noise_lod.json', 'w'), indent=1)
    print(json.dumps(R['K_window'], indent=0)); print(json.dumps(R['candidates']))

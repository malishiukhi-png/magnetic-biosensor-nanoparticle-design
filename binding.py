"""Target depletion and bridging aggregation for the practical CoFe2O4 design -> binding.json, fig_binding.png
Builds on model3.py (same particles, shell, layer and temperature as analysis.py)."""
import json, numpy as np
from model3 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

T0 = 298.15; DELTA = 12e-9; LAYER = 5e-9; MAT_ = 'CoFe2O4'; D0 = 12e-9; SIG = 0.1
KD = 1e-9                     # dissociation constant (M)
NS = 20                       # functional binding sites per particle (assumed; 5-50 tested)
DF = 1.8                      # fractal dimension used for the hydrodynamic size of k-particle clusters
NA = 6.02214076e23
ETA = eta_water(T0)
KS = 8 * kB * T0 / (3 * ETA)  # Smoluchowski rate constant for Brownian coagulation of equal spheres (m3/s)

# --- size distribution (same grid and weights as model3.spectrum) ---------------------------------
x = np.linspace(np.log(D0) - 4 * SIG, np.log(D0) + 4 * SIG, 301); dd = np.exp(x)
g = np.exp(-(x - np.log(D0))**2 / (2 * SIG**2)); wmass = g * dd**3; wmass /= wmass.sum()
m = MAT[MAT_]
Vmean_number = np.sum(g * np.pi * dd**3 / 6) / g.sum()           # number-averaged core volume
chi0 = mu0 * m['Ms']**2 * (np.pi * dd**3 / 6) / (3 * kB * T0) / m['rho']


def particle_molar(cm_ug_per_ml):
    """Molar concentration of particles for a magnetic-material mass concentration in ug/mL."""
    n = cm_ug_per_ml * 1e-3 / (m['rho'] * Vmean_number)          # particles per m3
    return n / NA / 1e3                                            # mol/L


def occupancy(Lt, Pt, Ns=NS, Kd=KD):
    """Fraction of occupied sites with ligand depletion (exact 1:1 site binding)."""
    St = Ns * Pt
    b = Lt + St + Kd
    C = (b - np.sqrt(np.maximum(b * b - 4 * Lt * St, 0))) / 2
    return C / St


def spectrum_k(f, theta, pk):
    """chi''(f) for coverage theta and mass fractions pk[k-1] of particles in k-particle clusters."""
    f = np.atleast_1d(f)
    dh = dd + 2 * DELTA + 2 * LAYER * theta
    tN, tB1, _ = taus(dd, dh, MAT_, T0)
    out = np.zeros(len(f))
    for k, p in enumerate(pk, start=1):
        if p < 1e-12:
            continue
        tB = tB1 * k**(3 / DF)
        t = tN * tB / (tN + tB)
        wt = np.outer(2 * np.pi * f, t)
        out += p * (wmass * chi0 * wt / (1 + wt**2)).sum(axis=1)
    return out


FG = np.logspace(1, 7, 3001)


def fpeak(theta, pk):
    s = spectrum_k(FG, theta, pk); i = int(np.argmax(s))
    y0, y1, y2 = np.log(s[i - 1:i + 2]); den = y0 - 2 * y1 + y2
    off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
    return float(np.exp(np.log(FG[i]) + off * (np.log(FG[i + 1]) - np.log(FG[i]))))


def cluster_fractions(theta, Pt, t_inc, alpha, Ns=NS, kmax=12):
    """Mass fractions in k-mers after Smoluchowski aggregation with bridging efficiency
    eps = alpha * 2 theta (1 - theta); bonds per particle capped by bound targets per particle."""
    n = Pt * 1e3 * NA
    eps = alpha * 2 * theta * (1 - theta)
    tau = eps * KS * n * t_inc / 2
    b = Ns * theta                                   # bound targets per particle
    if b < 1:
        tau = min(tau, b / (1 - b))
    k = np.arange(1, kmax + 1)
    pk = k * tau**(k - 1) / (1 + tau)**(k + 1)
    pk[-1] += max(0.0, 1 - pk.sum())
    return pk, tau / (1 + tau)


free = np.zeros(12); free[0] = 1.0
F0 = fpeak(0.0, free); FR = 3 * F0          # single-frequency read-out on the high-frequency flank
X0 = spectrum_k(np.array([FR]), 0.0, free)[0]


def signals(theta, pk):
    fp = fpeak(theta, pk)
    xs = spectrum_k(np.array([FR]), theta, pk)[0]
    return 100 * (1 - fp / F0), 100 * (1 - xs / X0)


R = dict(F0=round(F0, 1), FR=round(FR, 1), KS=KS, Ns=NS, Kd_nM=KD * 1e9, Df=DF)
Lgrid = np.logspace(-13, -5, 161)       # target concentration (M)
CMS = (1, 3, 10, 30, 100)
NOISE = 0.5                              # % relative noise of the read-out


def lod(sig, L):
    thr = 3 * NOISE
    i = np.argmax(sig >= thr)
    if sig.max() < thr:
        return None
    if i == 0:
        return float(L[0])
    return float(np.exp(np.interp(thr, [sig[i - 1], sig[i]], np.log([L[i - 1], L[i]]))))


# 1) depletion only (monolayer growth, no bridging)
dep = {}
curves = {}
for cm in CMS:
    Pt = particle_molar(cm)
    th = occupancy(Lgrid, Pt)
    sf, sx = np.array([signals(t, free) for t in th]).T
    dep[str(cm)] = dict(particle_nM=round(Pt * 1e9, 3), sites_nM=round(NS * Pt * 1e9, 2),
                        lod_fp_pM=round(lod(sf, Lgrid) * 1e12, 1), lod_chi_pM=round(lod(sx, Lgrid) * 1e12, 1))
    curves['dep_%s' % cm] = dict(fp=sf.tolist(), chi=sx.tolist())
# no-depletion reference (target in excess, as in the original analysis)
th_x = Lgrid / (Lgrid + KD)
sf_x, sx_x = np.array([signals(t, free) for t in th_x]).T
R['excess'] = dict(lod_fp_pM=round(lod(sf_x, Lgrid) * 1e12, 1), lod_chi_pM=round(lod(sx_x, Lgrid) * 1e12, 1))
curves['excess'] = dict(fp=sf_x.tolist(), chi=sx_x.tolist())
R['depletion'] = dep

# 2) sites per particle (10 ug/mL)
R['Ns_sens'] = {}
for ns in (5, 20, 50):
    Pt = particle_molar(10)
    th = occupancy(Lgrid, Pt, Ns=ns)
    sf = np.array([signals(t, free)[0] for t in th])
    R['Ns_sens'][str(ns)] = round(lod(sf, Lgrid) * 1e12, 1)

# 3) bridging aggregation (10 ug/mL, 10 min incubation)
brid = {}
CM_B = 10; Pt = particle_molar(CM_B); th = occupancy(Lgrid, Pt)
for alpha in (1e-4, 1e-3, 1e-2):
    sf, sx, bonds = [], [], []
    for t in th:
        pk, bp = cluster_fractions(t, Pt, 600.0, alpha)
        a, b = signals(t, pk); sf.append(a); sx.append(b); bonds.append(bp)
    sf, sx = np.array(sf), np.array(sx)
    ipk = int(np.argmax(sx))
    brid[str(alpha)] = dict(lod_fp_pM=round(lod(sf, Lgrid) * 1e12, 1) if lod(sf, Lgrid) else None,
                            lod_chi_pM=round(lod(sx, Lgrid) * 1e12, 1) if lod(sx, Lgrid) else None,
                            max_chi_signal=round(float(sx.max()), 1), L_at_max_nM=round(float(Lgrid[ipk] * 1e9), 2),
                            chi_at_1uM=round(float(np.interp(-6, np.log10(Lgrid), sx)), 1),
                            bonds_at_100pM=round(float(np.interp(-10, np.log10(Lgrid), bonds)), 4))
    curves['brid_%g' % alpha] = dict(fp=sf.tolist(), chi=sx.tolist(), bonds=bonds)
R['bridging_10ugml_10min'] = brid
# incubation time needed for bridging to reach the cap at the LOD concentration (alpha = 1e-3)
R['tinc'] = {}
for tinc in (60, 300, 600, 1800):
    sx = []
    for t in th:
        pk, _ = cluster_fractions(t, Pt, tinc, 1e-3)
        sx.append(signals(t, pk)[1])
    R['tinc'][str(tinc)] = round(lod(np.array(sx), Lgrid) * 1e12, 1)
# bridging LOD versus particle concentration (alpha = 1e-3, 10 min)
R['brid_vs_cm'] = {}
for cm in (1, 3, 10, 30, 100):
    Ptc = particle_molar(cm); thc = occupancy(Lgrid, Ptc); sf, sx = [], []
    for t in thc:
        pk, _ = cluster_fractions(t, Ptc, 600.0, 1e-3)
        a, b = signals(t, pk); sf.append(a); sx.append(b)
    R['brid_vs_cm'][str(cm)] = dict(lod_fp_pM=round(lod(np.array(sf), Lgrid) * 1e12, 1), lod_chi_pM=round(lod(np.array(sx), Lgrid) * 1e12, 1))
# drift of the flank read-out per K (25 C) for comparison with the temperature effect on f_p
x25 = spectrum_k(np.array([FR]), 0.0, free)[0]
TT = 299.15
dh = dd + 2 * DELTA
tN, tB1, _ = taus(dd, dh, MAT_, TT)
t = tN * tB1 / (tN + tB1); wt = 2 * np.pi * FR * t
chi0T = chi0 * T0 / TT
x26 = (wmass * chi0T * wt / (1 + wt**2)).sum()
R['flank_drift_per_K'] = round(100 * (x26 / x25 - 1), 2)
json.dump(dict(R=R, L=Lgrid.tolist(), curves=curves), open('binding.json', 'w'))
print(json.dumps(R, indent=1))


# ---- Figure 5 -------------------------------------------------------------------------------------
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 10, 'axes.labelweight': 'bold',
                     'axes.linewidth': 0.8, 'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.top': True, 'ytick.right': True, 'xtick.minor.visible': True, 'ytick.minor.visible': True,
                     'legend.frameon': False, 'legend.fontsize': 8, 'savefig.dpi': 300, 'lines.linewidth': 1.8})
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']; INK = '#0b0b0b'; MUTED = '#898781'; CMM = 1 / 2.54
fig, axs = plt.subplots(1, 2, figsize=(16 * CMM, 6.4 * CMM))
Lx = Lgrid * 1e9
ax = axs[0]
ax.semilogx(Lx, curves['excess']['fp'], color=INK, ls=':', label='Target excess')
for k, cm in enumerate((1, 10, 100)):
    ax.semilogx(Lx, curves['dep_%d' % cm]['fp'], color=C[k], ls=['-', '--', '-.'][k], label='%d µg mL$^{-1}$' % cm)
ax.axhline(1.5, color=MUTED, lw=1); ax.text(1.3e-4, 2.4, '3 × 0.5% noise', fontsize=7, color='#52514e', ha='left')
ax.set_xlim(1e-4, 1e4); ax.set_ylim(0, 55)
ax.set_xlabel('Target concentration (nM)'); ax.set_ylabel('Peak shift (%)')
ax.legend(loc='upper left', bbox_to_anchor=(0.0, 0.88)); ax.text(0.03, 0.95, '(a)', transform=ax.transAxes, fontweight='bold', va='top')
ax = axs[1]
ax.semilogx(Lx, curves['dep_10']['chi'], color=INK, ls=':', label='No bridging')
for k, a in enumerate((1e-4, 1e-3, 1e-2)):
    ax.semilogx(Lx, curves['brid_%g' % a]['chi'], color=C[k], ls=['-', '--', '-.'][k], label='α = 10$^{%d}$' % round(np.log10(a)))
ax.axhline(1.5, color=MUTED, lw=1)
ax.set_xlim(1e-4, 1e4); ax.set_ylim(0, 105)
ax.set_xlabel('Target concentration (nM)'); ax.set_ylabel('Decrease of χ″ at 3$f_0$ (%)')
ax.legend(loc='upper left', bbox_to_anchor=(0.0, 0.88)); ax.text(0.03, 0.95, '(b)', transform=ax.transAxes, fontweight='bold', va='top')
fig.tight_layout(); fig.savefig('fig_binding.png'); plt.close(fig)

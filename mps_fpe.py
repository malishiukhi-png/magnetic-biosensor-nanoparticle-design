"""Nonlinear harmonic read-out: Brownian rotation of blocked dipoles in a sinusoidal field (Fokker-Planck equation in Legendre form).
Tested against the measured harmonic ratio of Zhong et al. (J. Magn. Magn. Mater. 517 (2021) 167408) and applied to the designs of this work.
Writes mps_fpe.json."""
import json
import numpy as np
from scipy.linalg import solve_banded

kB = 1.380649e-23; MU0 = 4e-7 * np.pi


def water_eta(T):
    """Dynamic viscosity of water, Pa s (Vogel form, as a function of T in K)."""
    return 2.414e-5 * 10**(247.8 / (T - 140.0))


def harmonics(tau, xi, f, nl=None, ncyc=None, nstep=600):
    """Complex amplitudes M1, M3 (in units of the saturation moment) in the periodic steady state.
    d a_n/dt = n(n+1)/(2 tau) * ( -a_n + xi(t) * (a_{n-1}/(2n-1) - a_{n+1}/(2n+3)) ),  a_0 = 1/2,  M/m = 2 a_1 / 3."""
    if nl is None:
        nl = int(max(24, 3 * xi + 16))
    if ncyc is None:                       # enough cycles to reach the periodic state; the field accelerates the relaxation
        ncyc = int(np.ceil(4 * tau * f / (1 + xi / 3))) + 4
    n = np.arange(1, nl + 1, dtype=float); r = n * (n + 1) / (2.0 * tau)
    w = 2 * np.pi * f; dt = 1.0 / f / nstep; a = np.zeros(nl)
    L, D, U = r / (2 * n - 1), -r, -r / (2 * n + 3)     # sub-diagonal, diagonal (relaxation), super-diagonal; the first two scale with xi
    tr = []
    for k in range(ncyc * nstep):
        x0, x1 = xi * np.cos(w * k * dt), xi * np.cos(w * (k + 1) * dt)
        # Crank-Nicolson
        rhs = a + 0.5 * dt * (D * a)
        rhs[1:] += 0.5 * dt * x0 * L[1:] * a[:-1]; rhs[:-1] += 0.5 * dt * x0 * U[:-1] * a[1:]
        rhs[0] += 0.5 * dt * (x0 + x1) * L[0] * 0.5                 # a_0 = 1/2
        ab = np.zeros((3, nl))
        ab[0, 1:] = -0.5 * dt * x1 * U[:-1]; ab[1] = 1 - 0.5 * dt * D; ab[2, :-1] = -0.5 * dt * x1 * L[1:]
        a = solve_banded((1, 1), ab, rhs)
        if k >= (ncyc - 1) * nstep:
            tr.append(2.0 * a[0] / 3.0)
    tr = np.array(tr); F = np.fft.rfft(tr) / nstep * 2
    return F[1], F[3]


def tau_B(dh, eta, T):
    return 3 * eta * (np.pi * dh**3 / 6) / (kB * T)


def lognormal_nodes(mu, s, n=25):
    x = np.linspace(mu - 3.5 * s, mu + 3.5 * s, n); w = np.exp(-(x - mu)**2 / (2 * s * s)); return np.exp(x), w / w.sum()


def ratio(m, B, f, eta, T, mu, s, scale=1.0):
    """M3/M1 of an ensemble with equal moments and a log-normal hydrodynamic diameter (nm); scale multiplies every diameter."""
    d, w = lognormal_nodes(mu, s); xi = m * B / (kB * T); M1 = M3 = 0j
    for di, wi in zip(d, w):
        h1, h3 = harmonics(tau_B(di * scale * 1e-9, eta, T), xi, f)
        M1 += wi * h1; M3 += wi * h3
    return abs(M3) / abs(M1), abs(M1), abs(M3)


R = {}
# ---- Validation: rigid-dipole limits ----------------------------------------------------------------------------
# weak field: M1 must equal the Debye response (xi/3)/(1 + i w tau) in magnitude
tau0 = 1e-3; f0 = 1 / (2 * np.pi * tau0); h1, _ = harmonics(tau0, 0.05, f0)
R['check_debye'] = dict(M1=abs(h1), analytic=0.05 / 3 / np.sqrt(2))
# quasi-static: ratio of the third to first harmonic of the Langevin function at xi = 5
xi = 5.0; t = np.linspace(0, 1, 4096, endpoint=False); lang = lambda x: 1 / np.tanh(x) - 1 / x
y = lang(xi * np.cos(2 * np.pi * t + 1e-9)); Y = np.fft.rfft(y) / 4096 * 2
h1, h3 = harmonics(1e-7, xi, 1.0, nstep=800, ncyc=3)
R['check_langevin'] = dict(ratio=abs(h3) / abs(h1), analytic=abs(Y[3]) / abs(Y[1]))
# ---- Zhong et al. 2021: BNF-80, streptavidin coated, 10 mT, moment and distribution as given in their Section 2.2 ---
T_z = 295.15; eta_z = water_eta(T_z); m_z = 8.5e-18; B_z = 10e-3; s_z = 0.3
freqs = [180, 300, 460, 620, 720, 840, 1020, 1418]
# measured ratio without IgG at 6.67 nM particles, read from their Fig. 5a (uncertainty about 0.01)
meas = {180: 0.347, 300: 0.283, 460: 0.226, 620: 0.186, 720: 0.165, 840: 0.141, 1020: 0.116, 1418: 0.076}
# the mean of ln(d) is fixed once by requiring the Debye peak of the distribution to give the measured tau = 0.82 ms (no fit to the harmonics)
tau_z = 0.82e-3; mu_z = np.log((tau_z * kB * T_z / (3 * eta_z) * 6 / np.pi)**(1 / 3) * 1e9)
R['zhong'] = dict(T=T_z, eta_mPas=eta_z * 1e3, m=m_z, B=B_z, sigma=s_z, median_dh_nm=float(np.exp(mu_z)),
                  f=freqs, measured=[meas[f] for f in freqs], model=[ratio(m_z, B_z, f, eta_z, T_z, mu_z, s_z)[0] for f in freqs])
# binding: the change of tau_B measured by AC susceptibility predicts the change of the ratio with no adjustable parameter (1.00 nM particles, 8.33 nM IgG, 180 Hz)
scale = (1.31 / 0.86)**(1 / 3); mu_b = np.log((0.86e-3 * kB * T_z / (3 * eta_z) * 6 / np.pi)**(1 / 3) * 1e9)
r0 = ratio(m_z, B_z, 180.0, eta_z, T_z, mu_b, s_z)[0]; r1 = ratio(m_z, B_z, 180.0, eta_z, T_z, mu_b, s_z, scale)[0]
R['zhong_binding'] = dict(scale=float(scale), model_R0=r0, model_R1=r1, model_dR=r0 - r1, measured_R0=0.337, measured_R1=0.247, measured_dR=0.337 - 0.247)
# numerical convergence at the most demanding point
R['check_convergence'] = dict(f=1418.0, R_default=R['zhong']['model'][-1], R_finer=None)
d_, w_ = lognormal_nodes(mu_z, s_z, 49); xi_ = m_z * B_z / (kB * T_z); M1 = M3 = 0j
for di, wi in zip(d_, w_):
    h1, h3 = harmonics(tau_B(di * 1e-9, eta_z, T_z), xi_, 1418.0, nstep=1500)
    M1 += wi * h1; M3 += wi * h3
R['check_convergence']['R_finer'] = abs(M3) / abs(M1)
# ---- Designs of this work: CoFe2O4 cores (blocked at these fields), 12 nm capture shell, 5 nm bound layer, log-normal core diameter, sigma = 0.1 ----
T = 298.15; eta = water_eta(T); Ms = 425e3; rho = 5290.0; shell = 12.0; layer = 5.0; s = 0.1
R['design'] = {}


def design(d0, B, f, bound):
    d, w = lognormal_nodes(np.log(d0), s, 15)
    w = w / w.sum(); M1 = M3 = 0j; M = 0.0
    for di, wi in zip(d, w):
        V = np.pi * (di * 1e-9)**3 / 6; m = Ms * V; xi = m * B / (kB * T)
        dh = (di + 2 * shell + (2 * layer if bound else 0)) * 1e-9
        h1, h3 = harmonics(tau_B(dh, eta, T), xi, f)
        M1 += wi * m * h1; M3 += wi * m * h3; M += wi * rho * V
    return abs(M3) / abs(M1), abs(M1) / M, abs(M3) / M


for B in (5e-3, 10e-3):
    for d0 in (12.0, 16.0, 18.0):
        f0 = 1 / (2 * np.pi * tau_B((d0 + 2 * shell) * 1e-9, eta, T)); best = None; curve = []
        for x in np.geomspace(0.1, 2.0, 12):
            f = f0 * x
            r0, a1, a3 = design(d0, B, f, False); r1, _, a3b = design(d0, B, f, True)
            row = dict(f=f, x=float(x), R0=r0, R1=r1, dR=r0 - r1, rel=(r0 - r1) / r0, A1=a1, A3=a3, A3_bound=a3b, fom=(r0 - r1) * a3 / r0)
            curve.append(row)
            if best is None or row['fom'] > best['fom']:
                best = row
        R['design']['%d_%d' % (d0, B * 1e3)] = dict(f0=f0, best=best, curve=curve)
    ref = R['design']['16_%d' % (B * 1e3)]['best']
    for d0 in (12, 16, 18):
        R['design']['%d_%d' % (d0, B * 1e3)]['fom_rel'] = R['design']['%d_%d' % (d0, B * 1e3)]['best']['fom'] / ref['fom']


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float('%.5g' % o)
    return o


json.dump(clean(R), open('mps_fpe.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in clean(R).items() if k != 'design'}, indent=0))
for k, v in R['design'].items():
    b = v['best']; print(k, 'f0 %.0f f %.0f R0 %.3f dR/R %.3f A3 %.3g fom_rel %.2f' % (v['f0'], b['f'], b['R0'], b['rel'], b['A3'], v['fom_rel']))

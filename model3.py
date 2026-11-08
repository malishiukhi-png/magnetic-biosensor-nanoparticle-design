"""AC-susceptibility (Debye) model of Brownian-relaxation biosensing with magnetic nanoparticles."""
import numpy as np
kB = 1.380649e-23; mu0 = 4e-7 * np.pi
MAT = {'Fe3O4':   dict(Ms=446e3, K=13e3,  rho=5180.0),     # magnetite (as in paper-1 model)
       'CoFe2O4': dict(Ms=425e3, K=180e3, rho=5290.0)}     # cobalt ferrite (bulk-like, assumed)
TAU0 = 1e-9

def eta_water(T):
    """Viscosity of water (Pa s), Vogel-type fit: 2.414e-5 * 10**(247.8/(T-140))."""
    return 2.414e-5 * 10 ** (247.8 / (T - 140.0))

def taus(d, dh, mat, T=298.15, eta=None, K=None):
    m = MAT[mat]; K = m['K'] if K is None else K; eta = eta_water(T) if eta is None else eta
    V = np.pi * d ** 3 / 6; VH = np.pi * dh ** 3 / 6
    tN = TAU0 * np.exp(np.minimum(K * V / (kB * T), 700)); tB = 3 * eta * VH / (kB * T)
    return tN, tB, tN * tB / (tN + tB)

def spectrum(f, d0, delta, mat, sigma=0.0, cover=0.0, layer=5e-9, T=298.15, eta=None, K=None, n=301):
    """chi'' (arbitrary units, per unit mass of magnetic material) for a log-normal size distribution.
    Hydrodynamic diameter dh = d + 2*delta + 2*layer*cover (cover: fractional surface coverage by the target)."""
    m = MAT[mat]
    if sigma > 0:
        x = np.linspace(np.log(d0) - 4 * sigma, np.log(d0) + 4 * sigma, n); d = np.exp(x)
        g = np.exp(-(x - np.log(d0)) ** 2 / (2 * sigma ** 2)); w = g * d ** 3; w /= w.sum()   # volume (mass) weights
    else:
        d = np.array([d0]); w = np.array([1.0])
    dh = d + 2 * delta + 2 * layer * cover
    _, _, t = taus(d, dh, mat, T, eta, K)
    chi0 = mu0 * m['Ms'] ** 2 * (np.pi * d ** 3 / 6) / (3 * kB * T) / m['rho']      # initial susceptibility per mass (x volume fraction)
    wt = np.outer(2 * np.pi * np.atleast_1d(f), t)
    return (w * chi0 * wt / (1 + wt ** 2)).sum(axis=1)

def peak(d0, delta, mat, fgrid=None, **kw):
    """Frequency (Hz) and height of the chi'' maximum (log-parabolic refinement on a log grid)."""
    f = np.logspace(-1, 9, 4001) if fgrid is None else fgrid
    s = spectrum(f, d0, delta, mat, **kw); i = int(np.argmax(s))
    if 0 < i < len(f) - 1:
        y0, y1, y2 = np.log(s[i - 1:i + 2]); den = y0 - 2 * y1 + y2
        off = 0.5 * (y0 - y2) / den if den != 0 else 0.0
        lf = np.log(f[i]) + off * (np.log(f[i + 1]) - np.log(f[i]))
        return float(np.exp(lf)), float(s[i])
    return float(f[i]), float(s[i])

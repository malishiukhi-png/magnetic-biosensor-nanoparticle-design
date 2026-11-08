"""Comparison of the Brownian-peak relation and of the binding-induced shift with published AC-susceptibility data -> compare_exp.json
Values taken from the literature are marked in the source list of the manuscript (Table 4)."""
import json, numpy as np
kB = 1.380649e-23
T = 298.15; ETA = 0.890e-3


def f_of_d(d, eta=ETA, T=T):
    return kB * T / (3 * eta * np.pi * d**3 / 6) / (2 * np.pi)


def d_of_tau(tau, eta=ETA, T=T):
    V = tau * kB * T / (3 * eta)
    return (6 * V / np.pi)**(1 / 3)


def d_of_f(f, **k):
    return d_of_tau(1 / (2 * np.pi * f), **k)


R = {}
# Enpuku et al. (2009): free markers, DLS size 120 nm, chi'' peak near 380 Hz
R['enpuku'] = dict(dls_nm=120, f_meas=380, f_model=round(float(f_of_d(120e-9)), 0), d_from_peak_nm=round(float(d_of_f(380) * 1e9), 1))
# Chung et al. (2004): 210 Hz -> 120 Hz on binding; reported increase of the hydrodynamic radius about 10 nm
d0 = d_of_f(210); d1 = d_of_f(120)
R['chung'] = dict(f0=210, f1=120, shift_meas=round(100 * (1 - 120 / 210), 1), d0_nm=round(float(d0 * 1e9), 1), d1_nm=round(float(d1 * 1e9), 1),
                  layer_inferred_nm=round(float((d1 - d0) / 2 * 1e9), 1), shift_model_10nm=round(float(100 * (1 - (d0 / (d0 + 20e-9))**3)), 1))
# Zhong et al. (2021, ACS Sens.): tau_B 0.742 -> 1.292 ms for functionalized BNF-80 with increasing mimic-virus concentration
d0 = d_of_tau(0.742e-3); d1 = d_of_tau(1.292e-3)
R['zhong_acs'] = dict(tau0_ms=0.742, tau1_ms=1.292, ratio=round(1.292 / 0.742, 2), shift_meas=round(100 * (1 - 0.742 / 1.292), 1),
                      d0_nm=round(float(d0 * 1e9), 1), d1_nm=round(float(d1 * 1e9), 1), layer_inferred_nm=round(float((d1 - d0) / 2 * 1e9), 1),
                      k_equiv=round(float((1.292 / 0.742)**(1.8 / 3)), 2))
# Zhong et al. (2021, JMMM): tau_B 0.85 -> 1.20 ms at 5 nM IgG when the particle concentration falls from 6.67 to 0.67 nM
d0 = d_of_tau(0.85e-3); d1 = d_of_tau(1.20e-3)
R['zhong_jmmm'] = dict(tau_high_c_ms=0.85, tau_low_c_ms=1.20, ratio=round(1.2 / 0.85, 2), d_high_c_nm=round(float(d0 * 1e9), 1), d_low_c_nm=round(float(d1 * 1e9), 1),
                       target_per_particle=[round(5 / 6.67, 2), round(5 / 0.67, 1)])
# model: full-binding shift of a 5 nm layer on the design of this work and on a 130 nm multicore particle
R['model'] = dict(shift_5nm_on_130nm=round(100 * (1 - (130 / 140)**3), 1), shift_5nm_on_40nm=round(100 * (1 - (40 / 50)**3), 1))
json.dump(R, open('compare_exp.json', 'w'), indent=1); print(json.dumps(R, indent=1))

# Design of magnetic nanoparticles for Brownian-relaxation biosensing

Code for the paper "Design of Magnetic Nanoparticles for Brownian-Relaxation Biosensing: Material, Size, and Temperature Constraints".
Requirements: Python 3, NumPy, Matplotlib. No external data are needed.

- `model3.py`      Neel/Brownian relaxation times, Debye AC susceptibility of a log-normal ensemble, peak-frequency finder, water viscosity.
- `analysis.py`    Usable size windows, binding shifts, polydispersity, temperature and plasma effects, Langmuir dose-response and LOD -> results.json
- `shell.py`       Shell-thickness trade-off -> shell.json (Figure 3b)
- `window_frac.py` Mass fraction of a log-normal sample inside the usable window -> window_frac.json (Figure 4a); run after analysis.py
- `extra3.py`      Thinnest Fe3O4 shell with a usable window (11.2 nm), amplitude ratio of the practical design (0.23), its mass inside the window -> extra3.json
- `figs3.py`       Figures 1-4

Run in order: python analysis.py ; python shell.py ; python window_frac.py ; python extra3.py ; python figs3.py   (about 1 min in total)
- `binding.py`     Target depletion (exact 1:1 site binding), Smoluchowski bridging aggregation with a hook effect, peak-shift and
                   single-frequency (3 f0) read-outs, LOD versus particle concentration -> binding.json and Figure 5 (about 1 min)

Additional requirement for `binding.py`: none beyond NumPy and Matplotlib. Run after the scripts above: python binding.py
- `sensitivity_fe3o4.py` Fe3O4 usable window versus effective anisotropy (8-30 kJ m-3) and chaining threshold (0.5-2),
                   thicker shells, CoFe2O4 reference, and multicore-particle estimates -> sensitivity_fe3o4.json (Table 3, Section 3.3)

Revised analysis (submitted version; the practical CoFe2O4 design is now a 16 nm median):
- `noise_lod.py`   Limit of detection with a noise that depends on signal amplitude and particle concentration, CoFe2O4 window for K = 40-250 kJ/m3 (cubic barrier K1 V/4) -> noise_lod.json
- `compare_exp.py` Brownian peak and binding shift compared with published AC-susceptibility data -> compare_exp.json
- `fig_lod.py`     Figure 6 from noise_lod.json. Run after noise_lod.py; `binding.py` imports constants from noise_lod.py.
- `mps_fpe.py`     Harmonic read-out: Brownian Fokker-Planck equation in Legendre form, test against the measured harmonic ratio of Zhong et al. (2021), comparison of the designs -> mps_fpe.json (Sections 2.8 and 3.9)

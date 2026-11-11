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
- `ms_sensitivity.py` Usable window of Fe3O4 for reduced saturation magnetization (bulk 480, 390, and 350 kA m-1; Section 3.1),
                   with the lowest anisotropy giving a usable window and the change of the binding signal; CoFe2O4 at 300 kA m-1
                   -> ms_sensitivity.json (a few seconds): python ms_sensitivity.py


## Additions of the revised manuscript

| File | Content |
|---|---|
| `noise_lod.py` | Signal-dependent read-out noise (Eq. (6)), limits of detection versus particle concentration, windows for the cubic barrier *K*1/4 -> `noise_lod.json` (Sections 2.7, 3.1, 3.7; Table 4) |
| `fig_lod.py` | Figure 6 from `noise_lod.json` |
| `compare_exp.py` | Comparison of the Brownian relation and the binding model with published AC-susceptibility data -> `compare_exp.json` (Section 3.8; Table 5) |
| `literature_data.csv` | Published values used in Section 3.8 and Table 5 |
| `binding.py`, `figs3.py` | Label and legend positions of Figures 3a and 5a changed (no change in the data) |

Run: `python noise_lod.py; python fig_lod.py; python compare_exp.py` (numpy, scipy, matplotlib).

Software used for the published results: Python 3.13, NumPy 2.5, SciPy 1.18, Matplotlib 3.11.

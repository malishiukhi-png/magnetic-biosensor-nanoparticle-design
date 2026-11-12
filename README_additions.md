# Additions for the 16 nm CoFe2O4 design (Table 3 and Section 3.7)

These scripts repeat the reference calculations of `analysis.py` and `binding.py` for the recommended
16 nm CoFe2O4 design (median 16 nm, sigma = 0.1, 12 nm shell). They change only the median size
(`DP = 16.0` and `D0 = 16e-9`) and the output file names; everything else is identical to the archived scripts.

| File | Content |
|---|---|
| `analysis16.py` | Copy of `analysis.py` with `DP = 16.0` -> `results16.json` (peak frequency, full-binding shift, temperature coefficient, plasma effect, LOD with target excess) |
| `binding16.py` | Copy of `binding.py` with `D0 = 16e-9` (figure not saved) -> `binding16.json` (LOD with depletion, bridging, hook effect) |
| `runK.py` | Runs either script with a different CoFe2O4 anisotropy, e.g. `OUT=results16_K45.json python runK.py 45e3 analysis16.py` |
| `*_K45.json` | Same calculations for K = 45 kJ m-3 (cubic barrier K1/4) |

Copy these files into the repository folder (next to `model3.py`) and run:

```
python analysis16.py
python binding16.py
OUT=results16_K45.json python runK.py 45e3 analysis16.py
OUT=binding16_K45.json python runK.py 45e3 binding16.py
```

Values used in the manuscript (K = 180 kJ m-3; K = 45 kJ m-3 differs by less than 1%):
f_p = 6.8 kHz, full-binding shift 48.0%, LOD with target excess (0.5% noise) 21.1 pM,
LOD with depletion at 10 ug/mL 0.60 nM, bridging LOD (alpha = 1e-3, 10 min, 10 ug/mL, 3 f0) 51 pM,
temperature coefficient 2.64% per degC, plasma effect 37-47%, hook maximum at about 18 nM.
Lambda (0.32), chi'' amplitude (0.55) and mass inside the window (98.3% / 97.8%) are taken from the
archived `noise_lod.json`.

Checked with Python 3.13 / NumPy 2.5 / SciPy 1.18 / Matplotlib 3.11; running the unmodified `analysis.py`
and `binding.py` reproduces the archived `results.json` and `binding.json` exactly.

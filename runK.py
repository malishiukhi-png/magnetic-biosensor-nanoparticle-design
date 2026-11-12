"""Run a script with a modified CoFe2O4 anisotropy: python runK.py <K_J_per_m3> <script>"""
import sys, runpy, model3
model3.MAT['CoFe2O4']['K'] = float(sys.argv[1])
runpy.run_path(sys.argv[2], run_name='__main__')

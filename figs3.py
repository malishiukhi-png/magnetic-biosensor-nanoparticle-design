import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 10, 'axes.labelweight': 'bold', 'axes.linewidth': 0.8,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
                     'xtick.minor.visible': True, 'ytick.minor.visible': True, 'legend.frameon': False, 'legend.fontsize': 8.5,
                     'savefig.dpi': 300, 'lines.linewidth': 1.6})
C = {'Fe3O4': '#2a78d6', 'CoFe2O4': '#eb6834'}; LAB = {'Fe3O4': 'Fe$_3$O$_4$', 'CoFe2O4': 'CoFe$_2$O$_4$'}
INK = '#0b0b0b'; MUTED = '#898781'; GREEN = '#1baf7a'; CM = 1 / 2.54
R = json.load(open('results.json')); S = json.load(open('shell.json')); WF = json.load(open('window_frac.json'))
def tag(ax, t, x=0.03, y=0.96): ax.text(x, y, t, transform=ax.transAxes, fontweight='bold', va='top', ha='left')

# ---------- Fig 1: concept and workflow ----------
plt.rcParams.update({'font.size': 8.5})
fig = plt.figure(figsize=(16 * CM, 6.8 * CM)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 160); ax.set_ylim(0, 68); ax.axis('off'); ax.set_aspect('equal')
def particle(x, y, bound):
    ax.add_patch(Circle((x, y), 9.5 if bound else 7.0, fc='#fdeee7' if bound else 'white', ec='#d95926' if bound else MUTED, lw=1, ls='--' if bound else '-'))
    ax.add_patch(Circle((x, y), 7.0, fc='#e8f1fc', ec=C['Fe3O4'], lw=1)); ax.add_patch(Circle((x, y), 3.2, fc='#52514e', ec='none'))
    for a in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        ax.plot([x + 3.2 * np.cos(a), x + 7.0 * np.cos(a)], [y + 3.2 * np.sin(a), y + 7.0 * np.sin(a)], color=C['Fe3O4'], lw=0.8)
        if bound: ax.add_patch(Circle((x + 8.3 * np.cos(a), y + 8.3 * np.sin(a)), 1.2, fc='#d95926', ec='none'))
particle(14, 42, False); particle(48, 42, True)
ax.annotate('', xy=(36, 42), xytext=(24, 42), arrowprops=dict(arrowstyle='-|>', color=INK, lw=1.1, mutation_scale=9))
ax.text(30, 46.5, '+ target', ha='center', fontsize=7.5, style='italic')
ax.text(14, 30, 'core d, shell δ\n$d_{h}$ = d + 2δ', ha='center', va='top', fontsize=7.5)
ax.text(48, 30, '$d_{h}$ + 2Δθ\nslower $τ_B$', ha='center', va='top', fontsize=7.5)
# spectrum sketch
f = np.logspace(2, 6, 200); x0 = 70; w = 36
def curve(fp, c, ls):
    y = (f / fp) / (1 + (f / fp) ** 2); xx = x0 + (np.log10(f) - 2) / 4 * w; ax.plot(xx, 30 + 26 * y / 0.5, color=c, ls=ls, lw=1.3)
curve(1e4, C['Fe3O4'], '-'); curve(5e3, '#d95926', '--')
ax.plot([x0, x0 + w], [30, 30], color=INK, lw=0.8); ax.plot([x0, x0], [30, 58], color=INK, lw=0.8)
ax.text(x0 + w / 2, 27, 'log frequency', ha='center', va='top', fontsize=7.5); ax.text(x0 - 2, 44, "χ''", ha='right', va='center', fontsize=8.5, style='italic')
ax.text(x0 + w / 2, 62, 'Peak shifts to lower f', ha='center', fontsize=7.8, fontweight='bold')
# constraints
bx = 113
for i, (t, b) in enumerate((('Brownian dominance', '$τ_N$ > 9 $τ_B$'), ('No chaining', 'dipolar λ < 1'), ('Polydispersity', 'mass in window'), ('Read-out drift', 'temperature, matrix'))):
    y = 54 - i * 13.5
    ax.add_patch(FancyBboxPatch((bx, y), 45, 11, boxstyle='round,pad=0,rounding_size=1.5', lw=1, ec='#52514e', fc='#f3f2ef'))
    ax.text(bx + 22.5, y + 7.6, t, ha='center', va='center', fontsize=7.8, fontweight='bold'); ax.text(bx + 22.5, y + 3.2, b, ha='center', va='center', fontsize=7.4)
ax.text(30, 10, 'Design variables: material (Fe$_3$O$_4$, CoFe$_2$O$_4$), core size d, shell δ, size dispersion σ\nOutputs: usable size window, signal shift, limit of detection',
        ha='left', va='center', fontsize=7.6, color='#52514e')
fig.savefig('fig3_concept.png'); plt.close(fig)
plt.rcParams.update({'font.size': 10})

# ---------- Fig 2: usable window ----------
fig, axs = plt.subplots(1, 2, figsize=(16 * CM, 6.6 * CM))
for mat in ('Fe3O4', 'CoFe2O4'):
    cv = R['curves'][mat]; d = np.array(cv['d'])
    axs[0].plot(d, cv['b'], color=C[mat], label=LAB[mat]); axs[1].semilogy(d, cv['lam'], color=C[mat], label=LAB[mat])
    w = R['window'][mat]['window']
    if w: axs[0].axvspan(w[0], w[1], color=C[mat], alpha=0.12, lw=0)
axs[0].axhline(0.9, color=MUTED, lw=1, ls='--'); wc = R['window']['CoFe2O4']['window']; wm = R['window']['Fe3O4']['window']
axs[0].text(13.5, 0.5, 'CoFe$_2$O$_4$ window\n%.1f–%.1f nm' % tuple(wc), ha='center', fontsize=7.5, color='#a33d14')
axs[0].annotate('Fe$_3$O$_4$ window\n%.1f–%.1f nm' % tuple(wm), (19.8, 0.95), (25.5, 0.45), fontsize=7.5, color='#1d5ea8', ha='center', arrowprops=dict(arrowstyle='-', color='#1d5ea8', lw=0.7)); axs[0].set_xlim(5, 30); axs[0].set_ylim(-0.02, 1.05)
axs[0].set_xlabel('Core diameter (nm)'); axs[0].set_ylabel('Brownian fraction τ$_N$/(τ$_N$+τ$_B$)'); axs[0].legend(loc='lower right'); tag(axs[0], '(a)', 0.03, 0.8)
axs[1].axhline(1, color=MUTED, lw=1, ls='--'); axs[1].set_xlim(5, 30); axs[1].set_ylim(1e-4, 1e2)
axs[1].set_xlabel('Core diameter (nm)'); axs[1].set_ylabel('Dipolar coupling λ'); axs[1].legend(loc='lower right'); tag(axs[1], '(b)')
fig.tight_layout(w_pad=1.0); fig.savefig('fig3_window.png'); plt.close(fig)

# ---------- Fig 3: spectra and shell trade-off ----------
fig, axs = plt.subplots(1, 2, figsize=(16 * CM, 6.6 * CM))
ax = axs[0]
for mat in ('Fe3O4', 'CoFe2O4'):
    sp = R['spectra'][mat]; f = np.array(sp['f']); a = np.array(sp['free']); b = np.array(sp['bound']); n = a.max()
    ax.semilogx(f, a / n, color=C[mat], label='%s, %.1f nm' % (LAB[mat], sp['d0'])); ax.semilogx(f, b / n, color=C[mat], ls='--')
ax.plot([], [], color=MUTED, label='Free'); ax.plot([], [], color=MUTED, ls='--', label='Fully bound')
ax.set_xlim(1e2, 1e6); ax.set_ylim(0, 1.1); ax.set_xlabel('Frequency (Hz)'); ax.set_ylabel("Normalized χ''"); ax.legend(loc='upper right', fontsize=7.5); tag(ax, '(a)')
ax = axs[1]
for mat in ('Fe3O4', 'CoFe2O4'):
    rows = [r for r in S[mat] if r['shift'] is not None]
    ax.plot([r['delta'] for r in rows], [r['shift'] for r in rows], 'o-', ms=3, color=C[mat], label=LAB[mat])
ax.axvspan(5, 11.2, color=C['Fe3O4'], alpha=0.08, lw=0); ax.text(8.2, 20, 'No usable\nFe$_3$O$_4$ size', ha='center', fontsize=7.5, color=C['Fe3O4'])
ax.set_xlim(5, 30); ax.set_ylim(15, 75); ax.set_xlabel('Shell thickness δ (nm)'); ax.set_ylabel('Peak shift, full binding (%)'); ax.legend(loc='upper right'); tag(ax, '(b)', 0.85, 0.12)
fig.tight_layout(w_pad=1.0); fig.savefig('fig3_spectra.png'); plt.close(fig)

# ---------- Fig 4: polydispersity, dose-response and drift ----------
fig, axs = plt.subplots(1, 2, figsize=(16 * CM, 6.6 * CM))
ax = axs[0]
for mat in ('Fe3O4', 'CoFe2O4'):
    ax.plot([r['sigma'] for r in WF[mat]], [r['frac'] for r in WF[mat]], 'o-', color=C[mat], ms=4, label=LAB[mat])
ax.plot([r['sigma'] for r in WF['CoFe2O4_K45']], [r['frac'] for r in WF['CoFe2O4_K45']], 's--', color=C['CoFe2O4'], ms=4, mfc='white', mew=1.1, label='CoFe$_2$O$_4$, $K$ = 45')
ax.set_xlim(0.03, 0.32); ax.set_ylim(0, 105); ax.set_xlabel('Size dispersion σ'); ax.set_ylabel('Mass in usable window (%)'); ax.legend(loc='center right'); tag(ax, '(a)', 0.03, 0.6)
ax = axs[1]
th = np.array(R['shift_vs_cover']['theta']); sh = np.array(R['shift_vs_cover']['shift'])
c = np.logspace(-3, 2, 300); theta = c / (1 + c); y = np.interp(theta, th, sh)
ax.loglog(c, np.maximum(y, 1e-3), color=C['CoFe2O4'], label='Binding signal')
ax.axhline(1.5, color=MUTED, ls='--', lw=1); ax.text(60, 1.5 * 0.82, '3 × 0.5% noise', fontsize=7.5, color='#52514e', ha='right', va='top')
ax.axhline(3.0, color=MUTED, ls=':', lw=1); ax.text(60, 3.0 * 1.12, '3 × 1% noise', fontsize=7.5, color='#52514e', ha='right', va='bottom')
dr = abs(R['temp_per_K']); ax.axhline(dr, color=GREEN, lw=1.2); ax.text(1.3e-3, dr * 1.12, '1 °C drift', fontsize=7.5, color='#127a55', va='bottom')
for noise, mk in (('0.5', 'o'), ('1.0', 's')):
    l = R['lod'][noise]['lod_over_Kd']; ax.plot([l], [3 * float(noise)], mk, color=INK, ms=4.5, mfc='white', mew=1.2)
ax.set_ylim(0.2, 80); ax.set_xlim(1e-3, 1e2); ax.set_xlabel('Target concentration / $K_d$'); ax.set_ylabel('Peak shift (%)'); tag(ax, '(b)', 0.04, 0.96)
fig.tight_layout(w_pad=1.0); fig.savefig('fig3_detect.png'); plt.close(fig)
print('ok')

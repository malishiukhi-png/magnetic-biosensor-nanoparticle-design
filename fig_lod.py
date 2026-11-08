"""Figure 6: limit of detection versus particle concentration with the concentration-dependent noise (noise_lod.json)."""
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 10, 'axes.labelweight': 'bold', 'axes.linewidth': 0.8,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True, 'xtick.minor.visible': True,
                     'ytick.minor.visible': True, 'legend.frameon': False, 'legend.fontsize': 8, 'savefig.dpi': 300, 'lines.linewidth': 1.8})
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']; INK = '#0b0b0b'; CM = 1 / 2.54
R = json.load(open('noise_lod.json'))['lod']
DES = (('CFO12_K180', 'CoFe$_2$O$_4$, 12 nm', C[0], '-', 'o'), ('CFO16_K180', 'CoFe$_2$O$_4$, 16 nm', C[1], '--', 's'),
       ('CFO18_K45', 'CoFe$_2$O$_4$, 18 nm', C[2], '-.', '^'), ('FO18_K20', 'Fe$_3$O$_4$, 18.5 nm ($K$ = 20)', C[3], ':', 'D'))
fig, axs = plt.subplots(1, 2, figsize=(16 * CM, 7.2 * CM))
for ax, mode, key, tag in ((axs[0], 'monolayer', 'lod_peak', '(a)'), (axs[1], 'bridging_1e-3', 'lod_flank', '(b)')):
    r0 = R['CFO16_K180'][mode]
    lo = np.array([np.nan if x[key] is None else x[key] * 1e9 for x in r0['1e-07']]); hi = np.array([np.nan if x[key] is None else x[key] * 1e9 for x in r0['1e-05']])
    cm = np.array([x['cm'] for x in r0['1e-06']])
    ok = ~(np.isnan(lo) | np.isnan(hi))
    ax.fill_between(cm[ok], lo[ok], hi[ok], color=C[1], alpha=0.14, lw=0, label='16 nm, $n_χ$ = 10$^{-7}$–10$^{-5}$')
    for name, lab, col, ls, mk in DES:
        rows = R[name][mode]['1e-06']
        x = np.array([r['cm'] for r in rows if r[key] is not None]); y = np.array([r[key] * 1e9 for r in rows if r[key] is not None])
        ax.loglog(x, y, color=col, ls=ls, marker=mk, ms=4, mfc='white', mew=1.1, label=lab)
    ax.set_xlim(0.8, 1200); ax.set_xlabel('Particle concentration (µg mL$^{-1}$)'); ax.set_ylabel('Limit of detection (nM)')
    ax.text(0.04, 0.95, tag, transform=ax.transAxes, fontweight='bold', va='top')
axs[0].set_ylim(0.1, 300); axs[1].set_ylim(0.02, 30)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, bbox_to_anchor=(0.5, 0.0), fontsize=7.8, columnspacing=1.2, handlelength=2.6)
fig.tight_layout(rect=(0, 0.13, 1, 1), w_pad=1.2); fig.savefig('fig_lod.png'); plt.close(fig)

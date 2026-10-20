import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
fig, axs = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for ax, f, t in [(axs[0], 'day_stats.json', '(a) Day (mean of Terra 10:30 and Aqua 13:30)'),
                 (axs[1], 'night_stats.json', '(b) Night (mean of Terra 22:30 and Aqua 01:30)')]:
    r = json.load(open(f)); y = np.array([x['year'] for x in r])
    b = np.array([x['coef_built'] for x in r]); s = np.array([x['coef_sand'] for x in r]); g = np.array([x['coef_veg'] for x in r])
    for v, lab, c in [(b - s, 'Built-up minus sand', '#e8643a'), (b, 'Built-up minus rock', '#7a6a58'), (g - b, 'Vegetation minus built-up', '#2a9d5c')]:
        ax.plot(y, v, 'o-', color=c, lw=1.8, ms=4, label=lab)
        sl, ic, _, _ = stats.theilslopes(v, y); ax.plot(y, ic + sl * y, '--', color=c, lw=1)
    ax.axhline(0, color='k', lw=0.8); ax.set_title(t, loc='left', fontsize=10.5); ax.set_xlabel('Year')
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', color='#cccccc', lw=0.6)
axs[0].set_ylabel('Summer LST contrast (°C)')
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=9, ncol=3, loc='lower center')
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig('/home/claude/najran_paper/figures/fig_modis_night.png', dpi=250)

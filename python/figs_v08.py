"""v0.8 figures: drivers (4 panels), trends (3 panels) with data-availability
flags on the annual vegetation series, and LULC/fate maps with geographic ticks."""
exec(open('/home/claude/analysis/figs_v03.py').read().split("if __name__")[0])
import numpy as np, pandas as pd
from scipy import stats

ST = '/home/claude/analysis/stats/'
veg = pd.read_csv(ST + 'veg_annual.csv').sort_values('year')
ep = pd.read_csv('/tmp/claude-0/-home-claude/40128d67-1783-5f2f-811e-6fee803b10b9/scratchpad/acc/Najran_epoch_stats.csv').sort_values('year')
rain = pd.read_csv(ST + 'rain_annual.csv'); rain = rain[rain.year <= 2025]
air = pd.read_csv('/home/claude/analysis/rain_airport_2001_2014_digitised.csv')
gr = pd.read_csv('/home/claude/analysis/grace_clean.csv').dropna()
gl = pd.read_csv('/home/claude/analysis/gldas_clean.csv').dropna()
gl['year'] = gl.date.str[:4].astype(int)
BLUE, ORANGE, GREY = '#2a78d6', '#e8643a', '#444444'
EVENTS = [(2006.15, '2006 market crash'), (2008.5, '2008 wheat phase-out')]

# availability adjustment: V ~ a + b*post2008 + c*months  (sensitivity only)
X = np.c_[np.ones(len(veg)), (veg.year > 2008).astype(float), veg.n_months_with_data]
coef = np.linalg.lstsq(X, veg.veg_km2, rcond=None)[0]
veg['adj'] = veg.veg_km2 - coef[2] * (veg.n_months_with_data - 12)
full = veg.n_months_with_data >= 10

def veg_panel(ax, legend=True):
    ax.plot(veg.year, veg.veg_km2, '-', color=BLUE, lw=1.8, zorder=2)
    ax.scatter(veg.year[full], veg.veg_km2[full], s=26, color=BLUE, zorder=3, label='annual, ≥10 months of data')
    ax.scatter(veg.year[~full], veg.veg_km2[~full], s=30, facecolor='white', edgecolor=BLUE, lw=1.4, zorder=3,
               label='annual, <10 months (few scenes)')
    ax.plot(veg.year, veg.adj, ls=(0, (4, 2)), color='#7fb3e8', lw=1.3, label='adjusted for data availability')
    ax.plot(ep.year, ep.area_km2_veg, 's', color='#238b45', ms=6, mec='k', mew=0.5, zorder=4,
            label='epoch maps (3-yr composites)')
    for a, b in [(2000, 2008), (2009, 2026)]:
        s_ = veg[(veg.year >= a) & (veg.year <= b)]
        sl, ic, _, _ = stats.theilslopes(s_.veg_km2, s_.year)
        ax.plot([a, b], [ic + sl * a, ic + sl * b], '--', color=GREY, lw=1.1)
    ax.axvspan(2002.9, 2013.0, color='#bbbbbb', alpha=0.18, lw=0)
    ax.text(2008, 157, 'Landsat 7 SLC-off / sparse Landsat 5 era', ha='center', fontsize=7.5, color='#555555')
    ax.set_ylabel('Vegetated area (km²)'); ax.set_ylim(35, 165)
    if legend: ax.legend(fontsize=7.2, loc='upper right', frameon=False, ncol=1)

def style(ax):
    ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', color='#cccccc', lw=0.6)
    ax.set_xlim(1999.4, 2026.6)

def fig_drivers():
    fig, ax = plt.subplots(4, 1, figsize=(10, 13.2), sharex=True)
    ax[0].bar(rain.year, rain.rain_mm, color=BLUE, width=0.72, label='CHIRPS (study-area mean)')
    ax[0].plot(air.year, air.rain_fig_units / 10, 'o', color=ORANGE, label='Najran airport gauge')
    l = stats.linregress(rain.year, rain.rain_mm)
    ax[0].set_title(f'(a) Annual rainfall: no significant trend ({l.slope:+.1f} mm/yr, p = {l.pvalue:.2f})', loc='left')
    ax[0].set_ylabel('Rainfall (mm/yr)'); ax[0].legend(frameon=False, fontsize=8)
    ax[1].plot(gr.yr, gr.lwe_cm, color='#9cc3ee', lw=1)
    ax[1].plot(gr.yr, gr.lwe_cm.rolling(12, center=True, min_periods=6).mean(), color=BLUE, lw=2.2)
    l = stats.linregress(gr.yr, gr.lwe_cm); xx = np.array([gr.yr.min(), gr.yr.max()])
    ax[1].plot(xx, l.intercept + l.slope * xx, '--', color=GREY)
    ax[1].text(2016.5, -1.2, f'{l.slope:.2f} ± {l.stderr:.2f} cm/yr (p < 0.001)', fontsize=9, color=GREY)
    ax[1].set_title('(b) GRACE/GRACE-FO terrestrial water storage (includes pumping)', loc='left'); ax[1].set_ylabel('TWS anomaly (cm)')
    m = gl.groupby('year').gws_mm.min()
    ax[2].plot(gl.yr, gl.gws_mm, color='#9cc3ee', lw=1, label='monthly')
    ax[2].plot(m.index + 0.5, m.values, 'o-', color=BLUE, lw=2.2, ms=4, label='annual minimum')
    ax[2].set_title('(c) GLDAS-2.2 modelled groundwater storage', loc='left'); ax[2].set_ylabel('Groundwater storage (mm)')
    ax[2].legend(frameon=False, fontsize=8, loc='upper right'); ax[2].set_ylim(340, 415)
    veg_panel(ax[3])
    ax[3].set_title('(d) Vegetated area in Wadi Najran, with Pettitt segments (Sen fits)', loc='left'); ax[3].set_xlabel('Year')
    for a in ax:
        style(a)
        for x, t in EVENTS: a.axvline(x, color='#888888', ls=':', lw=1)
    for x, t in EVENTS: ax[0].text(x + (-0.15 if x < 2007 else 0.15), 205, t, fontsize=7.5, color='#555555',
                                   ha='right' if x < 2007 else 'left')
    ax[0].set_ylim(0, 215)
    fig.tight_layout(); fig.savefig(OUT + 'fig_drivers_rain_grace_veg.png', dpi=250); plt.close(fig)

def fig_trends():
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.0))
    ax[0].plot(ep.year, ep.urban_core_km2, 's-', color=ORANGE, lw=2.2, label='Urban core (≥25% built, 500 m)')
    ax[0].plot(ep.year, ep.area_km2_built, 'o-', color=BLUE, lw=2.2, label='Built-up pixels')
    ax[0].set_title('(a) Urban expansion', loc='left'); ax[0].set_ylabel('Area (km²)'); ax[0].legend(frameon=False, fontsize=8.5)
    ax[0].set_xticks(ep.year)
    veg_panel(ax[1], legend=True); ax[1].axvline(2008.5, color='#888888', ls=':', lw=1)
    ax[1].set_title('(b) Annual vegetated area', loc='left')
    yrs = ep.year.values; w = 0.27; x = np.arange(len(yrs))
    for k, (c, nm, col) in enumerate([('sand', 'Sand', '#8c8a85'), ('rock', 'Rock', ORANGE), ('veg', 'Vegetation', BLUE)]):
        ax[2].bar(x + (k - 1) * w, ep[f'lst_C_{c}'] - ep.lst_C_built, w, color=col, label=nm)
    ax[2].axhline(0, color='k', lw=0.8); ax[2].set_xticks(x); ax[2].set_xticklabels(yrs)
    ax[2].set_title('(c) Summer daytime LST contrast', loc='left'); ax[2].set_ylabel('LST minus built-up LST (°C)')
    ax[2].legend(frameon=False, fontsize=8.5, ncol=3, loc='lower center'); ax[2].set_ylim(-4.5, 4)
    for a in ax[:2]: style(a)
    ax[2].spines[['top', 'right']].set_visible(False); ax[2].grid(axis='y', color='#cccccc', lw=0.6)
    fig.tight_layout(); fig.savefig(OUT + 'fig_preliminary_trends.png', dpi=250); plt.close(fig)

def fig_maps():
    cm = ListedColormap(['#d7301f', '#41ab5d', '#f2e3bd', '#3182bd', '#8c7e6e'])
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))
    for k, (y, t) in enumerate([(2000, '(a) Land cover 2000'), (2026, '(b) Land cover 2026 with approved plans')]):
        ax[k].imshow(lc[y], cmap=cm, vmin=-0.5, vmax=4.5, extent=ext, interpolation='nearest')
        ax[k].set_title(t, loc='left')
    for p in plans:
        for r in p: ax[1].plot(r[:, 0], r[:, 1], color='k', lw=0.6)
    v0 = lc[2000] == 1; f = np.full(lc[2000].shape, np.nan)
    f[v0 & (lc[2026] == 1)] = 0; f[v0 & (lc[2026] == 0)] = 1; f[v0 & np.isin(lc[2026], [2, 4])] = 2
    bg = np.where(mask, 0.0, np.nan)
    ax[2].imshow(bg, cmap=ListedColormap(['#f3f1ec']), extent=ext)
    ax[2].imshow(f, cmap=ListedColormap(['#1a9850', '#d7301f', '#e6a532']), vmin=-0.5, vmax=2.5, extent=ext, interpolation='nearest')
    ax[2].set_title('(c) Fate of 2000 vegetation by 2026', loc='left')
    for a in ax: km_axes(a, (388.4, 457.6), (1921.6, 1957.3))
    a1 = [Patch(color=c, label=l) for c, l in zip(['#d7301f', '#41ab5d', '#f2e3bd', '#8c7e6e'], ['Built-up', 'Vegetation', 'Sand / bare soil', 'Rock / gravel'])]
    ar = [np.nansum(f == k) * 0.0009 for k in range(3)]
    a2 = [Patch(color=c, label=f'{l} ({v:.1f} km²)') for c, l, v in zip(['#1a9850', '#d7301f', '#e6a532'], ['Still vegetated', 'Built over', 'Abandoned to bare'], ar)]
    fig.legend(handles=a1, loc='lower left', bbox_to_anchor=(0.04, -0.02), ncol=4, frameon=False, fontsize=9)
    fig.legend(handles=a2, loc='lower right', bbox_to_anchor=(0.98, -0.02), ncol=3, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(OUT + 'fig_maps_lulc_fate.png', dpi=250); plt.close(fig)
    return ar

if __name__ == '__main__':
    print('adjust coef', coef.round(2))
    fig_drivers(); fig_trends(); print('fate', [round(a, 1) for a in fig_maps()])
    fig_study_area(); fig_growth(); fig_vacant(); fig_heat(); fig_atlas(); print('v03 maps done')

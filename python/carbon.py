"""Carbon model (Monte Carlo) and scenario/carbon figures for manuscript v0.4.
Rooftop PV uses built-up area per zone until building footprints arrive
(ROOF_SOURCE switches to the footprint CSV automatically when present)."""
import os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
exec(open('/home/claude/analysis/figs_v03.py').read().split("if __name__")[0])

N = 20000
rng = np.random.default_rng(7)
OUTF = '/home/claude/najran_paper/figures/'
m = pd.read_csv('/home/claude/analysis/scenario_metrics.csv').set_index('scenario')
meta = json.load(open('/home/claude/analysis/scenario_meta.json'))
z = pd.read_csv('/home/claude/zones_out/Najran_zone_atlas_indicators.csv').set_index('zone_id')

# ---------------------------------------------------------- sourced constants
YIELD = 1721            # kWh/kWp/yr, Najran University station (Awan et al. 2018)
GEF_DNA = 0.568         # t CO2/MWh, Saudi DNA 2021 grid emission factor
RES_REGION_GWH = 2158   # Najran region residential electricity 2023 (GASTAT via Argaam)
POP_REGION, POP_CITY = 592300, 381431   # census 2022
RES_CITY_GWH = RES_REGION_GWH * POP_CITY / POP_REGION
POP40 = 490000

# ---------------------------------------------------------- rooftop PV
fp = '/home/claude/analysis/Najran_buildings_per_zone.csv'
if os.path.exists(fp):
    b = pd.read_csv(fp).set_index('zone_id')
    roof = b.roof_m2.reindex(z.index).fillna(0).values / 1e6               # km2, measured
    roof_src = 'footprints'
    roof_draw = lambda: roof
else:
    built = z.built26_km2.values
    roof_src = 'built-up proxy'
    roof_draw = lambda: built * 0.517 * rng.uniform(0.20, 0.40)              # km2 roofs
usable = rng.uniform(0.30, 0.60, N)
kwp_m2 = rng.uniform(0.15, 0.20, N)
yld = rng.normal(YIELD, 0.05 * YIELD, N) * rng.uniform(0.90, 0.97, N)     # soiling/availability
gef = rng.triangular(0.40, 0.55, 0.60, N)                                  # lifetime mean, decarbonising grid
res_city = RES_CITY_GWH * rng.uniform(0.9, 1.1, N)
roofs = np.array([roof_draw().sum() for _ in range(N)])                    # km2
cap_gw = roofs * 1e6 * usable * kwp_m2 / 1e6                                # GWp
gen_gwh = cap_gw * yld                                                      # GWh/yr (GWp*kWh/kWp = GWh)
co2_pv = gen_gwh * gef                                                      # kt/yr  (GWh * t/MWh = kt)
share = gen_gwh / res_city

# per-zone central estimate
zone_cap_mw = (z.built26_km2 * 0.517 * 0.30 if roof_src != 'footprints' else pd.Series(roof, z.index)) * 1e6 * 0.45 * 0.175 / 1e3
z['pv_mw_central'] = zone_cap_mw
z['pv_gwh_central'] = zone_cap_mw * YIELD * 0.935 / 1e3

# ---------------------------------------------------------- transport (2040, BAU vs Local Plan)
e_dens = rng.uniform(-0.12, -0.04, N)        # NRC 2009 range; Ewing & Cervero 2010 = -0.04
e_dist = rng.uniform(0.10, 0.30, N)          # distance to centre, Ewing & Cervero 2010 = 0.22
pc_t = rng.uniform(1.5, 2.5, N)              # t CO2 per capita, passenger travel
dln = (e_dens * np.log(m.loc['LocalPlan', 'density_pph'] / m.loc['BAU', 'density_pph']) +
       e_dist * np.log(m.loc['LocalPlan', 'mean_dist_all_built_km'] / m.loc['BAU', 'mean_dist_all_built_km']))
co2_tr = pc_t * POP40 * (1 - np.exp(dln)) / 1e3                             # kt/yr avoided by LP
vkt_red = 1 - np.exp(dln)

def q(a): return np.percentile(a, [5, 50, 95])
res = dict(roof_source=roof_src, roof_km2=q(roofs).tolist(), cap_gw=q(cap_gw).tolist(),
           gen_gwh=q(gen_gwh).tolist(), share_res=q(share).tolist(), co2_pv_kt=q(co2_pv).tolist(),
           res_city_gwh=RES_CITY_GWH, res_city_co2_kt=RES_CITY_GWH * GEF_DNA,
           vkt_reduction=q(vkt_red).tolist(), co2_transport_kt=q(co2_tr).tolist(),
           lp_target_mw=440, lp_target_gwh=0.44 * YIELD * 0.935)
json.dump(res, open('/home/claude/analysis/carbon_results.json', 'w'), indent=1)
z.round(3).to_csv('/home/claude/zones_out/Najran_zone_atlas_indicators.csv')
for k, v in res.items(): print(k, np.round(v, 3) if isinstance(v, list) else v)

# ---------------------------------------------------------- figure: scenarios 2040
S = np.load('/home/claude/analysis/scenarios_2040.npz')
c26 = S['c26']
fig = plt.figure(figsize=(13, 8.2))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.78])
titles = [('BAU', '(a) Trend (BAU)'), ('BND', '(b) Growth boundary only'), ('LP', '(c) Local Plan: boundary + density')]
cm = ListedColormap(['#ebe7e0', '#bdbdbd', '#d7301f', '#41ab5d'])
for k, (key, t) in enumerate(titles):
    a = fig.add_subplot(gs[0, k])
    sc = S[key]
    img = np.full(c26.shape, np.nan)
    img[mask] = 0
    img[(c26 == 0) & mask] = 1
    img[(sc == 0) & (c26 != 0) & mask] = 2
    img[(sc == 1) & mask] = 3
    a.imshow(img, cmap=cm, vmin=-0.5, vmax=3.5, extent=ext, interpolation='nearest')
    outline(a, list(zones.values()), color='k', lw=0.2)
    km_axes(a, (398, 457), (1922, 1957)); a.set_title(t, loc='left'); a.set_xlabel(''); a.set_ylabel('')
    if k == 0:
        a.legend(handles=[Patch(color=c, label=l) for c, l in zip(cm.colors[1:], ['built-up 2026', 'new 2026–2040', 'vegetation'])],
                 loc='lower right', fontsize=7)
sn = ['BAU', 'Boundary', 'LocalPlan']; lab = ['Trend', 'Boundary', 'Local Plan']; col = ['#d7301f', '#fdae61', '#2b8cbe']
mets = [('new_km2', 'New built-up land\n2026–2040 (km²)'), ('farm_lost_km2', 'Farmland built over (km²)'),
        ('density_pph', 'Gross density 2040\n(persons per ha)')]
for k, (mc, t) in enumerate(mets):
    a = fig.add_subplot(gs[1, k])
    v = m.loc[sn, mc].values
    a.bar(lab, v, color=col, edgecolor='k', lw=0.4)
    for i, x in enumerate(v): a.text(i, x * 1.02 + 0.3, f'{x:.1f}', ha='center', fontsize=8.5)
    a.set_title(t, fontsize=9.5); a.spines[['top', 'right']].set_visible(False)
    a.set_ylim(0, max(v) * 1.2 + 1)
fig.tight_layout(); fig.savefig(OUTF + 'fig_scenarios_2040.png', dpi=250); plt.close(fig)

# ---------------------------------------------------------- figure: carbon
fig, ax = plt.subplots(1, 3, figsize=(14, 4.3), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
top = z.sort_values('pv_gwh_central', ascending=False).head(15)
ax[0].barh([f'Zone {i}' for i in top.index][::-1], top.pv_gwh_central[::-1], color='#fe9929', edgecolor='k', lw=0.4)
ax[0].set_xlabel('Rooftop PV generation from footprints (GWh/yr)'); ax[0].set_title('(a) Fifteen zones with the largest potential', loc='left')
ax[0].spines[['top', 'right']].set_visible(False); ax[0].tick_params(labelsize=7.5)
lo, md, hi = q(gen_gwh)
bars = [('Rooftop PV\n(median, 90% range)', md, lo, hi, '#fe9929'),
        ('City residential\ndemand 2023', RES_CITY_GWH, None, None, '#636363'),
        ('Local Plan solar\ntarget (440 MW)', res['lp_target_gwh'], None, None, '#2b8cbe')]
for i, (l, v, a_, b_, c) in enumerate(bars):
    ax[1].bar(i, v, color=c, edgecolor='k', lw=0.4)
    if a_ is not None: ax[1].errorbar(i, v, yerr=[[v - a_], [b_ - v]], color='k', capsize=5)
    ax[1].text(i, (b_ or v) * 1.03, f'{v:,.0f}', ha='center', fontsize=8.5)
ax[1].set_xticks(range(3)); ax[1].set_xticklabels([b[0] for b in bars], fontsize=7.5)
ax[1].set_ylabel('GWh per year'); ax[1].set_title('(b) Rooftop PV vs demand and target', loc='left')
ax[1].spines[['top', 'right']].set_visible(False)
items = [('Rooftop PV\n(avoided)', q(co2_pv), '#fe9929'),
         ('Compact growth,\ntransport 2040 (avoided)', q(co2_tr), '#2b8cbe'),
         ('City residential\nelectricity (emitted)', [res['res_city_co2_kt']] * 3, '#636363')]
for i, (l, (lo_, md_, hi_), c) in enumerate(items):
    ax[2].barh(i, md_, color=c, edgecolor='k', lw=0.4)
    if hi_ > lo_: ax[2].errorbar(md_, i, xerr=[[md_ - lo_], [hi_ - md_]], color='k', capsize=4)
    ax[2].text(hi_ * 1.15, i, f'{md_:,.0f}', va='center', fontsize=8.5)
ax[2].set_xscale('log'); ax[2].set_xlim(5, 2e4); ax[2].set_yticks(range(3)); ax[2].set_yticklabels([x[0] for x in items], fontsize=7.5)
ax[2].invert_yaxis(); ax[2].set_xlabel('kt CO$_2$ per year (log scale; median and 90% range)')
ax[2].set_title('(c) Avoided and current emissions', loc='left'); ax[2].spines[['top', 'right']].set_visible(False)
fig.tight_layout(); fig.savefig(OUTF + 'fig_carbon.png', dpi=250); plt.close(fig)

"""New figures for manuscript v0.3 (study area, growth by year, Sankey,
vacant parcels, relative LST + NDVI trend, citizen zone atlas)."""
import sys, json, csv
sys.path.insert(0, '/home/claude/analysis')
from load import *
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Polygon as MPoly, Patch, PathPatch
from matplotlib.path import Path as MPath
from matplotlib.collections import PatchCollection
from scipy import ndimage

OUT = '/home/claude/najran_paper/figures/'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.titlesize': 10, 'axes.titleweight': 'bold'})
L, T, tr, mask, gt = load_all()
H, W = mask.shape
ext = [gt[0] / 1e3, (gt[0] + W * gt[1]) / 1e3, (gt[2] + H * gt[3]) / 1e3, gt[2] / 1e3]
lc = {y: np.where(mask, L[y][0].astype(float), np.nan) for y in YEARS}

# ---------- vector layers in UTM km
zj = json.load(open('/home/claude/zones_out/Najran_residential_zones.geojson'))
zones = {}
for f in zj['features']:
    ring = [s.wgs84_to_utm(lo, la) for lo, la in f['geometry']['coordinates'][0]]
    zones[f['properties']['zone_id']] = np.array(ring) / 1e3
_, pg = s.read_shp(SHP + 'plans.shp')
_, prow = s.read_dbf(SHP + 'plans.dbf')
plans = [[np.array(p) / 1e3 for p in g] for g in pg]
zs = pd.read_csv('/home/claude/zones_out/Najran_zones_stats_2000_2026.csv').set_index('zone_id')

def outline(ax, polys, **kw):
    for p in polys:
        ax.plot(p[:, 0], p[:, 1], **kw)

def geo_ticks(ax, step=0.1):
    """Label UTM-km axes with geographic coordinates (deg E / deg N)."""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    if x1 - x0 < 25: step = 0.05
    lon_a, lat_a = s.utm_to_wgs84(x0 * 1e3, y0 * 1e3); lon_b, lat_b = s.utm_to_wgs84(x1 * 1e3, y1 * 1e3)
    latm = (lat_a + lat_b) / 2; lonm = (lon_a + lon_b) / 2
    lons = np.arange(np.ceil(lon_a / step) * step, lon_b, step)
    lats = np.arange(np.ceil(lat_a / (step / 2)) * (step / 2), lat_b, step / 2)
    xt = [s.wgs84_to_utm(lo, latm)[0] / 1e3 for lo in lons]
    yt = [s.wgs84_to_utm(lonm, la)[1] / 1e3 for la in lats]
    ax.set_xticks(xt); ax.set_xticklabels([(f'{lo:.2f}°E' if step < 0.1 else f'{lo:.1f}°E') for lo in lons])
    ax.set_yticks(yt); ax.set_yticklabels([f'{la:.2f}°N' for la in lats])
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)

def km_axes(ax, xl=None, yl=None):
    if xl: ax.set_xlim(xl)
    if yl: ax.set_ylim(yl)
    ax.set_aspect('equal'); ax.tick_params(labelsize=7)
    geo_ticks(ax)
    ax.set_xlabel('Longitude', fontsize=7.5)
    ax.set_ylabel('Latitude', fontsize=7.5)

def scalebar(ax, x, y, km=5):
    ax.plot([x, x + km], [y, y], 'k-', lw=3, solid_capstyle='butt')
    ax.text(x + km / 2, y + 0.6, f'{km} km', ha='center', fontsize=7)

def north(ax, x, y):
    ax.annotate('N', xy=(x, y), xytext=(x, y - 3), ha='center', fontsize=9, weight='bold',
                arrowprops=dict(arrowstyle='-|>', lw=1.2, color='k'))

CITY = ([398, 457], [1922, 1957])
LCC = ['#d7301f', '#41ab5d', '#f2e3bd', '#3182bd', '#9c8f80']

# ====================================================== 1. study area
def fig_study_area():
    fig, ax = plt.subplots(figsize=(10, 5.6))
    cm = ListedColormap(['#e8a29a', '#7fc58d', '#f5ecd6', '#9ecae1', '#c9c1b6'])
    ax.imshow(lc[2026], cmap=cm, vmin=-0.5, vmax=4.5, extent=ext, interpolation='nearest')
    zpolys = list(zones.values())
    outline(ax, zpolys, color='#333333', lw=0.5)
    for p in plans:
        for r in p:
            ax.fill(r[:, 0], r[:, 1], facecolor='none', edgecolor='#6a3d9a', lw=0.9, hatch='////', alpha=0.9)
    # AOI frame and labels (well-known locations only)
    ax.text(412, 1931.2, 'W a d i   N a j r a n', rotation=24, fontsize=9, style='italic', color='#00441b')
    e, n = s.wgs84_to_utm(44.419, 17.611)
    ax.plot(e / 1e3, n / 1e3, marker='^', ms=8, color='k')
    ax.text(e / 1e3 + 0.8, n / 1e3 + 0.5, 'Najran Airport', fontsize=8)
    ax.annotate('to Rub\' al Khali', xy=(457.4, 1938), xytext=(445, 1929), fontsize=8,
                arrowprops=dict(arrowstyle='->', color='k'))
    ax.annotate('granitic uplands', xy=(392, 1950), fontsize=8, color='#4d4d4d')
    scalebar(ax, 392, 1924); north(ax, 455, 1956)
    km_axes(ax, (388.4, 457.6), (1921.6, 1957.3))
    h = [Patch(color=c, label=l) for c, l in zip(['#e8a29a', '#7fc58d', '#f5ecd6', '#c9c1b6'],
         ['Built-up 2026', 'Vegetation 2026', 'Sand / bare soil', 'Rock / gravel'])]
    h += [plt.Line2D([], [], color='#333333', lw=0.8, label='Residential zones (63)'),
          Patch(facecolor='none', edgecolor='#6a3d9a', hatch='////', label='Approved plans (east)')]
    ax.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, -0.11), fontsize=7.5, frameon=False, ncol=6)
    fig.savefig(OUT + 'fig_study_area.png', dpi=250, bbox_inches='tight'); plt.close(fig)

# ====================================================== 2. growth by year built
def fig_growth():
    first = np.full((H, W), np.nan)
    for i, y in enumerate(YEARS[::-1]):
        first = np.where(lc[y] == 0, i, first)        # earliest epoch overwrites
    first = np.where(np.isnan(first), np.nan, 5 - first)  # 0=2000 ... 5=2026
    cols = ['#3f007d', '#6a51a3', '#2171b5', '#41ab5d', '#fdae61', '#d7191c']
    labs = ['before 2000', '2000–05', '2005–10', '2010–15', '2015–20', '2020–26']
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[2.3, 1]))
    bg = np.where(np.isnan(lc[2026]), np.nan, np.where(lc[2026] == 1, 1, 0))
    ax[0].imshow(bg, cmap=ListedColormap(['#f3f1ec', '#cfe8cf']), extent=ext, interpolation='nearest')
    ax[0].imshow(first, cmap=ListedColormap(cols), vmin=-0.5, vmax=5.5, extent=ext, interpolation='nearest')
    outline(ax[0], list(zones.values()), color='#555555', lw=0.3)
    scalebar(ax[0], 400, 1924); north(ax[0], 455, 1956)
    km_axes(ax[0], *CITY)
    ax[0].legend(handles=[Patch(color=c, label=l) for c, l in zip(cols, labs)] +
                 [Patch(color='#cfe8cf', label='vegetation 2026')], loc='lower right', fontsize=7.5,
                 title='Built-up first mapped', title_fontsize=7.5, framealpha=0.9)
    ax[0].set_title('(a) Where and when Najran grew')
    area = [np.nansum(first == k) * 0.0009 for k in range(6)]
    yrs = [5, 5, 5, 5, 6]
    rate = [a / n for a, n in zip(area[1:], yrs)]
    ax[1].bar(labs[1:], rate, color=cols[1:], edgecolor='k', lw=0.4)
    for i, v in enumerate(rate):
        ax[1].text(i, v + 0.2, f'{v:.1f}', ha='center', fontsize=8)
    ax[1].set_ylabel('New built-up land (km² per year)')
    ax[1].tick_params(axis='x', rotation=30, labelsize=7.5)
    ax[1].spines[['top', 'right']].set_visible(False)
    ax[1].set_title('(b) Annual rate of new built-up land')
    fig.tight_layout(); fig.savefig(OUT + 'fig_growth_by_epoch.png', dpi=250); plt.close(fig)
    return area

_DOM = None
def urban_domain():
    global _DOM
    if _DOM is None:
        m = np.zeros((H, W), bool)
        for z in zones.values(): m |= rasterize([z * 1e3], gt, (H, W))
        for p in plans: m |= rasterize([r * 1e3 for r in p], gt, (H, W))
        _DOM = m & mask
    return _DOM

def zone_masks():
    return {z: rasterize([p * 1e3], gt, (H, W)) & mask for z, p in zones.items()}

# ====================================================== 3. Sankey 2000 -> 2026
def fig_sankey():
    names = ['Built-up', 'Vegetation', 'Sand / soil', 'Rock / gravel']
    cls = [0, 1, 2, 4]
    dom = urban_domain()
    a, b = np.where(dom, lc[2000], np.nan), np.where(dom, lc[2026], np.nan)
    M = np.array([[np.nansum((a == i) & (b == j)) * 0.0009 for j in cls] for i in cls])
    fig, ax = plt.subplots(figsize=(8, 5.2))
    gap = 0.03; scale = 1 / M.sum() * (1 - gap * 3)
    xL, xR, w = 0.0, 1.0, 0.06
    def stack(tot):
        pos, y = [], 1.0
        for t in tot:
            h = t * scale; pos.append([y - h, y]); y -= h + gap
        return pos
    pl, pr = stack(M.sum(1)), stack(M.sum(0))
    offL = [p[1] for p in pl]; offR = [p[1] for p in pr]
    # draw small flows first, keep order by source then target
    for i in range(4):
        for j in range(4):
            h = M[i, j] * scale
            if h <= 0: continue
            y0t, y1t = offL[i], offR[j]
            offL[i] -= h; offR[j] -= h
            xs = [xL + w, 0.5, 0.5, xR]
            verts = [(xL + w, y0t), (0.5, y0t), (0.5, y1t), (xR, y1t), (xR, y1t - h),
                     (0.5, y1t - h), (0.5, y0t - h), (xL + w, y0t - h), (xL + w, y0t)]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
                     MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.CLOSEPOLY]
            alpha = 0.35 if i == j else 0.75
            ax.add_patch(PathPatch(MPath(verts, codes), fc=LCC[cls[i]], ec='none', alpha=alpha))
    for k in range(4):
        for x, p, tot in [(xL, pl[k], M.sum(1)[k]), (xR, pr[k], M.sum(0)[k])]:
            ax.add_patch(plt.Rectangle((x, p[0]), w, p[1] - p[0], fc=LCC[cls[k]], ec='k', lw=0.5))
        ax.text(xL - 0.01, np.mean(pl[k]), f'{names[k]}\n{M.sum(1)[k]:.0f} km²', ha='right', va='center', fontsize=8.5)
        ax.text(xR + w + 0.01, np.mean(pr[k]), f'{names[k]}\n{M.sum(0)[k]:.0f} km²', ha='left', va='center', fontsize=8.5)
    ax.text(xL + w / 2, 1.04, '2000', ha='center', weight='bold'); ax.text(xR + w / 2, 1.04, '2026', ha='center', weight='bold')
    ax.text(0.53, 1.04, f'into built-up: rock {M[3,0]:.0f}, sand {M[2,0]:.0f}, vegetation {M[1,0]:.0f} km²',
            ha='center', fontsize=8.5, color='#b2182b')
    ax.set_xlim(-0.25, 1.3); ax.set_ylim(-0.02, 1.08); ax.axis('off')
    ax.text(0.53, -0.04, f'Urban domain: 63 residential zones + approved plans ({np.sum(dom)*0.0009:.0f} km²). Light bands = unchanged.',
            ha='center', fontsize=7.5, color='#555555')
    fig.savefig(OUT + 'fig_sankey_2000_2026.png', dpi=250, bbox_inches='tight'); plt.close(fig)
    return M

# ====================================================== 4. vacant parcels
def fig_vacant():
    _, gpar = s.read_shp(SHP + 'parcels.shp')
    _, rpar = s.read_dbf(SHP + 'parcels.dbf')
    built = (lc[2026] == 0)
    bdil = ndimage.binary_dilation(built, np.ones((3, 3)))   # tolerate 1-pixel offset
    cx, cy, st = [], [], []
    for g in gpar:
        if not g: continue
        p = np.vstack(g); x, y = p.mean(0)
        c = int((x - gt[0]) / gt[1]); r = int((y - gt[2]) / gt[3])
        if 0 <= r < H and 0 <= c < W:
            cx.append(x / 1e3); cy.append(y / 1e3); st.append(2 if built[r, c] else (1 if bdil[r, c] else 0))
    cx, cy, st = map(np.array, (cx, cy, st))
    fig, ax = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw=dict(width_ratios=[1.7, 1]))
    bg = np.where(np.isnan(lc[2026]), np.nan, np.where(lc[2026] == 0, 1, 0))
    ax[0].imshow(bg, cmap=ListedColormap(['#f4f2ee', '#d9d9d9']), extent=ext, interpolation='nearest')
    for p in plans:
        for r in p: ax[0].plot(r[:, 0], r[:, 1], color='#6a3d9a', lw=0.7)
    for k, c, l in [(0, '#f16913', 'vacant'), (1, '#fdd0a2', 'within 30 m of built-up'), (2, '#08519c', 'built')]:
        m = st == k
        ax[0].scatter(cx[m], cy[m], s=0.6, c=c, label=f'{l} ({m.sum():,})', lw=0)
    km_axes(ax[0], (439, 454), (1940.5, 1956)); scalebar(ax[0], 440, 1941.5, 2)
    ax[0].legend(markerscale=12, fontsize=7.5, loc='lower right', title='Approved parcels, 2026', title_fontsize=7.5)
    ax[0].set_title('(a) Approved plots with and without buildings')
    # per-plan bars
    pnames = {'مخطط الدارة': 'Al-Dara', 'مخطط موقع الجامعة': 'University site', 'مخطط حي الجامعة': 'University district',
              'مخطط النهضة السكنى': 'Al-Nahda', 'مخطط شرق المطار - أ': 'East Airport A', 'مخطط شرق المطار - ب': 'East Airport B',
              'مخطط شرق المطار - ج': 'East Airport C', 'مخطط الروضة': 'Al-Rawdah', 'هئية الاسكان العام': 'Gen. Housing Auth.',
              'مخطط شرق نجران الجديد': 'New East Najran'}
    pid_name = {r['PlanNO']: pnames.get(r['PlansName']) for r in prow}
    _, rpar2 = s.read_dbf(SHP + 'parcels.dbf')
    keys = [r['PlanNO'] for g, r in zip(gpar, rpar2) if g]
    df = pd.DataFrame({'plan': [pid_name.get(k) for k in keys][:len(st)], 'st': st})
    t = df.dropna().groupby('plan').st.agg(n='size', vac=lambda v: (v == 0).sum())
    t['pct'] = 100 * t.vac / t.n; t = t.sort_values('pct')
    ax[1].barh(t.index, t.pct, color='#f16913', edgecolor='k', lw=0.4)
    for i, (p, n) in enumerate(zip(t.pct, t.n)):
        ax[1].text(p + 1, i, f'{p:.0f}% of {n:,}', va='center', fontsize=7.5)
    ax[1].set_xlim(0, 115); ax[1].set_xlabel('Parcels with no building detected (%)')
    ax[1].spines[['top', 'right']].set_visible(False); ax[1].tick_params(labelsize=7.5)
    ax[1].set_title('(b) Vacant share by plan')
    fig.tight_layout(); fig.savefig(OUT + 'fig_vacant_parcels.png', dpi=250); plt.close(fig)
    return st, t

# ====================================================== 5. relative LST and NDVI trend
def fig_heat():
    rel = {}
    for y in [2000, 2026]:
        nat = np.isin(lc[y], [2, 4])
        rel[y] = T[y] - np.nanmedian(T[y][nat])
    fig, axs = plt.subplots(2, 2, figsize=(12, 7.6))
    ax = axs.ravel()
    panels = [(rel[2000], 'RdYlBu_r', -6, 6, '(a) Relative summer LST, 2000', 'LST − desert median (°C)'),
              (rel[2026], 'RdYlBu_r', -6, 6, '(b) Relative summer LST, 2026', 'LST − desert median (°C)'),
              (rel[2026] - rel[2000], 'PuOr_r', -5, 5, '(c) Change in relative LST, 2000–2026', '°C (positive = warmer)'),
              (np.where(mask, tr[0], np.nan), 'BrBG', -0.06, 0.06, '(d) NDVI trend, 2000–2026', 'Sen slope (NDVI per decade)')]
    for a, (img, cm, v0, v1, t, lab) in zip(ax, panels):
        im = a.imshow(img, cmap=cm, vmin=v0, vmax=v1, extent=ext)
        outline(a, list(zones.values()), color='k', lw=0.25)
        km_axes(a, *CITY); a.set_title(t, loc='left')
        cb = plt.colorbar(im, ax=a, shrink=0.72, pad=0.01); cb.set_label(lab, fontsize=7.5); cb.ax.tick_params(labelsize=7)
    scalebar(ax[0], 400, 1924)
    fig.tight_layout(); fig.savefig(OUT + 'fig_heat_relative_ndvi_trend.png', dpi=250); plt.close(fig)
    return rel

# ====================================================== 6. citizen atlas (4 priorities per zone)
ROOF_FRAC, USABLE, KWP_M2, UA_BUILT = 0.30, 0.50, 0.15, 0.517   # assumptions, see text

def zone_table():
    zm = zone_masks()
    rel = {}
    for y in [2000, 2026]:
        nat = np.isin(lc[y], [2, 4]); rel[y] = T[y] - np.nanmedian(T[y][nat])
    _, gpar = s.read_shp(SHP + 'parcels.shp')
    built = lc[2026] == 0
    vac_xy = []
    for g in gpar:
        if not g: continue
        x, y = np.vstack(g).mean(0); c = int((x - gt[0]) / gt[1]); r = int((y - gt[2]) / gt[3])
        if 0 <= r < H and 0 <= c < W and not ndimage.binary_dilation(built[max(r-1,0):r+2, max(c-1,0):c+2]).any():
            vac_xy.append((x / 1e3, y / 1e3))
    vac_xy = np.array(vac_xy)
    rows = []
    for z, m in zm.items():
        n = m.sum(); a = n * 0.0009
        b26 = np.sum(m & (lc[2026] == 0)) * 0.0009
        v00 = np.sum(m & (lc[2000] == 1)) * 0.0009; v26 = np.sum(m & (lc[2026] == 1)) * 0.0009
        open26 = np.sum(m & np.isin(lc[2026], [2, 4])) * 0.0009
        inside = MPath(zones[z]).contains_points(vac_xy).sum() if len(vac_xy) else 0
        rows.append(dict(zone_id=z, area_km2=a, built26_km2=b26, veg00_km2=v00, veg26_km2=v26,
                         open26_km2=open26, vacant_parcels=inside,
                         lst26=np.nanmean(T[2026][m]), rel_lst26=np.nanmean(rel[2026][m]),
                         d_rel_lst=np.nanmean(rel[2026][m] - rel[2000][m]),
                         pv_mwp=b26 * UA_BUILT * 1e6 * ROOF_FRAC * USABLE * KWP_M2 / 1e3))
    d = pd.DataFrame(rows).set_index('zone_id')
    import os
    fpc = '/home/claude/analysis/Najran_buildings_per_zone.csv'
    if os.path.exists(fpc):   # measured footprints (Script 03)
        fp = pd.read_csv(fpc).set_index('zone_id')
        d['roof_km2'] = fp.roof_m2.reindex(d.index).fillna(0) / 1e6
        d['n_bldg'] = fp.n_bldg.reindex(d.index).fillna(0)
        d['pv_mwp'] = d.roof_km2 * 1e6 * 0.45 * 0.175 / 1e3
    d['veg26_pct'] = 100 * d.veg26_km2 / d.area_km2; d['built26_pct'] = 100 * d.built26_km2 / d.area_km2
    zsc = lambda v: (v - v.mean()) / v.std()
    d['heat_idx'] = zsc(d.rel_lst26) - zsc(d.veg26_pct) + 0.5 * zsc(d.built26_pct)
    d['heat_cls'] = pd.qcut(d.heat_idx, 3, labels=['Lower', 'Medium', 'High'])
    loss = d.veg00_km2 - d.veg26_km2
    d['farm_cls'] = np.select(
        [(d.veg26_pct >= 20) & (loss / d.veg00_km2.clip(lower=1e-6) >= 0.25),
         d.veg26_pct >= 20,
         (d.veg00_km2 / d.area_km2 >= 0.2)],
        ['Protect now', 'Keep and monitor', 'Largely lost'], 'Little farmland')
    return d, vac_xy

def choropleth(ax, d, col, cmap=None, cats=None, colors=None, vmin=None, vmax=None, label=''):
    ps = [MPoly(zones[z], closed=True) for z in d.index]
    if cats is not None:
        cmap_ = ListedColormap(colors); vals = d[col].map({c: i for i, c in enumerate(cats)}).astype(float)
        pc = PatchCollection(ps, cmap=cmap_, edgecolor='k', lw=0.35); pc.set_array(vals.values); pc.set_clim(-0.5, len(cats) - 0.5)
    else:
        pc = PatchCollection(ps, cmap=cmap, edgecolor='k', lw=0.35); pc.set_array(d[col].values); pc.set_clim(vmin, vmax)
    ax.add_collection(pc); return pc

def fig_atlas():
    d, vac = zone_table()
    d.round(3).to_csv('/home/claude/zones_out/Najran_zone_atlas_indicators.csv')
    bg = np.where(np.isnan(lc[2026]), np.nan, np.where(lc[2026] == 1, 1, 0))
    fig, axs = plt.subplots(2, 2, figsize=(13, 8.4)); ax = axs.ravel()
    for a in ax:
        a.imshow(bg, cmap=ListedColormap(['#f5f3ef', '#d6ead3']), extent=ext, interpolation='nearest')
    # (a) heat & greening
    choropleth(ax[0], d, 'heat_cls', cats=['Lower', 'Medium', 'High'], colors=['#fee8c8', '#fdbb84', '#d7301f'])
    ax[0].legend(handles=[Patch(color=c, label=l) for c, l in zip(['#d7301f', '#fdbb84', '#fee8c8'], ['High', 'Medium', 'Lower'])],
                 title='Greening & shade priority', loc='lower right', fontsize=7.5, title_fontsize=7.5)
    ax[0].set_title('(a) Heat and greenery: where to plant and shade first', loc='left')
    # (b) farmland
    cats = ['Protect now', 'Keep and monitor', 'Largely lost', 'Little farmland']
    cols = ['#00441b', '#74c476', '#bf812d', '#f0f0f0']
    choropleth(ax[1], d, 'farm_cls', cats=cats, colors=cols)
    ax[1].legend(handles=[Patch(color=c, label=f'{l} ({(d.farm_cls == l).sum()})') for c, l in zip(cols, cats)],
                 title='Wadi farmland status', loc='lower right', fontsize=7.5, title_fontsize=7.5)
    ax[1].set_title('(b) Farmland protection in Wadi Najran', loc='left')
    # (c) infill housing
    pc = choropleth(ax[2], d, 'open26_km2', cmap='Blues', vmin=0, vmax=np.percentile(d.open26_km2, 95))
    ax[2].scatter(vac[:, 0], vac[:, 1], s=0.4, c='#e6550d', lw=0, label=f'vacant approved parcel ({len(vac):,})')
    cb = plt.colorbar(pc, ax=ax[2], shrink=0.7, pad=0.01); cb.set_label('Unbuilt non-farm land in zone (km²)', fontsize=7.5)
    ax[2].legend(markerscale=14, loc='lower right', fontsize=7.5)
    ax[2].set_title('(c) Future housing: infill before new subdivisions', loc='left')
    # (d) rooftop solar
    pc = choropleth(ax[3], d, 'pv_mwp', cmap='YlOrBr', vmin=0, vmax=np.percentile(d.pv_mwp, 95))
    cb = plt.colorbar(pc, ax=ax[3], shrink=0.7, pad=0.01); cb.set_label('Rooftop PV potential from building footprints (MWp)', fontsize=7.5)
    ax[3].set_title('(d) Rooftop solar potential per zone', loc='left')
    for a in ax:
        km_axes(a, (398, 457), (1922, 1957)); a.set_xlabel(''); a.set_ylabel('')
    scalebar(ax[0], 400, 1924); north(ax[0], 455, 1956)
    fig.tight_layout(); fig.savefig(OUT + 'fig_citizen_atlas.png', dpi=250); plt.close(fig)
    return d, vac

if __name__ == '__main__':
    fig_study_area(); print('study area')
    a = fig_growth(); print('growth', [round(x, 1) for x in a])
    M = fig_sankey(); print('sankey\n', M.round(1))
    st, t = fig_vacant(); print('parcels', len(st), np.bincount(st)); print(t)
    fig_heat(); print('heat; trend band0 range', np.nanpercentile(tr[0], [2, 50, 98]))
    fig_atlas(); print('atlas')

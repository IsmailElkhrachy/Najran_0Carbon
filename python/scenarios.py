"""ML-CA growth scenarios for Najran to 2040.
Transition model: gradient boosting on non-built(t) -> built(t+1), trained on
2010->2015 and 2015->2020, validated on 2020->2026 (figure of merit).
Two scenarios to 2040: BAU (trend demand, unconstrained) and Local Plan
(population-based demand at target density, inside the urban domain, no
farmland). Outputs a raster npz, a CSV of metrics and a figure."""
import sys, json
sys.path.insert(0, '/home/claude/analysis')
import numpy as np, pandas as pd
from scipy import ndimage
from sklearn.ensemble import HistGradientBoostingClassifier
exec(open('/home/claude/analysis/figs_v03.py').read().split("if __name__")[0])

rng = np.random.default_rng(42)
PIX = 0.0009
dom = urban_domain()
plans_m = np.zeros((H, W), bool)
for p in plans: plans_m |= rasterize([r * 1e3 for r in p], gt, (H, W))
xs, ys = pixel_xy(gt, (H, W)); XX, YY = np.meshgrid(xs, ys)
core00 = lc[2000] == 0
cy, cx = np.argwhere(ndimage.binary_opening(core00, iterations=2)).mean(0)
# city centre = centroid of the largest connected built-up patch in 2000
lab, n = ndimage.label(core00)
big = lab == (np.bincount(lab.ravel())[1:].argmax() + 1)
cy, cx = np.argwhere(big).mean(0)
CX, CY = xs[int(cx)], ys[int(cy)]
dist_c = np.hypot(XX - CX, YY - CY) / 1e3

import tiflzw as _tl
_ter, _ = _tl.read_geotiff('/home/claude/analysis/gee03/Najran_terrain_30m.tif')
ELEV, SLOPE = np.nan_to_num(_ter[0], nan=1300.0), np.nan_to_num(_ter[1], nan=0.0)

def features(cls):
    b = (cls == 0)
    d_b = ndimage.distance_transform_edt(~b) * 0.03
    f5 = ndimage.uniform_filter(b.astype(float), 5)
    f21 = ndimage.uniform_filter(b.astype(float), 21)
    f61 = ndimage.uniform_filter(b.astype(float), 61)
    veg = cls == 1
    d_v = ndimage.distance_transform_edt(~veg) * 0.03
    F = np.stack([d_b, f5, f21, f61, d_v, dist_c, plans_m.astype(float), dom.astype(float),
                  (cls == 1).astype(float), (cls == 2).astype(float), (cls == 4).astype(float),
                  SLOPE, ELEV], -1)
    return F
FN = ['dist_built', 'built_5', 'built_21', 'built_61', 'dist_veg', 'dist_centre', 'in_plan',
      'in_domain', 'is_veg', 'is_sand', 'is_rock', 'slope', 'elevation']

def cls_of(y): return np.where(mask, L[y][0], 255).astype(int)

# ---------------------------------------------------------- training
Xs, ys_ = [], []
for y0, y1 in [(2010, 2015), (2015, 2020)]:
    c0, c1 = cls_of(y0), cls_of(y1)
    cand = mask & (c0 != 0) & (c0 != 3)
    F = features(c0)
    pos = np.argwhere(cand & (c1 == 0)); neg = np.argwhere(cand & (c1 != 0))
    neg = neg[rng.choice(len(neg), min(len(neg), 6 * len(pos)), replace=False)]
    for idx, lab_ in [(pos, 1), (neg, 0)]:
        Xs.append(F[idx[:, 0], idx[:, 1]]); ys_.append(np.full(len(idx), lab_))
X = np.vstack(Xs); Y = np.concatenate(ys_)
clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.08, max_leaf_nodes=48, random_state=0)
clf.fit(X, Y)

def allocate(cls, n_new, allowed, F=None):
    F = features(cls) if F is None else F
    cand = allowed & (cls != 0) & (cls != 3) & mask
    idx = np.argwhere(cand)
    p = clf.predict_proba(F[idx[:, 0], idx[:, 1]])[:, 1]
    take = idx[np.argsort(-p)[:n_new]]
    out = cls.copy(); out[take[:, 0], take[:, 1]] = 0
    return out, p

# ---------------------------------------------------------- validation 2020 -> 2026
c20, c26 = cls_of(2020), cls_of(2026)
n_obs = int(np.sum(mask & (c20 != 0) & (c26 == 0)))
sim26, _ = allocate(c20, n_obs, np.ones_like(mask))
obs_ch = mask & (c20 != 0) & (c26 == 0); sim_ch = mask & (c20 != 0) & (sim26 == 0)
hits = np.sum(obs_ch & sim_ch); misses = np.sum(obs_ch & ~sim_ch); fa = np.sum(~obs_ch & sim_ch)
fom = hits / (hits + misses + fa)
# null model: nearest to existing built-up
d20 = ndimage.distance_transform_edt(c20 != 0)
cand = mask & (c20 != 0) & (c20 != 3); idx = np.argwhere(cand)
take = idx[np.argsort(d20[cand])[:n_obs]]
null = np.zeros_like(mask); null[take[:, 0], take[:, 1]] = True
fom_null = np.sum(obs_ch & null) / (np.sum(obs_ch | null))
# relaxed (within 90 m) agreement
obs_d = ndimage.binary_dilation(obs_ch, iterations=3)
fom_relax = np.sum(sim_ch & obs_d) / sim_ch.sum()

# ---------------------------------------------------------- scenarios to 2040
POP26_LP = 380000 * (490000 / 380000) ** (4 / 18)       # Local Plan trajectory
dPOP = 490000 - POP26_LP
rate_bau = (259.2 - 187.8) / 6                              # km2/yr, 2020-2026
dem_bau = int(round(rate_bau * 14 / PIX))
DENS_LP = 33                                                # persons/ha, Local Plan core target
dem_lp = int(round(dPOP / DENS_LP / 100 / PIX))             # ha -> km2 -> pixels
farm = c26 == 1
res = {}
for name, n_new, allowed in [('BAU', dem_bau, np.ones_like(mask)),
                             ('Boundary', dem_bau, dom & ~farm),
                             ('LocalPlan', dem_lp, dom & ~farm)]:
    c = c26.copy()
    for step in range(2):                                   # 2033, 2040
        c, _ = allocate(c, n_new // 2 if step == 0 else n_new - n_new // 2, allowed)
    res[name] = c
np.savez_compressed('/home/claude/analysis/scenarios_2040.npz', BAU=res['BAU'].astype(np.uint8),
                    BND=res['Boundary'].astype(np.uint8), LP=res['LocalPlan'].astype(np.uint8), c26=c26.astype(np.uint8))

# ---------------------------------------------------------- metrics
CONTRAST = {1: 3.0, 2: -2.9, 4: 1.5}   # dLST (built - former cover), mean of epochs (Table 2)
def pairdist(b, n=6000):
    idx = np.argwhere(b); idx = idx[rng.choice(len(idx), n, replace=False)]
    P = np.c_[xs[idx[:, 1]], ys[idx[:, 0]]] / 1e3
    a, bb = P[:n // 2], P[n // 2:]
    return np.mean(np.hypot(*(a[:, None, :] - bb[None, :, :]).transpose(2, 0, 1)))
rows = []
for name, c in list(res.items()) + [('2026', c26)]:
    new = mask & (c26 != 0) & (c == 0)
    built = mask & (c == 0)
    row = dict(scenario=name, new_km2=new.sum() * PIX,
               farm_lost_km2=np.sum(new & (c26 == 1)) * PIX,
               from_sand_km2=np.sum(new & (c26 == 2)) * PIX,
               from_rock_km2=np.sum(new & (c26 == 4)) * PIX,
               outside_domain_km2=np.sum(new & ~dom) * PIX,
               mean_dist_new_km=dist_c[new].mean() if new.any() else np.nan,
               built_total_km2=built.sum() * PIX,
               density_pph=490000 / (built.sum() * 100 * PIX) if name != '2026' else POP26_LP / (built.sum() * 100 * PIX),
               mean_pair_dist_km=pairdist(built),
               mean_dist_all_built_km=dist_c[built].mean(),
               surface_heat_index=sum(np.sum(new & (c26 == k)) * PIX * v for k, v in CONTRAST.items()))
    rows.append(row)
m = pd.DataFrame(rows).set_index('scenario')
m.to_csv('/home/claude/analysis/scenario_metrics.csv')
imp = pd.Series(np.abs(clf.predict_proba(X[:1])[0]))  # placeholder to keep API simple
print(f'centre E{CX:.0f} N{CY:.0f}; train n={len(Y)} pos={Y.sum()}')
print(f'validation 2020-26: n_obs={n_obs} ({n_obs*PIX:.1f} km2) FoM={fom:.3f} null FoM={fom_null:.3f} relaxed-90m precision={fom_relax:.3f}')
print(f'POP26_LP={POP26_LP:.0f} dPOP={dPOP:.0f} dem_bau={dem_bau*PIX:.1f} km2 dem_lp={dem_lp*PIX:.1f} km2')
print(m.round(2).T)
json.dump(dict(fom=fom, fom_null=fom_null, fom_relax=fom_relax, n_obs_km2=n_obs * PIX, CX=CX, CY=CY,
               pop26=POP26_LP, dpop=dPOP, dem_bau=dem_bau * PIX, dem_lp=dem_lp * PIX),
          open('/home/claude/analysis/scenario_meta.json', 'w'), indent=1)

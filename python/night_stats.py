import sys; sys.path.insert(0, '/home/claude/tools')
import numpy as np, tifffile as T, json
from tiflzw import read_geotiff
night, gt = read_geotiff('night.tif')  # gt = (x0, dx, y0, dy)
NB, NR, NC = night.shape
R = '/home/claude/Najran_LowCarbon_ArcMap/01_Rasters/'
def agg(path, vals, nodata=255):
    tf = T.TiffFile(path); pg = tf.pages[0]; a = pg.asarray(); ny, nx = a.shape
    tp = pg.tags['ModelTiepointTag'].value; sc = pg.tags['ModelPixelScaleTag'].value
    x = tp[3] + sc[0] * (np.arange(nx) + 0.5); y = tp[4] - sc[1] * (np.arange(ny) + 0.5)
    ci = ((x - gt[0]) // 1000).astype(int); ri = ((gt[2] - y) // 1000).astype(int)
    CI, RI = np.meshgrid(ci, ri); ok = (CI >= 0) & (CI < NC) & (RI >= 0) & (RI < NR)
    idx = RI * NC + CI
    out = {}
    valid = ok & (a != nodata) if vals is not None else ok & np.isfinite(a) & (a > -1000)
    tot = np.bincount(idx[valid], minlength=NR * NC).reshape(NR, NC)
    if vals is None:
        s = np.bincount(idx[valid], weights=a[valid].astype(float), minlength=NR * NC).reshape(NR, NC)
        return s / np.maximum(tot, 1), tot
    for k, v in vals.items():
        out[k] = np.bincount(idx[valid & (a == v)], minlength=NR * NC).reshape(NR, NC) / np.maximum(tot, 1)
    return out, tot
elev, _ = agg(R + 'SRTM_elevation_m.tif', None)
EP = [2000, 2005, 2010, 2015, 2020, 2026]
F = {e: agg(R + f'LULC_{e}.tif', {'B': 0, 'G': 1, 'S': 2, 'R': 4}) for e in EP}
years = list(range(2003, 2026)); res = []
for b, yr in enumerate(years):
    e = min(EP, key=lambda z: abs(z - yr)); f, tot = F[e]
    L = night[b]; m = np.isfinite(L) & (tot > 0.8 * (1000 / 30) ** 2) & (L > 15)
    X = np.c_[np.ones(m.sum()), f['B'][m], f['G'][m], f['S'][m], elev[m] / 100]
    beta, *_ = np.linalg.lstsq(X, L[m], rcond=None)
    pred = X @ beta; r2 = 1 - ((L[m] - pred)**2).sum() / ((L[m] - L[m].mean())**2).sum()
    cls = {k: float(np.nanmean(L[m & (f[k] >= th)])) if (m & (f[k] >= th)).sum() >= 3 else None
           for k, th in [('B', 0.5), ('G', 0.4), ('S', 0.8), ('R', 0.8)]}
    n = {k: int((m & (f[k] >= th)).sum()) for k, th in [('B', 0.5), ('G', 0.4), ('S', 0.8), ('R', 0.8)]}
    res.append(dict(year=yr, epoch=e, n=int(m.sum()), coef_built=beta[1], coef_veg=beta[2], coef_sand=beta[3],
                    coef_elev100=beta[4], r2=r2, cls=cls, ncls=n))
for r in res:
    print(r['year'], r['epoch'], r['n'], 'B %.2f G %.2f S %.2f elev %.2f R2 %.2f' % (r['coef_built'], r['coef_veg'], r['coef_sand'], r['coef_elev100'], r['r2']),
          {k: (round(v, 2) if v else None) for k, v in r['cls'].items()}, r['ncls'])
json.dump(res, open('night_stats.json', 'w'), indent=1, default=float)

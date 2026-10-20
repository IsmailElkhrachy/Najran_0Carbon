"""Blind stratified accuracy for the 2000 map (Olofsson et al. 2014).
Input: Najran_VAL2000_label.csv returned by the author (ref_class in B,V,G,S,R)."""
import numpy as np, pandas as pd, json, sys
key = pd.read_csv('/home/claude/Najran_labeling_kit/Najran_VAL_key_2000.csv')[['pid', 'map_class']]
lab = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else '/home/claude/Najran_VAL2000_label.csv', encoding='utf-8-sig')
d = key.merge(lab[['pid', 'ref_class']], on='pid'); d['ref'] = d.ref_class.str.strip().str.upper()
d['map'] = d.map_class.map({0: 'B', 1: 'G', 2: 'S', 4: 'R'})
assert d.ref.isin(list('BVGSR')).all(), d[~d.ref.isin(list('BVGSR'))]
area = {'B': 85.93, 'G': 120.77, 'S': 456.85, 'R': 1780.67}
M = list('BGSR'); A = sum(area.values()); W = np.array([area[m] / A for m in M])
n = np.array([(d['map'] == m).sum() for m in M])
out = {}
for label, f in [('surface_def', lambda x: 'S' if x == 'V' else x), ('urbanised_def', lambda x: 'B' if x == 'V' else x)]:
    r = d.ref.map(f); N = np.array([[((d['map'] == m) & (r == c)).sum() for c in M] for m in M])
    fr = N / n[:, None]; p = W[:, None] * fr
    oa = np.trace(p); se = np.sqrt(sum(W[i]**2 * fr[i, i] * (1 - fr[i, i]) / (n[i] - 1) for i in range(4)))
    ua = np.diag(fr); pa = np.diag(p) / p.sum(0)
    adj = p.sum(0) * A; se_a = np.sqrt(((W[:, None]**2) * fr * (1 - fr) / (n[:, None] - 1)).sum(0)) * A
    out[label] = dict(N=N.tolist(), OA=round(100 * oa, 1), OA_ci=round(196 * se, 1),
                      UA=dict(zip(M, np.round(100 * ua, 1).tolist())), PA=dict(zip(M, np.round(100 * pa, 1).tolist())),
                      area_adj={m: (round(a, 1), round(1.96 * s, 1)) for m, a, s in zip(M, adj, se_a)})
out['n_V'] = int((d.ref == 'V').sum())
print(json.dumps(out, indent=1)); json.dump(out, open('/home/claude/analysis/stats/accuracy_2000.json', 'w'), indent=1)

"""Workflow flowchart for the Najran low-carbon paper (vector PDF + PNG)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

plt.rcParams.update({'font.family': 'DejaVu Sans'})
fig, ax = plt.subplots(figsize=(16, 10.5))
ax.set_xlim(0, 160); ax.set_ylim(0, 105); ax.axis('off')

C = dict(data='#e8eef6', proc='#f3f3f3', rq1='#fde0dd', rq2='#e5f5e0', rq3='#fff2cc',
         rq4='#deebf7', out='#efe3f7')
EC = dict(rq1='#cb181d', rq2='#238b45', rq3='#d79b00', rq4='#2171b5', out='#6a3d9a')
B = {}

def box(key, x, y, w, h, text, fc, ec='#555555', fs=11.8, bold_first=True, lw=0.9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.25,rounding_size=1.2',
                                fc=fc, ec=ec, lw=lw))
    lines = text.split('\n')
    if bold_first:
        ax.text(x + w / 2, y + h - 1.6, lines[0], ha='center', va='top', fontsize=fs, weight='bold')
        if len(lines) > 1:
            ax.text(x + w / 2, y + h - 4.6, '\n'.join(lines[1:]), ha='center', va='top',
                    fontsize=fs - 1.6, linespacing=1.18)
    else:
        ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs)
    B[key] = (x, y, w, h)

def arrow(a, b, side_a='r', side_b='l', color='#666666', lw=1.0, ya=0.5, yb=0.5):
    xa, ya_, wa, ha = B[a]; xb, yb_, wb, hb = B[b]
    pa = (xa + wa, ya_ + ha * ya) if side_a == 'r' else (xa + wa / 2, ya_)
    pb = (xb, yb_ + hb * yb) if side_b == 'l' else (xb + wb / 2, yb_ + hb)
    ax.annotate('', xy=pb, xytext=pa, arrowprops=dict(arrowstyle='-|>', color=color, lw=lw,
                shrinkA=2, shrinkB=2, connectionstyle='arc3,rad=0'))

# column headers
for x, w, t in [(1, 33, '1  Data'), (40, 43, '2  Processing & models'),
                (89, 33, '3  Results (RQs)'), (128, 31, '4  Planning outputs')]:
    ax.text(x + w / 2, 102.5, t, ha='center', va='center', fontsize=15, weight='bold')
    ax.plot([x, x + w], [100.3, 100.3], color='#333333', lw=1.2)

# ---------------- 1 data
D = [('landsat', 'Landsat 5/7/8/9 (C2 L2)\nreflectance + LST, 30 m\n2000–2026, six epochs'),
     ('aux', 'Ancillary layers\nSRTM DEM · WorldCover\nGHSL built-up'),
     ('vhr', 'Google Earth VHR\n240 blind reference\npoints'),
     ('hydro', 'Water records\nCHIRPS · GRACE/-FO\nGLDAS-2.2 · station'),
     ('muni', 'Municipal & planning\nplans · 18,096 parcels\nLocal Plan · 63 zones'),
     ('bldg', 'Building footprints\n104,952 polygons\n(Microsoft, via VIDA)'),
     ('socio', 'Population & energy\nCensus · UN WPP\nGASTAT electricity'),
     ('param', 'Parameters & evidence\nPV yield · grid EF\nelasticities · testimony')]
y = 90.0
for k, t in D:
    box(k, 1, y, 33, 9.2, t, C['data'], ec='#4a6fa5')
    y -= 9.2 + 2.0

# ---------------- 2 processing
P = [('pre', 'Pre-processing (GEE)\ncloud mask · sensor harmonisation\nepoch composites + indices', 9.8),
     ('rf', 'Random Forest classification\n5 classes incl. rock/gravel\nfiltered labels · no-reversion rule', 9.8),
     ('acc', 'Accuracy & area estimation\nOlofsson estimators · 95% CI\nvacant subdivided land', 9.8),
     ('sdg', 'Land-use efficiency\nSDG 11.3.1 · m² per resident\nplans and parcels', 9.8),
     ('ndvi', 'Vegetation series & transitions\nMann–Kendall · Sen · Pettitt\nfate of 2000 farmland', 9.8),
     ('drv', 'Driver attribution\nrain · water storage\npolicy & market timeline', 8.6),
     ('lst', 'Surface heat\nLST · SUHII · relative LST\n63-zone statistics', 8.6),
     ('mlca', 'ML–CA scenarios to 2040\ngradient boosting · FoM\nTrend · Boundary · Plan', 9.8),
     ('mc', 'Carbon (Monte Carlo)\nrooftop PV\ncompact-growth travel', 8.6)]
y = 90.0
for k, t, h in P:
    box(k, 40, y, 43, 9.2, t, C['proc'])
    y -= 9.2 + 1.2

# ---------------- 3 results
R = [('rq1', 'RQ1 Urban expansion\nbuilt-up 86→259 km²\nLCRPGR 0.34→2.97\n288 km² vacant plots', 'rq1', 76.5),
     ('rq2', 'RQ2 Oasis decline\nvegetation −39%\nbreak after 2008\ngroundwater-driven', 'rq2', 57.5),
     ('rq3', 'RQ3 Surface heat\nbuilt < sand, > rock\nvegetation −3 °C\nfarm loss +1.3 °C', 'rq3', 38.5),
     ('rq4', 'RQ4 Priorities & carbon\nTrend 167 vs Plan 27 km²\nrooftop PV ≈1.8 GWp\n1.5 Mt CO₂/yr avoided', 'rq4', 16.0)]
for k, t, c, y in R:
    box(k, 89, y, 33, 15, t, C[c], ec=EC[c], lw=1.6, fs=12.2)

# ---------------- 4 outputs
O = [('o1', 'Plan verification\nNajran Local Plan 2040\nboundary · infill\ndensity', 70),
     ('o2', 'Citizen atlas\n63 zones\nheat & greening\nfarmland · housing\nrooftop solar', 44),
     ('o3', 'Low-carbon pathway\nprotect oasis\nrecharge aquifer\ninfill · rooftop PV', 18)]
for k, t, y in O:
    box(k, 128, y, 31, 18, t, C['out'], ec=EC['out'], lw=1.6, fs=12.2)

# ---------------- arrows: data -> processing
g = '#7a8fb0'
arrow('landsat', 'pre', color=g); arrow('aux', 'rf', color=g); arrow('vhr', 'acc', color=g)
arrow('hydro', 'drv', color=g); arrow('muni', 'sdg', color=g); arrow('muni', 'lst', color=g, ya=0.3)
arrow('bldg', 'mc', color=g); arrow('socio', 'sdg', color=g, ya=0.8, yb=0.3); arrow('socio', 'mc', color=g, ya=0.3, yb=0.6)
arrow('param', 'mc', color=g, yb=0.3); arrow('muni', 'mlca', color=g, ya=0.15, yb=0.7)
# processing chain (vertical)
for a, b in [('pre', 'rf'), ('rf', 'acc'), ('acc', 'sdg'), ('ndvi', 'drv'), ('drv', 'lst')]:
    arrow(a, b, side_a='b', side_b='t', color='#444444')
# processing -> results
arrow('acc', 'rq1', color=EC['rq1'], lw=1.4); arrow('sdg', 'rq1', color=EC['rq1'], lw=1.4, yb=0.3)
arrow('ndvi', 'rq2', color=EC['rq2'], lw=1.4, yb=0.7); arrow('drv', 'rq2', color=EC['rq2'], lw=1.4, yb=0.3)
arrow('lst', 'rq3', color=EC['rq3'], lw=1.4)
arrow('lst', 'rq4', color=EC['rq4'], lw=1.4, ya=0.2, yb=0.85)
arrow('mlca', 'rq4', color=EC['rq4'], lw=1.4, yb=0.6); arrow('mc', 'rq4', color=EC['rq4'], lw=1.4, yb=0.3)
# results -> outputs
arrow('rq1', 'o1', color=EC['out'], lw=1.4); arrow('rq2', 'o1', color=EC['out'], lw=1.4, yb=0.25)
arrow('rq2', 'o2', color=EC['out'], lw=1.4, ya=0.3, yb=0.8); arrow('rq3', 'o2', color=EC['out'], lw=1.4, yb=0.45)
arrow('rq4', 'o2', color=EC['out'], lw=1.4, ya=0.8, yb=0.15); arrow('rq4', 'o3', color=EC['out'], lw=1.4)


for ext in ['pdf', 'png']:
    fig.savefig(f'/home/claude/najran_paper/figures/fig_workflow.{ext}', dpi=300, bbox_inches='tight')

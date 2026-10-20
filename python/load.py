import sys, numpy as np
sys.path.insert(0,'/home/claude/tools'); import tiflzw, shpread as s
from matplotlib.path import Path
D='/tmp/claude-0/-home-claude/40128d67-1783-5f2f-811e-6fee803b10b9/scratchpad/v3/Najran_LowCarbon_GEE_v2/'
SHP='/tmp/claude-0/-home-claude/40128d67-1783-5f2f-811e-6fee803b10b9/scratchpad/shp/'
YEARS=[2000,2005,2010,2015,2020,2026]
def load_all():
    L={};T={}
    for y in YEARS:
        a,gt=tiflzw.read_geotiff(D+f'Najran_LULC_{y}.tif'); L[y]=a
        b,_=tiflzw.read_geotiff(D+f'Najran_LST_{y}.tif'); T[y]=b[0]
    tr,_=tiflzw.read_geotiff(D+'Najran_NDVI_trend.tif')
    mask=~np.isnan(T[2026])
    return L,T,tr,mask,gt
def pixel_xy(gt,shape):
    h,w=shape; xs=gt[0]+gt[1]*(np.arange(w)+0.5); ys=gt[2]+gt[3]*(np.arange(h)+0.5); return xs,ys
def rasterize(parts,gt,shape):
    xs,ys=pixel_xy(gt,shape)
    allp=np.vstack(parts); x0,y0=allp.min(0); x1,y1=allp.max(0)
    ci=np.where((xs>=x0)&(xs<=x1))[0]; ri=np.where((ys>=y0)&(ys<=y1))[0]
    m=np.zeros(shape,bool)
    if len(ci)==0 or len(ri)==0: return m
    X,Y=np.meshgrid(xs[ci],ys[ri]); pts=np.c_[X.ravel(),Y.ravel()]
    inside=np.zeros(len(pts),bool)
    for r in parts: inside^=Path(r).contains_points(pts)
    m[np.ix_(ri,ci)]=inside.reshape(len(ri),len(ci)); return m

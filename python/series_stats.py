import numpy as np, pandas as pd, json
from scipy import stats
v=pd.read_csv('veg_annual.csv')[['year','veg_km2','n_months_with_data']].sort_values('year')
rain=pd.read_csv('rain_annual.csv')[['year','rain_mm']]
gr=pd.read_csv('../grace_clean.csv'); gr['year']=gr.date.str[:4].astype(int)
gl=pd.read_csv('../gldas_clean.csv'); gl['year']=gl.date.str[:4].astype(int)
def mk(x):
    x=np.asarray(x,float); n=len(x); s=0
    for i in range(n-1): s+=np.sign(x[i+1:]-x[i]).sum()
    # ties correction
    _,c=np.unique(x,return_counts=True)
    var=(n*(n-1)*(2*n+5)-sum(t*(t-1)*(2*t+5) for t in c))/18
    z=(s-np.sign(s))/np.sqrt(var) if s!=0 else 0
    p=2*(1-stats.norm.cdf(abs(z)))
    tau=stats.kendalltau(np.arange(n),x)[0]
    slopes=[(x[j]-x[i])/(j-i) for i in range(n) for j in range(i+1,n)]
    return dict(n=n,S=int(s),Z=round(z,2),p=round(p,4),tau=round(tau,2),sen=round(np.median(slopes),2))
def pettitt(x):
    x=np.asarray(x,float); n=len(x)
    U=[sum(np.sign(x[i]-x[j]) for i in range(t+1) for j in range(t+1,n)) for t in range(n-1)]
    t=int(np.argmax(np.abs(U))); K=abs(U[t]); p=2*np.exp(-6*K**2/(n**3+n**2))
    return t,K,min(p,1)
y=v.veg_km2.values; yrs=v.year.values
res={'full':mk(y)}
t,K,p=pettitt(y); res['pettitt']=dict(break_after=int(yrs[t]),K=int(K),p=round(p,4))
res['pre']=mk(y[yrs<=yrs[t]]); res['post']=mk(y[yrs>yrs[t]])
res['pre_mean']=round(y[yrs<=yrs[t]].mean(),1); res['post_mean']=round(y[yrs>yrs[t]].mean(),1)
res['welch_pre_post']=[round(v,4) for v in stats.ttest_ind(y[yrs<=yrs[t]],y[yrs>yrs[t]],equal_var=False)]
res['rain_mk']=mk(rain.dropna().rain_mm.values); res['rain_mean']=round(rain.rain_mm.mean(),1)
# annual driver series
gra=gr.groupby('year').lwe_cm.mean(); glm=gl.groupby('year').gws_mm.min(); ra=rain.set_index('year').rain_mm
def corr(a,b,name):
    j=pd.concat([a,b],axis=1).dropna(); r,p=stats.pearsonr(j.iloc[:,0],j.iloc[:,1]); rs,ps=stats.spearmanr(j.iloc[:,0],j.iloc[:,1])
    res[name]=dict(n=len(j),r=round(r,2),p=round(p,4),rho=round(rs,2),p_rho=round(ps,4),years=f'{j.index.min()}-{j.index.max()}')
vs=v.set_index('year').veg_km2
corr(vs,ra,'veg_rain'); corr(vs,gra,'veg_grace'); corr(vs,glm,'veg_gldasmin')
corr(vs[vs.index<=2012],gra,'veg_grace_2002_2012')
# GRACE trend (monthly OLS)
g=gr.dropna(); sl=stats.linregress(g.yr,g.lwe_cm); res['grace_trend']=dict(slope_cm_yr=round(sl.slope,3),se=round(sl.stderr,3),p=sl.pvalue,n=len(g))
gm=glm[(glm.index>=2003)&(glm.index<=2025)]
for a,b in [(2003,2008),(2009,2025)]:
    s=gm[(gm.index>=a)&(gm.index<=b)]; l=stats.linregress(s.index,s.values); res[f'gldas_{a}_{b}']=dict(slope=round(l.slope,2),se=round(l.stderr,2),p=round(l.pvalue,4))
print(json.dumps(res,indent=1,default=float)); json.dump(res,open('series_stats.json','w'),indent=1,default=float)

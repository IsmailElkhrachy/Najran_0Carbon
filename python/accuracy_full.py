"""Stratified area-adjusted accuracy (Olofsson et al. 2014) for the 2026 map with a
separate reference stratum for vacant subdivided land (V)."""
import numpy as np, pandas as pd, json
A='/tmp/claude-0/-home-claude/40128d67-1783-5f2f-811e-6fee803b10b9/scratchpad/acc/'
k=pd.read_csv(A+'Najran_VAL_key_2026.csv'); l=pd.read_csv(A+'Najran_VAL_label_2026.csv'); r=pd.read_csv(A+'Najran_VAL_review_2026.csv')
d=k[['pid','map_class']].merge(l[['pid','ref_class']],on='pid')
vac=set(r.loc[r.new_class==0,'pid'])
d['ref6']=d.ref_class.astype(int).astype(str).map({'0':'B','1':'G','2':'S','4':'R'})
d.loc[d.pid.isin(vac),'ref6']='V'
# Al-Ghuwaila agricultural plan: agricultural title, no residential deeds, building
# permits suspended by the municipality (resident enquiry, Oct 2026) -> not
# committed urban land; reference class returned to sand.
GHUWAILA=['2026_156','2026_162','2026_131','2026_129']
d.loc[d.pid.isin(GHUWAILA),'ref6']='S'
d['map4']=d.map_class.map({0:'B',1:'G',2:'S',4:'R'})
area_map={'B':259.25,'G':68.79,'S':450.05,'R':1666.12}
M=['B','G','S','R']; Atot=sum(area_map.values()); W=np.array([area_map[m]/Atot for m in M])
n_i=np.array([(d.map4==m).sum() for m in M])
def matrix(refcols, refmap):
    N=np.array([[((d.map4==m)&(d.ref6.map(refmap)==c)).sum() for c in refcols] for m in M])
    return N
def olof(N, cols):
    f=N/n_i[:,None]; p=W[:,None]*f
    pj=p.sum(0); A_adj=pj*Atot
    se=np.sqrt(((W[:,None]**2)*f*(1-f)/(n_i[:,None]-1)).sum(0))*Atot
    return p,pj,A_adj,se
out={}
# scheme 6-ref columns for area estimation
cols6=['B','V','G','S','R']; N6=matrix(cols6,lambda x:x)
p6,pj6,A6,se6=olof(N6,cols6)
out['N6']=N6.tolist(); out['area6']={c:(round(a,1),round(1.96*s,1)) for c,a,s in zip(cols6,A6,se6)}
# urbanised = B+V
f=N6/n_i[:,None]; fu=f[:,0]+f[:,1]
se_u=np.sqrt(((W**2)*fu*(1-fu)/(n_i-1)).sum())*Atot
out['urbanised']=(round((A6[0]+A6[1]),1), round(1.96*se_u,1))
def acc(refmap, label):
    N=matrix(M, refmap); f=N/n_i[:,None]; p=W[:,None]*f
    oa=np.trace(p); se_oa=np.sqrt(sum(W[i]**2*f[i,i]*(1-f[i,i])/(n_i[i]-1) for i in range(4)))
    ua=np.diag(f); se_ua=np.sqrt(ua*(1-ua)/(n_i-1))
    pa=np.diag(p)/p.sum(0)
    # PA SE (Olofsson 2014 eq. 7) via areas
    Nj=np.array([sum(area_map[M[i]]/n_i[i]*N[i,j] for i in range(4)) for j in range(4)])
    se_pa=[]
    for j in range(4):
        t1=area_map[M[j]]**2*(1-pa[j])**2*ua[j]*(1-ua[j])/(n_i[j]-1)
        t2=pa[j]**2*sum(area_map[M[i]]**2*f[i,j]*(1-f[i,j])/(n_i[i]-1) for i in range(4) if i!=j)
        se_pa.append(np.sqrt(t1+t2)/Nj[j])
    # 3-class (bare = S+R)
    g3={'B':0,'G':1,'S':2,'R':2}
    p3=np.zeros((3,3))
    for i in range(4):
        for j in range(4): p3[g3[M[i]],g3[M[j]]]+=p[i,j]
    out[label]=dict(N=N.tolist(),OA=round(100*oa,1),OA_ci=round(196*se_oa,1),UA=dict(zip(M,np.round(100*ua,1).tolist())),
        UA_ci=dict(zip(M,np.round(196*se_ua,1).tolist())),PA=dict(zip(M,np.round(100*pa,1).tolist())),
        PA_ci=dict(zip(M,np.round(196*np.array(se_pa),1).tolist())),OA3=round(100*np.trace(p3),1),p=np.round(p,4).tolist())
acc(lambda x: 'B' if x in('B','V') else x, 'urbanised_def')     # V counted as built
acc(lambda x: 'S' if x=='V' else x, 'surface_def')               # V counted as bare
print(json.dumps(out,indent=1)); json.dump(out,open('accuracy_full.json','w'),indent=1)

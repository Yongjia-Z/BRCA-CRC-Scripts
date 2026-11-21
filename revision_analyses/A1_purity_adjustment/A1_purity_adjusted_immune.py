import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
ROOT=Path(sys.argv[1]); OUT=Path(sys.argv[2])
d=pd.read_csv(Path(sys.argv[3]) if len(sys.argv)>3 else ROOT/'core_data/remap_reanalysis/figure5-8/figure6/figure6_crc13_scores_with_tumor_state_mcp.tsv',sep='\t')
d=d.drop(columns=[c for c in ['context'] if c in d.columns])
pur=pd.read_excel(OUT/'Aran2015_SupData1.xlsx',skiprows=3)
pur=pur[pur['Cancer type'].eq('BRCA')][['Sample ID','CPE','ESTIMATE','ABSOLUTE','LUMP','IHC']].rename(columns={'Sample ID':'sample_id'})
d['sid16']=d.sample_id.str[:16]
d=d.merge(pur.drop_duplicates('sample_id'),left_on='sid16',right_on='sample_id',how='left',suffixes=('','_pur'))
ctx=d.SUBTYPE.map({'BRCA_LumA':'Luminal','BRCA_LumB':'Luminal','BRCA_Her2':'HER2','BRCA_Basal':'Basal/TNBC','BRCA_Normal':'Normal-like'}).fillna('Unknown')
d['context']=ctx
scores=['pan_CRC13_score','luminal_CRC13_score','her2_CRC13_score','basal_TNBC_CRC13_score']
mcp=['T cells','CD8 T cells','Cytotoxic lymphocytes','NK cells','B lineage','Monocytic lineage','Myeloid dendritic cells','Neutrophils','Endothelial cells','Fibroblasts']
def resid(y,X):
    X=np.column_stack([np.ones(len(y)),X]); b=np.linalg.lstsq(X,y,rcond=None)[0]; return y-X@b
dummies=pd.get_dummies(d.context,drop_first=True).astype(float).to_numpy()
models={'unadjusted':None,'context':dummies}
ok=d.CPE.notna().to_numpy()
rows=[]
for label,cov in [('unadjusted',None),('context_adjusted',dummies),('purity_adjusted',d[['CPE']].to_numpy()),('context_plus_purity_adjusted',None)]:
    sub=ok if label in ('purity_adjusted','context_plus_purity_adjusted') else np.ones(len(d),bool)
    for s in scores:
        for m in mcp:
            y=d.loc[sub,s].to_numpy(float); x=d.loc[sub,m].to_numpy(float)
            if label=='context_adjusted': y,x=resid(y,dummies[sub]),resid(x,dummies[sub])
            elif label=='purity_adjusted': C=d.loc[sub,['CPE']].to_numpy(float); y,x=resid(y,C),resid(x,C)
            elif label=='context_plus_purity_adjusted': C=np.column_stack([dummies[sub],d.loc[sub,'CPE'].to_numpy(float)]); y,x=resid(y,C),resid(x,C)
            r,p=spearmanr(y,x); rows.append({'model':label,'CRC_score':s,'MCPcounter_population':m,'n':int(sub.sum()),'spearman_rho':r,'p':p})
res=pd.DataFrame(rows)
res['BH_q_within_model']=res.groupby('model').p.transform(lambda p: multipletests(p,method='fdr_bh')[1])
res.to_csv(OUT/'A1_purity_adjusted_MCPcounter_correlations.tsv',sep='\t',index=False)
pp=[]
for v in scores+mcp:
    r,p=spearmanr(d.loc[ok,v],d.loc[ok,'CPE']); pp.append({'variable':v,'n':int(ok.sum()),'rho_with_CPE':r,'p':p})
pp=pd.DataFrame(pp); pp.to_csv(OUT/'A1_scores_vs_purity.tsv',sep='\t',index=False)
pd.DataFrame({'context':d.context,'has_CPE':ok}).value_counts().reset_index().to_csv(OUT/'A1_sample_counts.tsv',sep='\t',index=False)
w=res.pivot_table(index=['CRC_score','MCPcounter_population'],columns='model',values='spearman_rho')[['unadjusted','context_adjusted','purity_adjusted','context_plus_purity_adjusted']]
print('n with CPE',ok.sum(),'of',len(d)); print(pp.round(3).to_string(index=False)); print(w.round(3).to_string())

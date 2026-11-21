import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, kruskal
from statsmodels.stats.multitest import multipletests
from lifelines import CoxPHFitter
ROOT=Path(sys.argv[1]); OUT=Path(sys.argv[2]); OUT.mkdir(parents=True,exist_ok=True)
A=ROOT/'source_tables'; C=ROOT/'core_data'; R=C/'remap_reanalysis'
REV=ROOT/'revision_tables'
bh=lambda p: multipletests(p,method='fdr_bh')[1]
log={}
f=pd.read_csv(A/'figure1/figure1a_all_tf_scores.csv')
cols=['RNA','TF_activity_regulon','ATAC_motif','ATAC_footprint','ATAC_peak_to_gene_accessibility','context_matched_H3K27ac_SE_MAX','context_matched_direct_binding_MAX']
core=f[cols].fillna(0).sum(1); s=f.max_state_score_0_5; sp=f.state_spread_max_minus_second; sup=f.Supportive_score_0_3
high=(core>=5.5)&(s>=2.6)
prov=((core>=4.75)&(s>=3)&(sup>=.5))|((core>=4.75)&(s>=3.5)&(sp>=1))|((core+sup>=6)&(s>=2.4))
sel=(high|prov)&(s>=2.35)
log['n_TF_universe']=len(f); log['formal']=f.TF[sel].tolist(); log['high_conf']=f.TF[sel&high].tolist()
log['core_score_matches_stored']=bool(np.allclose(core,f.Core_discovery_score_v2_0_7,atol=1e-6))
rs=pd.read_csv(REV/'R1_5_selection_stage_random_sets.tsv',sep='\t')
log['random_sets']={'n':len(rs),'median_n_passing':float(rs.n_passing_original_rules.median()),'max_n_passing':int(rs.n_passing_original_rules.max()),'mean_core_random_median':float(rs.mean_core_score.median()),'observed_mean_core_formal':float(core[sel].mean()),'frac_random_mean_core_ge_observed':float((rs.mean_core_score>=core[sel].mean()).mean())}
w=pd.read_csv(REV/'R1_5_weight_sensitivity_all_TFs.tsv',sep='\t')
log['weight_freq_formal']=w.set_index('TF').loc[log['formal'],'selection_frequency_2000_weight_perturbations'].round(3).to_dict()
log['weight_freq_nonformal_ge_0.5']=w[~w.TF.isin(log['formal'])&(w.selection_frequency_2000_weight_perturbations>=.5)][['TF','selection_frequency_2000_weight_perturbations']].values.tolist()
pk=pd.read_csv(C/'figure3/raw/figure3d__all_local_peak_RNA_correlations.tsv',sep='\t')
pk['q_all']=bh(pk.peak_RNA_p)
d3=pd.read_csv(A/'figure3/figure3d.csv')
m=d3.merge(pk[['TF','Peak_Name','q_all']],left_on=['TF','best_local_peak'],right_on=['TF','Peak_Name'],how='left')
m['p_bonf_locus']=(m.best_local_peak_RNA_p*m.n_local_peaks).clip(upper=1); m['q_bonf_then_BH']=bh(m.p_bonf_locus)
log['fig3d']={'n_peak_tests':len(pk),'n_TF_tested':int(pk.TF.nunique()),'n_selected':len(m),'pass_allpeak_BH':int((m.q_all<.05).sum()),'pass_bonf_then_BH':int((m.q_bonf_then_BH<.05).sum())}
m.to_csv(OUT/'fig3d_peak_selection_recomputed.tsv',sep='\t',index=False)
cox=pd.read_csv(A/'figure8b/figure8b_tcga_os_pfs_cox_ph_truncated_results.tsv',sep='\t'); cox['q_BH12']=bh(cox.pvalue)
mi=pd.read_csv(A/'figure8b/figure8b_tcga_cox_model_input.tsv',sep='\t')
rec=[]
for (st,ep),g in mi.groupby(['state_id','endpoint']):
    x=g[['time_truncated','event_truncated','score_z']].dropna()
    fit=CoxPHFitter().fit(x,'time_truncated','event_truncated'); r=fit.summary.loc['score_z']
    rec.append({'state_id':st,'endpoint':ep,'n_re':len(x),'events_re':int(x.event_truncated.sum()),'hr_recomputed':r['exp(coef)'],'p_recomputed':r.p})
rec=pd.DataFrame(rec); cx=cox.merge(rec,on=['state_id','endpoint'])
log['cox_unadj_max_abs_HR_diff']=float((cx.hr_per_1sd-cx.hr_recomputed).abs().max())
clin=pd.read_csv(ROOT/'datasets/TCGA/results/TCGA_CRC_clinical_analysis_dataset.tsv',sep='\t').drop_duplicates('sample_id')
adj=mi.merge(clin[['sample_id','AGE','AJCC_PATHOLOGIC_TUMOR_STAGE_clean','SUBTYPE']],on='sample_id',how='left')
st_=adj.AJCC_PATHOLOGIC_TUMOR_STAGE_clean.astype(str)
adj['late_stage']=st_.map(lambda v:0 if v in ['I','IA','IB','II','IIA','IIB'] else (1 if v in ['III','IIIA','IIIB','IIIC','IV'] else np.nan))
adj['age_10y']=pd.to_numeric(adj.AGE,errors='coerce')/10; adj['LumB']=adj.SUBTYPE.eq('BRCA_LumB').astype(int)
ar=[]
for (st,ep),g in adj.groupby(['state_id','endpoint']):
    cv=['score_z','age_10y','late_stage']+(['LumB'] if st.startswith('Luminal') else [])
    x=g[['time_truncated','event_truncated',*cv]].dropna(); fit=CoxPHFitter().fit(x,'time_truncated','event_truncated'); r=fit.summary.loc['score_z']
    ar.append({'state_id':st,'endpoint':ep,'n_adj':len(x),'events_adj':int(x.event_truncated.sum()),'hr_adj':r['exp(coef)'],'lo_adj':r['exp(coef) lower 95%'],'hi_adj':r['exp(coef) upper 95%'],'p_adj':r.p})
ar=pd.DataFrame(ar); ar['q_adj_BH12']=bh(ar.p_adj)
cxa=cx.merge(ar,on=['state_id','endpoint'])
cxa.to_csv(OUT/'cox_original_BH_and_adjusted_recomputed.tsv',sep='\t',index=False)
log['cox']=cxa[['state_id','endpoint','n','events','n_re','hr_per_1sd','pvalue','q_BH12','n_adj','hr_adj','lo_adj','hi_adj','p_adj','q_adj_BH12']].round(4).to_dict('records')
ez=pd.read_csv(REV/'R1_6_DepMap_expression_z_51_models.tsv',sep='\t').set_index('ModelID')
ge=pd.read_csv(REV/'R1_6_DepMap_gene_effect_51_models.tsv',sep='\t').set_index('ModelID')
nul=pd.read_csv(REV/'R1_6_DepMap_empirical_null_replicates.tsv',sep='\t')
ct=pd.read_csv(REV/'R1_6_DepMap_permutation_controls.tsv',sep='\t')
dd=[]
for _,r in ct.iterrows():
    gs=r.genes.split(';'); rho=spearmanr(ez.loc[:,gs].mean(1),ge.loc[ez.index,gs].mean(1)).statistic
    nr=nul[nul.state_id==r.state_id].random_same_set_rho.to_numpy(); p=(1+np.sum(nr<=rho))/(len(nr)+1)
    dd.append({'state_id':r.state_id,'n_genes':len(gs),'rho_recomputed':rho,'rho_stored':r.reconstructed_rho_51_line_zscore,'null_median':np.median(nr),'emp_p_recomputed':p,'emp_p_stored':r.empirical_p_one_sided_matched_gene_sets,'p_labelperm_stored':r.empirical_p_two_sided_label_permutation})
dd=pd.DataFrame(dd); dd['q_matched']=bh(dd.emp_p_recomputed); dd['q_labelperm']=bh(dd.p_labelperm_stored)
dd.to_csv(OUT/'depmap_matched_null_recomputed.tsv',sep='\t',index=False); log['depmap']=dd.round(4).to_dict('records')
for nm,rel,grp,val in [('6A','figure6/figure6a.csv',['score','SUBTYPE'],'value'),('8A','figure8/figure8a.csv',['state_id','pam50'],'state_score')]:
    d=pd.read_csv(A/rel); t=d.groupby(grp)[val].agg(['count','median']).reset_index(); t.to_csv(OUT/f'fig{nm}_group_medians.tsv',sep='\t',index=False)
    top=t.loc[t.groupby(grp[0])['median'].idxmax()]; log[f'fig{nm}_top_group']=dict(zip(top[grp[0]],top[grp[1]]))
    log[f'fig{nm}_n']=d.drop_duplicates('sample_id').groupby(grp[1]).size().to_dict()
from collections import defaultdict
ext=pd.read_csv(REV/'R1_7_external_state_scores.tsv',sep='\t')
tt=[]
for (co,stt),g in ext.groupby(['cohort','state_id']):
    pr=g[g.group.isin(['Luminal','HER2-enriched','Basal-like','Normal-like'])]
    H,p=kruskal(*[x.state_score.dropna() for _,x in pr.groupby('group')]); med=pr.groupby('group').state_score.median()
    tt.append({'cohort':co,'state_id':stt,'n':len(pr),'H':H,'p':p,'top_group':med.idxmax(),**{f'median_{k}':v for k,v in med.items()}})
tt=pd.DataFrame(tt); tt.to_csv(OUT/'external_projection_recomputed.tsv',sep='\t',index=False); log['external']=tt[['cohort','state_id','n','top_group']].values.tolist()
log['external_n']={f'{a}|{b}':int(v) for (a,b),v in ext.drop_duplicates(['cohort','sample_id']).groupby(['cohort','group']).size().items()}
occ=pd.read_csv(REV/'R1_1_occupancy_14.tsv',sep='\t'); log['self_occupancy_nonzero']=occ[occ.self_occupancy_value>0].regulator_TF.tolist()
json.dump(log,open(OUT/'A0_verify_summary.json','w'),ensure_ascii=False,indent=1,default=str)
print(json.dumps(log,ensure_ascii=False,indent=1,default=str))

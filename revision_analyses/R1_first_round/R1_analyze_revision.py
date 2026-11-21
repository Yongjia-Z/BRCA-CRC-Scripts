from pathlib import Path
import sys, json, hashlib, warnings
W=Path(__file__).resolve().parent
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

ROOT=Path(sys.argv[1]) if len(sys.argv)>1 else W.parents[2]; A=ROOT/'source_tables'
C=ROOT/'core_data'; R=C/'remap_reanalysis'
OUT=Path(sys.argv[2]) if len(sys.argv)>2 else ROOT/'revision_tables'; OUT.mkdir(parents=True,exist_ok=True)
SEED=20261026; rng=np.random.default_rng(SEED); inputs=[]; summary={'seed':SEED}
def read(p,**kw):
    p=Path(p);inputs.append(p)
    return pd.read_csv(p,sep='\t' if p.suffix=='.tsv' else ',',**kw)
def save(d,n):d.to_csv(OUT/n,sep='\t',index=False,encoding='utf-8')
def bh(p):
    a=np.asarray(p,float);q=np.full(len(a),np.nan);ix=np.flatnonzero(np.isfinite(a));o=ix[np.argsort(a[ix])]
    q[o]=np.minimum.accumulate((a[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1].clip(0,1)
    return q
def report(k,v):summary[k]=v;print(k,json.dumps(v,ensure_ascii=False),flush=True)

f=read(A/'figure1/figure1a_all_tf_scores.csv').fillna({'active_contexts':''})
cols=['RNA','TF_activity_regulon','ATAC_motif','ATAC_footprint','ATAC_peak_to_gene_accessibility','context_matched_H3K27ac_SE_MAX','context_matched_direct_binding_MAX']
mat=f[cols].fillna(0).to_numpy(float);s=f['max_state_score_0_5'].to_numpy();spread=f['state_spread_max_minus_second'].to_numpy();sup=f['Supportive_score_0_3'].to_numpy()
def select(core,offset=0):
    high=(core>=5.5+offset)&(s>=2.6)
    provisional=((core>=4.75+offset)&(s>=3)&(sup>=.5))|((core>=4.75+offset)&(s>=3.5)&(spread>=1))|((core+sup>=6+offset)&(s>=2.4))
    active=(s>=2.8)|((core>=5.5+offset)&(s>=2.6))|((core>=4.5+offset)&(s>=2.35))
    active=active|(s>=2.35)
    return (high|provisional)&active
base=select(mat.sum(1));formal=f.loc[base,'TF'].tolist();report('formal',formal)
assert set(formal)==set(['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1'])
sens=[];membership=[]
for j,col in enumerate(cols):
    core=(mat.sum(1)-mat[:,j])*7/6
    hit=select(core);names=f.loc[hit,'TF'].tolist()
    sens.append({'analysis':'leave_one_discovery_layer_out_rescaled','omitted_layer':col,'n_selected':int(hit.sum()),'original_retained':int((hit&base).sum()),'lost':';'.join(f.loc[base&~hit,'TF']),'new':';'.join(f.loc[~base&hit,'TF'])})
    membership.append(hit)
save(pd.DataFrame(sens),'R1_5_leave_one_layer_out.tsv');report('leave_one_layer_out',sens)
for off in [-.5,-.25,0,.25,.5]:
    hit=select(mat.sum(1),off);sens.append({'analysis':'core_total_threshold_shift','omitted_layer':str(off),'n_selected':int(hit.sum()),'original_retained':int((hit&base).sum()),'lost':';'.join(f.loc[base&~hit,'TF']),'new':';'.join(f.loc[~base&hit,'TF'])})
save(pd.DataFrame(sens),'R1_5_sensitivity_summary.tsv')
weight_hits=[]
for i in range(2000):
    w=rng.uniform(.8,1.2,7);w=w*7/w.sum();weight_hits.append(select(mat@w))
freq=np.mean(weight_hits,axis=0)
freqtab=f[['TF','confidence_tier','Core_discovery_score_v2_0_7','max_state_score_0_5','active_contexts']].copy();freqtab['selection_frequency_2000_weight_perturbations']=freq
save(freqtab,'R1_5_weight_sensitivity_all_TFs.tsv');report('formal_weight_frequencies',dict(zip(f.loc[base,'TF'],freq[base])))
save(f[f.TF.isin(['XBP1','ZNF703'])], 'R1_5_XBP1_ZNF703_selection_evidence.tsv')
report('exclusions',f[f.TF.isin(['XBP1','ZNF703'])][['TF',*cols,'Core_discovery_score_v2_0_7','max_state_score_0_5','confidence_tier','active_contexts']].to_dict('records'))
null=[]
for i in range(2000):
    ix=rng.choice(len(f),13,replace=False)
    null.append({'replicate':i+1,'mean_core_score':float(mat[ix].sum(1).mean()),'n_passing_original_rules':int(base[ix].sum())})
save(pd.DataFrame(null),'R1_5_selection_stage_random_sets.tsv')
report('null_scope','2000 random 13-TF sets from 1639; identical stored scoring and selection rules only. Does not rerun downstream raw-omics pipeline or correct targeted acquisition bias. No inferential enrichment claim.')

peaks=read(C/'figure3/raw/figure3d__all_local_peak_RNA_correlations.tsv')
peaks['revision_q_all_saved_peak_tests']=bh(peaks.peak_RNA_p)
save(peaks,'R1_6_all_local_peak_tests_BH.tsv')
current=read(A/'figure3/figure3d.csv')
best=current.merge(peaks[['TF','Peak_Name','revision_q_all_saved_peak_tests']],left_on=['TF','best_local_peak'],right_on=['TF','Peak_Name'],how='left',validate='one_to_one')
best['revision_p_bonferroni_within_locus']=(best.best_local_peak_RNA_p*best.n_local_peaks).clip(upper=1)
best['revision_q_BH_after_within_locus_Bonferroni']=bh(best.revision_p_bonferroni_within_locus)
save(best,'R1_6_peak_selection_audit.tsv')
report('peak_audit',{'all_tests':len(peaks),'source_TFs':int(peaks.TF.nunique()),'selected_TFs':len(best),'all_peak_BH_pass':int((best.revision_q_all_saved_peak_tests<.05).sum()),'within_locus_then_BH_pass':int((best.revision_q_BH_after_within_locus_Bonferroni<.05).sum()),'max_difference_published_global_q':float(np.max(np.abs(best.revision_q_all_saved_peak_tests-best.best_local_peak_RNA_FDR_global)))})

counts=[]
for panel,rel in [('1B','figure1/figure1b_sample_counts.csv'),('1C','figure1/figure1c_sample_counts.csv'),('1F','figure1/figure1f_sample_counts.csv')]:
    d=read(A/rel)
    grouped=d.set_index('SUBTYPE')['n_samples'] if 'n_samples' in d else d.groupby('SUBTYPE').size()
    for group,n in grouped.items():counts.append({'panel':panel,'group':group,'n':int(n)})
clin=read(ROOT/'datasets/TCGA/results/TCGA_CRC_clinical_analysis_dataset.tsv').drop_duplicates('sample_id')
for panel,rel,group,value in [('6A','figure6/figure6a.csv',['score','SUBTYPE'],'value'),('8A','figure8/figure8a.csv',['state_id','pam50'],'state_score')]:
    d=read(A/rel);save(d.groupby(group)[value].agg(['count','median']).reset_index(),f'R1_4_Figure{panel}_medians.tsv')
    normal=d[d[group[1]].eq('Normal')].drop_duplicates('sample_id').merge(clin[['sample_id','sample_type','SUBTYPE']],on='sample_id',suffixes=('','_clinical'))
    save(normal[['sample_id','sample_type','SUBTYPE_clinical' if 'SUBTYPE_clinical' in normal else 'SUBTYPE']],f'R1_4_Figure{panel}_Normal_label_identity.tsv')
    report('Normal_label_'+panel,{'n':len(normal),'sample_type':normal.sample_type.value_counts().to_dict(),'barcode_type':normal.sample_id.str[13:15].value_counts().to_dict()})
sc=read(A/'figure7/figure7a.csv')
sc=sc.drop_duplicates('barcode')
sc_counts=sc.groupby(['celltype_major','subtype']).agg(cells=('barcode','size'),patients=('orig.ident','nunique')).reset_index();save(sc_counts,'R1_7_single_cell_counts.tsv')
report('single_cell',{'all_cells':len(sc),'patients':int(sc['orig.ident'].nunique()),'types':sc.groupby('celltype_major').size().to_dict()})
save(pd.DataFrame(counts),'R1_7_figure_sample_counts.tsv');report('sample_counts',counts)

cox=read(A/'figure8b/figure8b_tcga_os_pfs_cox_ph_truncated_results.tsv');cox['BH_q_across_12']=bh(cox.pvalue);save(cox,'R1_6_Cox_original_with_BH.tsv')
model=read(A/'figure8b/figure8b_tcga_cox_model_input.tsv')
adj=model.merge(clin[['sample_id','AGE','AJCC_PATHOLOGIC_TUMOR_STAGE_clean','SUBTYPE']],on='sample_id',how='left',validate='many_to_one')
stage=adj.AJCC_PATHOLOGIC_TUMOR_STAGE_clean.astype(str)
adj['late_stage']=stage.map(lambda x:0 if x in ['I','IA','IB','II','IIA','IIB'] else (1 if x in ['III','IIIA','IIIB','IIIC','IV'] else np.nan))
adj['age_10y']=pd.to_numeric(adj.AGE,errors='coerce')/10
adj['LumB']=adj.SUBTYPE.eq('BRCA_LumB').astype(int)
res=[]
for (state,endpoint),g in adj.groupby(['state_id','endpoint']):
    cov=['score_z','age_10y','late_stage']+(['LumB'] if state.startswith('Luminal') else [])
    x=g[['time_truncated','event_truncated',*cov]].dropna()
    row={'state_id':state,'endpoint':endpoint,'n':len(x),'events':int(x.event_truncated.sum()),'covariates':';'.join(cov),'stage_encoding':'I/II vs III/IV; X/unknown excluded','status':'ok'}
    try:
        fit=CoxPHFitter(penalizer=0).fit(x,'time_truncated','event_truncated')
        r=fit.summary.loc['score_z'];ph=proportional_hazard_test(fit,x,time_transform='rank').summary
        row.update(hr=float(r['exp(coef)']),lo=float(r['exp(coef) lower 95%']),hi=float(r['exp(coef) upper 95%']),p=float(r.p),score_ph_p=float(ph.loc['score_z','p']))
    except Exception as e:row['status']=str(e)
    res.append(row)
out=pd.DataFrame(res);out['BH_q_across_12']=bh(out.p);save(out,'R1_6_Cox_age_stage_intrinsic_subtype_adjusted.tsv');save(adj,'R1_6_Cox_adjusted_model_input.tsv')
report('cox_BH',cox[['state_id','endpoint','BH_q_across_12']].to_dict('records'))
report('cox_adjusted',out.to_dict('records'))

occ=read(A/'figure4/figure4c.csv');report('occupancy_sources',occ.source_table.unique().tolist())
for source,part in occ.groupby('source_table'):
    if 'regulator_TF' in part and part.regulator_TF.notna().any():
        part=part[part.regulator_TF.isin(formal)]
        rr=[]
        for _,r in part.iterrows():
            targets=[g for g in formal if g!=r.regulator_TF and pd.notna(r.get(g)) and float(r[g])>0]
            rr.append({'regulator_TF':r.regulator_TF,'source':source,'self_occupancy_value':r.get(r.regulator_TF),'n_outgoing_to_other_candidates':len(targets),'targets':';'.join(targets)})
        save(pd.DataFrame(rr),'R1_1_occupancy_'+str(len(summary))+'.tsv');report('occupancy_'+str(len(summary)),rr)

sig=read(R/'figure5-8/figure8/figure8i_scrna_guided_crc_state_signature_table.tsv');sig=sig[sig.recommended_for_bulk_projection.eq(True)]
stored=read(R/'figure5-8/figure7/figure8c_scrna_state_depmap_breast_expression_scores.tsv')
genes=set(f.TF)|set(';'.join(sig.signature_genes).split(';'))
def depmatrix(path):
    headers=pd.read_csv(path,nrows=0).columns
    mapping={c:c.split(' (')[0] for c in headers if c.split(' (')[0] in genes}
    id_col='ModelID' if 'ModelID' in headers else headers[0]
    flags=['IsDefaultEntryForModel'] if 'IsDefaultEntryForModel' in headers else []
    d=read(path,usecols=[id_col,*flags,*mapping]).rename(columns={id_col:'ModelID',**mapping})
    if flags:d=d[d.IsDefaultEntryForModel.eq('Yes')].drop(columns=flags)
    return d.groupby('ModelID').mean(numeric_only=True)
expr=depmatrix(ROOT/'datasets/DepMap/data/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv')
dep=depmatrix(ROOT/'datasets/DepMap/data/CRISPRGeneEffect.csv')
stored=stored[stored.ModelID.isin(expr.index.intersection(dep.index))].copy()
expr=expr.reindex(stored.ModelID);dep=dep.reindex(stored.ModelID)
common=[g for g in f.TF if g in expr and g in dep and expr[g].notna().all() and dep[g].notna().all() and expr[g].std()>0]
expr_z=(expr-expr.mean())/expr.std(ddof=0)
save(expr_z.reset_index(),'R1_6_DepMap_expression_z_51_models.tsv');save(dep.reset_index(),'R1_6_DepMap_gene_effect_51_models.tsv')
bin_expr=pd.qcut(expr[common].mean().rank(method='first'),5,labels=False)
bin_dep=pd.qcut(dep[common].mean().rank(method='first'),5,labels=False)
pool=pd.DataFrame({'gene':common,'expr_bin':bin_expr.values,'dep_bin':bin_dep.values})
save(pool,'R1_6_DepMap_null_matching_pool.tsv')
ctrl=[];nullrows=[]
for _,row in sig.iterrows():
    gs=[g for g in row.signature_genes.split(';') if g in expr_z and g in dep]
    x=stored[row.state_id].to_numpy(float);y=dep[gs].mean(axis=1).to_numpy(float)
    observed=float(spearmanr(x,y).statistic)
    reconstructed=float(spearmanr(expr_z[gs].mean(axis=1),y).statistic)
    nullrho=[];perm=[]
    choices=[]
    for g in gs:
        eb=int(np.clip(np.searchsorted(np.quantile(expr[common].mean(),[.2,.4,.6,.8]),expr[g].mean()),0,4))
        db=int(np.clip(np.searchsorted(np.quantile(dep[common].mean(),[.2,.4,.6,.8]),dep[g].mean()),0,4))
        dist=(pool.expr_bin-eb).abs()+(pool.dep_bin-db).abs()
        candidates=pool.loc[(dist==dist.min())&~pool.gene.isin(gs),'gene'].tolist()
        if not candidates:candidates=[g for g in common if g not in gs]
        choices.append(candidates)
    for b in range(2000):
        selected=[]
        for candidates in choices:
            available=[g for g in candidates if g not in selected]
            selected.append(str(rng.choice(available or candidates)))
        nr=float(spearmanr(expr_z[selected].mean(axis=1),dep[selected].mean(axis=1)).statistic)
        nullrho.append(nr)
        perm.append(float(spearmanr(x,rng.permutation(y)).statistic))
        nullrows.append({'state_id':row.state_id,'replicate':b+1,'random_same_set_rho':nr,'model_label_permuted_rho':perm[-1]})
    ctrl.append({'state_id':row.state_id,'n_models':len(x),'n_genes':len(gs),'observed_rho_mean_GeneEffect':observed,'reconstructed_rho_51_line_zscore':reconstructed,'null_median_rho':float(np.median(nullrho)),'empirical_p_one_sided_matched_gene_sets':(1+np.sum(np.asarray(nullrho)<=reconstructed))/2001,'empirical_p_two_sided_label_permutation':(1+np.sum(np.abs(perm)>=abs(observed)))/2001,'genes':';'.join(gs)})
ct=pd.DataFrame(ctrl);ct['BH_q_matched_gene_sets_6']=bh(ct.empirical_p_one_sided_matched_gene_sets);ct['BH_q_label_permutation_6']=bh(ct.empirical_p_two_sided_label_permutation)
save(ct,'R1_6_DepMap_permutation_controls.tsv');save(pd.DataFrame(nullrows),'R1_6_DepMap_empirical_null_replicates.tsv');report('depmap_controls',ct.to_dict('records'))
report('depmap_null_pool_n',len(common))

manifest=[]
for p in dict.fromkeys(inputs):
    h=hashlib.sha256()
    with p.open('rb') as fh:
        for block in iter(lambda:fh.read(8*1024*1024),b''):h.update(block)
    manifest.append({'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
save(pd.DataFrame(manifest),'input_manifest_sha256.tsv')
(OUT/'analysis_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
print('DONE',flush=True)

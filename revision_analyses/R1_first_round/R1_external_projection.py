from pathlib import Path
import sys, json, hashlib
W=Path(__file__).resolve().parent
import pandas as pd, numpy as np
from scipy.stats import kruskal
P=Path(sys.argv[1]) if len(sys.argv)>1 else W.parents[2];D=P/'core_data/remap_reanalysis/figure5-8'
OUT=Path(sys.argv[2]) if len(sys.argv)>2 else P/'revision_tables';OUT.mkdir(parents=True,exist_ok=True)
sig=pd.read_csv(D/'figure8/figure8i_scrna_guided_crc_state_signature_table.tsv',sep='\t');sig=sig[sig.recommended_for_bulk_projection.eq(True)]
outs=[];cov=[];tests=[];means=[];manifest=[]
for cohort,expr_file,meta_file,id_col,group_col in [
('METABRIC',D/'figure7/raw/METABRIC_scrna_state_signature_genes_zscores_wide.tsv',P/'datasets/METABRIC/data/metabric_sample_clinical_wide.tsv','sampleId','CLAUDIN_SUBTYPE'),
('GSE96058',D/'figure7/raw/GSE96058_scrna_state_signature_genes_expression.tsv',P/'datasets/GSE96058/results/GSE96058_sample_annotation.tsv','sample_id','pam50_subtype')]:
    expr=pd.read_csv(expr_file,sep='\t',index_col=0);meta=pd.read_csv(meta_file,sep='\t').set_index(id_col)
    if cohort=='METABRIC':
        mapping=pd.read_csv(P/'datasets/METABRIC/data/metabric_samples.tsv',sep='\t')
        patient=pd.read_csv(P/'datasets/METABRIC/data/metabric_patient_clinical_wide.tsv',sep='\t')
        sub=mapping[['sampleId','patientId']].merge(patient[['patientId','CLAUDIN_SUBTYPE']],on='patientId').set_index('sampleId')
        meta=meta.join(sub[['CLAUDIN_SUBTYPE']])
    if 'isRepl' in meta:meta=meta[meta.isRepl.astype(str).str.lower().eq('false')]
    common=expr.columns.intersection(meta.index);meta=meta.loc[common].copy();expr=expr[common]
    z=expr.sub(expr.mean(axis=1),axis=0).div(expr.std(axis=1,ddof=0).replace(0,np.nan),axis=0)
    group=meta[group_col].replace({'LumA':'Luminal','LumB':'Luminal','Her2':'HER2-enriched','HER2':'HER2-enriched','Basal':'Basal-like','Normal':'Normal-like'})
    for _,s in sig.iterrows():
        genes=s.signature_genes.split(';');present=[g for g in genes if g in z.index]
        vals=z.loc[present].mean(axis=0)
        out=pd.DataFrame({'cohort':cohort,'sample_id':common,'group':group.to_numpy(),'state_id':s.state_id,'state_score':vals.to_numpy()});outs.append(out)
        cov.append({'cohort':cohort,'state_id':s.state_id,'n_samples_with_expression':len(common),'n_genes_present':len(present),'n_genes_requested':len(genes),'missing_genes':';'.join(set(genes)-set(present))})
        for g,part in out.groupby('group'):means.append({'cohort':cohort,'state_id':s.state_id,'group':g,'n':len(part),'median':part.state_score.median()})
        primary=out[out.group.isin(['Luminal','HER2-enriched','Basal-like','Normal-like'])]
        x=[p.state_score.dropna().to_numpy() for _,p in primary.groupby('group')];stat,p=kruskal(*x)
        tests.append({'cohort':cohort,'state_id':s.state_id,'n_classified':len(primary),'kruskal_H':stat,'p':p})
    for path in [expr_file,meta_file]:manifest.append({'path':path.relative_to(P).as_posix(),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
for data,name in [(pd.concat(outs),'R1_7_external_state_scores.tsv'),(pd.DataFrame(cov),'R1_7_external_gene_coverage.tsv'),(pd.DataFrame(means),'R1_7_external_state_medians.tsv'),(pd.DataFrame(tests),'R1_7_external_state_tests.tsv'),(pd.DataFrame(manifest),'external_projection_inputs_sha256.tsv')]:
    data.to_csv(OUT/name,sep='\t',index=False)
    if name!='R1_7_external_state_scores.tsv':print(name);print(data.to_string(index=False))

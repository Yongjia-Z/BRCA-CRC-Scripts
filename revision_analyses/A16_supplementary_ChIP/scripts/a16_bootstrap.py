from a16_common import *
import pandas as pd
for d in ['raw/metadata','raw/peaks','raw/signals','derived/hg38','results','logs','tests']:(A16/d).mkdir(parents=True,exist_ok=True)
save_tsv(A16/'results/A16_promoter_windows.tsv',promoter_table())
cu=pd.read_csv(ROOT/'data/dataset_lists/A14_PGR_curation.tsv',sep='\t');kept=cu[cu.keep.astype(str).str.lower().eq('true')];rows=[]
for r in kept.itertuples(index=False):
    f=A14/'raw_chipatlas_PGR'/f'{r.srx}.05.bed';n,bad,h=peak_hits(f)
    rows.append({'srx':r.srx,'cell_line':r.cell_type,'n_peaks':n,'bad_rows':bad,'self_peaks':h['PGR'],'path':str(f.relative_to(ROOT)),'sha256':sha(f),**h})
save_tsv(A16/'results/A16_baseline_PGR_by_dataset.tsv',rows)
stored=pd.read_csv(A14/'A14_chipatlas_PGR_by_dataset.tsv',sep='\t').set_index('srx');now=pd.DataFrame(rows).set_index('srx')
equal=all((now.loc[stored.index,t].to_numpy()==stored[t].to_numpy()).all() for t in C16+['PGR'])
assert len(rows)==77 and sum(r['self_peaks']>0 for r in rows)==22 and equal
score=pd.read_csv(PRIOR/'A7_Basal_TNBC_H3K27ac_layer/A7_all_tf_scores_corrected_basal_H3K27ac.csv').set_index('TF');assert len(score)==1639
status={'baseline':'submitted','PGR_datasets':len(rows),'PGR_with_self_peak':sum(r['self_peaks']>0 for r in rows),'PGR_total_self_peaks':sum(r['self_peaks'] for r in rows),'PGR_all_promoter_counts_match_A14':equal,'score_table_TFs':len(score),'checked_utc':stamp()}
save_json(A16/'results/A16_baseline_verification.json',status);print(json.dumps(status,ensure_ascii=False),flush=True)

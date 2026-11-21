from a16_common import *
import pandas as pd, numpy as np, re

rows=[];coordinate_audit=[]
legacy=promoter_table();zero=[{**r,'start':r['start']-1,'end':r['end']-1} for r in legacy]

def assay(rec):
    f=Path(rec['path']);f=f if f.is_absolute() else ROOT/f
    n,bad,h=peak_hits(f,legacy);n2,bad2,h2=peak_hits(f,zero)
    assert n==n2 and bad==bad2
    changes=[t for t in h if h[t]!=h2[t]]
    coordinate_audit.append({'dataset_id':rec['dataset_id'],'TF':rec['TF'],'path':str(f.relative_to(ROOT)),'n_changed_promoters_1bp_conversion':len(changes),'changed_promoters':';'.join(changes),'legacy_counts':json.dumps({t:h[t] for t in changes}),'BED0_counts':json.dumps({t:h2[t] for t in changes})})
    if n==0 and rec['source']=='ChIP-Atlas supplemental':
        rec={**rec,'decision':'unresolved_no_peaks','use_strict':False}
    return {**rec,'path':str(f.relative_to(ROOT)),'sha256':sha(f),'bytes':f.stat().st_size,'n_peaks':n,'malformed_rows':bad,'self_peaks':h[rec['TF']],'quality_flag':'no usable peak result; unresolved, not a negative assay' if n==0 else ('low peak count (<3000); not an exclusion or biological-negative criterion' if n<3000 else ''),**h}

st=pd.read_csv(REMAP/'results/remap_breast_tier1_download_status.tsv',sep='\t').fillna('')
cu=pd.read_csv(ROOT/'data/dataset_lists/A14_remap_curation.tsv',sep='\t').fillna('');keep=set(cu[cu.keep.astype(str).str.lower().eq('true')].dataset_file)
for x in st.to_dict('records'):
    if x['dataset_file'] not in keep or not Path(x['local_peak_file']).exists():continue
    series=x['experiment']+' | GSE170028' if x['target_name']=='SPDEF' and x['experiment']=='ENCSR042GSX' else x['experiment']
    rows.append(assay({'dataset_id':x['dataset_file'],'TF':x['target_name'],'source':'ReMap baseline','path':x['local_peak_file'],'cell_line':x['biotype'],'context':CONTEXT.get(x['biotype'],''),'series':series,'gsm':'','srx':'','condition':x['biotype_modification'],'decision':'baseline_main','use_strict':True,'processing_note':'historical curated baseline; endogenous SPDEF tagging retained after review; ENCSR042GSX and GSE170028 identify the same experiment','url':x['download_url']}))

pgr=pd.read_csv(A16/'results/A16_baseline_PGR_by_dataset.tsv',sep='\t').fillna('')
pmeta=pd.read_csv(ROOT/'data/dataset_lists/A14_chipatlas_PGR_manifest.tsv',sep='\t').fillna('').set_index('srx')
for x in pgr.to_dict('records'):
    rows.append(assay({'dataset_id':x['srx'],'TF':'PGR','source':'PGR curated ChIP-Atlas','path':x['path'],'cell_line':x['cell_line'],'context':'Luminal','series':'','gsm':'','srx':x['srx'],'condition':pmeta.loc[x['srx'],'title'],'decision':'baseline_main','use_strict':True,'processing_note':'77 previously curated profiles, not newly collected','url':f"https://chip-atlas.dbcls.jp/data/hg38/eachData/bed05/{x['srx']}.05.bed"}))

tn=pd.read_csv(REMAP/'results/TNBC_direct_binding/tnbc_direct_binding_dataset_manifest.tsv',sep='\t').fillna('')
for x in tn[tn.usable_for_hg38_overlap.astype(str).str.lower().eq('true')].to_dict('records'):
    rows.append(assay({'dataset_id':x['file_name'],'TF':x['regulator_TF'],'source':'TNBC baseline author/ReMap','path':x['local_peak_file'],'cell_line':x['cell_line'],'context':'Basal/TNBC','series':x['dataset'],'gsm':'','srx':'','condition':'','decision':'baseline_main','use_strict':True,'processing_note':'same study/condition may be represented in another processing source; not independent','url':f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={x['dataset']}"}))

m=pd.read_csv(ROOT/'data/dataset_lists/A16_peak_file_manifest.tsv',sep='\t').fillna('')
for x in m[m.status.ne('failed')].to_dict('records'):
    rows.append(assay({'dataset_id':x['srx'],'TF':x['TF'],'source':'ChIP-Atlas supplemental','path':x['path'],'cell_line':x['cell_canonical'],'context':CONTEXT.get(x['cell_canonical'],''),'series':x['geo_series'],'gsm':x['gsm'],'srx':x['srx'],'condition':x['title'],'decision':x['decision'],'use_strict':x['decision'].startswith('main'),'processing_note':'q<1e-5; repeats of an existing GSE are source comparisons, not new independent studies','url':x['url']}))

lift=json.loads((A16/'results/A16_liftover_verification.json').read_text(encoding='utf-8'))
extra=ROOT/'data/dataset_lists/A16_extra_liftover_verification.json'
if extra.exists():lift+=json.loads(extra.read_text(encoding='utf-8'))
for x in lift:
    tf=Path(x['dest_path']).name.split('_')[0];gsm_match=re.search(r'GSM\d+',x['dest_path']);gsm=gsm_match.group(0) if gsm_match else ''
    if gsm:
        s=(A16/'raw/metadata/geo'/f'{gsm}.soft.txt').read_text(encoding='utf-8');series=' | '.join(re.findall(r'^!Sample_series_id = (.*)',s,re.M))
    else:series='GSE75070'
    rows.append(assay({'dataset_id':gsm or 'GSE75070_author_pooled_peaks','TF':tf,'source':'GEO author supplemental','path':x['dest_path'],'cell_line':'MCF-7','context':'Luminal','series':series,'gsm':gsm,'srx':'SRX1434852;SRX1434854' if not gsm else '', 'condition':'parental MCF-7; author peak table','decision':'main','use_strict':True,'processing_note':f"{x['source_build']} to hg38; >=95% aligned bases; author thresholds retained; same biological assay across processing sources counted separately only in source comparison",'url':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm or "GSE75070"}'}))

d=pd.DataFrame(rows)
d.to_csv(A16/'results/A16_occupancy_by_dataset.tsv',sep='\t',index=False,encoding='utf-8-sig')
save_tsv(A16/'results/A16_coordinate_1bp_sensitivity.tsv',coordinate_audit)

def mat(frame,mode):
    if mode=='peak_sum':return frame.groupby('TF')[TFS].sum().reindex(TFS).fillna(0)
    q=frame[TFS].gt(0);q['TF']=frame.TF.to_numpy()
    return q.groupby('TF')[TFS].sum().reindex(TFS).fillna(0)

sources={'ReMap baseline':d.source.eq('ReMap baseline'), 'PGR curated ChIP-Atlas':d.source.eq('PGR curated ChIP-Atlas'),'TNBC baseline author/ReMap':d.source.eq('TNBC baseline author/ReMap')}
baseline=mat(d[sources['ReMap baseline']],'peak_sum')+8*mat(d[sources['TNBC baseline author/ReMap']],'peak_sum')
baseline.loc['PGR']=mat(d[sources['PGR curated ChIP-Atlas']],'peak_sum').loc['PGR']
baseline=np.log1p(baseline)
stored=pd.read_csv(A14/'A14_occupancy_matrix_baseline.tsv',sep='\t').set_index('regulator_TF').loc[TFS,TFS]
err=float(np.max(np.abs(baseline.to_numpy()-stored.to_numpy())))
baseline.to_csv(A16/'results/A16_baseline_matrix_replay.tsv',sep='\t',encoding='utf-8-sig')
save_json(A16/'results/A16_baseline_matrix_check.json',{'max_abs_difference':err,'matches_tolerance_1e-6':err<1e-6,'n_files_recomputed':len(d[d.source.isin(sources)]),'all_13x13_entries_checked':True})
assert err<1e-6,err

summary=[]
for source in ['ChIP-Atlas supplemental','GEO author supplemental']:
    for stratum in ['main_unsorted','main_with_sorted','support']:
        if stratum=='main_unsorted':sel=d.source.eq(source)&d.decision.eq('main')
        elif stratum=='main_with_sorted':sel=d.source.eq(source)&d.use_strict
        else:sel=d.source.eq(source)&d.decision.str.startswith('support')
        f=d[sel];tag=('ChIPAtlas' if source.startswith('ChIP') else 'GEOauthor')+'_'+stratum
        for mode in ['peak_sum','profiles_with_peak']:
            mx=mat(f,mode);mx.loc[~mx.index.isin(f.TF)]=np.nan
            mx.index.name='regulator_TF';mx.to_csv(A16/'results'/f'A16_{tag}_{mode}.tsv',sep='\t',encoding='utf-8-sig')
        for tf,sub in f.groupby('TF'):
            summary.append({'source':source,'stratum':stratum,'TF':tf,'n_profiles':len(sub),'n_with_self_peak':int(sub.self_peaks.gt(0).sum()),'total_self_peaks':int(sub.self_peaks.sum()),'min_peaks':int(sub.n_peaks.min()),'max_peaks':int(sub.n_peaks.max()),'cells':';'.join(sorted(sub.cell_line.unique())),'series':';'.join(sorted(set(g.strip() for s in sub.series for g in re.split(r'[|;]',s) if g.strip()))),'status':'self-locus occupancy observed' if sub.self_peaks.gt(0).any() else 'no promoter peak in these profiles; no biological absence inference'})
save_tsv(A16/'results/A16_self_occupancy_summary.tsv',summary)

strict=d[d.use_strict].copy();u=mat(strict,'profiles_with_peak').gt(0).astype(float)
u.loc[~u.index.isin(strict.TF)]=np.nan;u.index.name='regulator_TF'
u.to_csv(A16/'results/A16_strict_evidence_union_binary.tsv',sep='\t',encoding='utf-8-sig')
new_pairs=[]
for a in TFS:
    for b in TFS:
        if a==b:continue
        supporting=strict[(strict.TF==a)&strict[b].gt(0)]
        if len(supporting) and stored.loc[a,b]==0:
            new_pairs.append({'regulator_TF':a,'target_TF':b,'supporting_sources':';'.join(sorted(supporting.source.unique())),'series':';'.join(sorted(supporting.series.unique())),'dataset_ids':';'.join(supporting.dataset_id),'binding_only':'not a full edge score or causal validation'})
save_tsv(A16/'results/A16_new_binding_pairs.tsv',new_pairs,fields=['regulator_TF','target_TF','supporting_sources','series','dataset_ids','binding_only'])
print(json.dumps({'total_dataset_source_rows':len(d),'baseline_matrix_error':err,'new_binding_pairs':len(new_pairs),'coordinate_1bp_changed_files':sum(r['n_changed_promoters_1bp_conversion']>0 for r in coordinate_audit)},ensure_ascii=False),flush=True)

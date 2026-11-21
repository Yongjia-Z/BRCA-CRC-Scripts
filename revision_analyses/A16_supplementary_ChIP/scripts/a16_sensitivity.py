from a16_common import *
import pandas as pd, numpy as np, networkx as nx, re

CTX=['Luminal','HER2','Basal/TNBC']

def call(scores,core,support):
    s=dict(zip(CTX,map(float,scores)));mx=max(s.values());sec=sorted(s.values())[-2]
    primary=[k for k,v in s.items() if v>=mx-.20 and v>=2.80]
    if not primary and core>=5.50 and mx>=2.60:primary=[max(s,key=s.get)]
    if not primary and core>=4.50 and mx>=2.35:primary=[max(s,key=s.get)]
    secondary=[k for k,v in s.items() if k not in primary and v>=2.35 and v>=.82*mx]
    if core>=5.50 and mx>=2.60:tier='High-confidence'
    elif (core>=4.75 and mx>=3 and support>=.50) or (core>=4.75 and mx>=3.5 and mx-sec>=1) or (core+support>=6 and mx>=2.4):tier='Provisional'
    else:tier='Supporting'
    return tier, ';'.join(primary), ';'.join(primary+secondary)

def serials(s):return set(re.findall(r'GSE\d+|ENCSR[A-Z0-9]+',str(s)))

if __name__=='__main__':
    a=pd.read_csv(PRIOR/'A7_Basal_TNBC_H3K27ac_layer/A7_all_tf_scores_corrected_basal_H3K27ac.csv').set_index('TF').loc[TFS]
    d=pd.read_csv(A16/'results/A16_occupancy_by_dataset.tsv',sep='\t').fillna('')
    d['use_strict']=d.use_strict.astype(str).str.lower().eq('true')
    d['count16']=d[C16].gt(0).sum(axis=1)
    base=d[d.source.eq('ReMap baseline')&d.use_strict];counts=base.groupby('TF').count16.sum().reindex(TFS).fillna(0)
    denom=float(counts.max());assert denom==1033,denom
    scale=lambda x: np.log1p(x.astype(float))/np.log1p(denom)
    old_b=a.Breast_context_direct_binding.copy();old_t=a.TNBC_direct_binding.copy()
    clean_b=scale(counts)
    pgr_count=d[d.source.eq('PGR curated ChIP-Atlas')].count16.sum();assert pgr_count==171
    author=d[d.source.eq('GEO author supplemental')&d.use_strict]
    prior_sources=d[d.use_strict&d.source.isin(['ReMap baseline','TNBC baseline author/ReMap','PGR curated ChIP-Atlas'])]
    seen={tf:set().union(*(serials(s) for s in prior_sources.loc[prior_sources.TF.eq(tf),'series'])) for tf in TFS}
    records=[];scenarios={};network_frames={}
    def add(name,b,t,note,frame=None):
        scenarios[name]=(b.copy(),t.copy(),note)
        if frame is not None:network_frames[name]=frame.copy()
    base_strict=d[d.use_strict&d.source.isin(['ReMap baseline','PGR curated ChIP-Atlas','TNBC baseline author/ReMap'])]
    add('S0_submitted',old_b,old_t,'unchanged manuscript baseline')
    add('S1_curated_baseline',clean_b,old_t,'34 prior exclusions applied; endogenous SPDEF tagging retained after review; original 16 target promoters',base_strict)
    tnbc_base=d[d.source.eq('TNBC baseline author/ReMap')]
    tnbc_counts=tnbc_base.groupby('TF').n_peaks.sum().reindex(TFS).fillna(0);tnbc_denom=float(tnbc_counts.max())
    overlap_review=[]
    for sorted_ok in [False,True]:
        tag='with_sorted' if sorted_ok else 'unsorted'
        ca=d[d.source.eq('ChIP-Atlas supplemental')&d.use_strict&(True if sorted_ok else d.decision.ne('main_condition_sorted'))&d.n_peaks.gt(0)].copy()
        extra=ca[[not bool(serials(r.series)&seen[r.TF]) for r in ca.itertuples()]].copy()
        for r in ca.itertuples():
            overlap_review.append({'scenario_stratum':tag,'TF':r.TF,'dataset_id':r.dataset_id,'series':r.series,'already_used_series':';'.join(sorted(serials(r.series)&seen[r.TF])),'counts_as_new_study_source':not bool(serials(r.series)&seen[r.TF]),'note':'same-series profiles are alternate processing or repeats, not new independent studies'})
        n=counts.add(extra.groupby('TF').count16.sum().reindex(TFS).fillna(0));n.loc['PGR']=pgr_count
        n.loc['PBX1']=author[author.TF.eq('PBX1')].count16.sum()
        extra_t=extra[extra.context.eq('Basal/TNBC')]
        nt=tnbc_counts.add(extra_t.groupby('TF').n_peaks.sum().reindex(TFS).fillna(0))
        nt_scale=np.log1p(nt)/np.log1p(max(tnbc_denom,float(nt.max())))
        tn_new=old_t.copy()
        for tf in extra_t.TF.unique():tn_new.loc[tf]=nt_scale.loc[tf]
        fr=pd.concat([base_strict,extra,author[author.TF.eq('PBX1')]],ignore_index=True)
        add(f'S2_new_studies_{tag}',scale(n),tn_new,'add only previously unrepresented GEO series; known same-series assays not counted twice; fixed clean-baseline breast denominator 1033',fr)
        nr=counts.copy()
        for tf in ca.TF.unique():
            sub=ca[ca.TF.eq(tf)];new_gses=set().union(*(serials(s) for s in sub.series))
            removed=base[base.TF.eq(tf)&base.series.map(lambda s: bool(serials(s)&new_gses))].count16.sum()
            nr.loc[tf]=counts.loc[tf]-removed+sub.count16.sum()
        nr.loc['PGR']=pgr_count;nr.loc['PBX1']=author[author.TF.eq('PBX1')].count16.sum()
        cat=ca[ca.context.eq('Basal/TNBC')];ct=cat.groupby('TF').n_peaks.sum()
        source_t=old_t.copy();source_denom=max(float(ct.max()) if len(ct) else 0,1)
        for tf in cat.TF.unique():source_t.loc[tf]=np.log1p(ct.loc[tf])/np.log1p(source_denom)
        def replace_old_group(r):
            sub=ca[ca.TF.eq(r.TF)]
            gs=set().union(*(serials(s) for s in sub.series)) if len(sub) else set()
            return bool(serials(r.series)&gs)
        fr3=base_strict[[not replace_old_group(r) for r in base_strict.itertuples()]]
        fr3=pd.concat([fr3,ca,author[author.TF.eq('PBX1')]],ignore_index=True)
        add(f'S3_source_replacement_{tag}',scale(nr),source_t,'replace matching study groups by reviewed Atlas profiles; TNBC peak totals max-scaled within that source; incomplete source coverage is not a fair all-TF rescreen',fr3)
        na=n.copy()
        for tf in author.TF.unique():na.loc[tf]=author[author.TF.eq(tf)].count16.sum()
        fr4=pd.concat([fr[~fr.TF.isin(author.TF.unique())],author],ignore_index=True)
        add(f'S4_author_peak_alternative_{tag}',scale(na),tn_new,'for PBX1/RUNX1/SPDEF replace by author peak tables; original author significance thresholds retained',fr4)
    for name,(b,t,note) in scenarios.items():
        for tf in TFS:
            oldmax=float(a.loc[tf,'Direct_binding_MAX']);newmax=max(float(b.loc[tf]),float(t.loc[tf]))
            core=float(a.loc[tf,'Core_discovery_score_v2_0_7'])-oldmax+newmax
            scores=[float(a.loc[tf,'state_score_'+ctx])+(float(b.loc[tf])-float(old_b.loc[tf]) if ctx!='Basal/TNBC' else float(t.loc[tf])-float(old_t.loc[tf])) for ctx in CTX]
            tier,primary,active=call(scores,core,float(a.loc[tf,'Supportive_score_0_3']))
            records.append({'scenario':name,'TF':tf,'old_core':a.loc[tf,'Core_discovery_score_v2_0_7'],'new_core':core,'old_binding':oldmax,'new_binding':newmax,'new_breast_binding':b.loc[tf],'new_TNBC_binding':t.loc[tf], 'old_tier':a.loc[tf,'confidence_tier'],'new_tier':tier,'primary':primary,'active':active,**{f'context_{ctx}':v for ctx,v in zip(CTX,scores)},'scope':note})
    out=pd.DataFrame(records);out.to_csv(A16/'results/A16_tier_sensitivity.tsv',sep='\t',index=False,encoding='utf-8-sig')
    save_tsv(A16/'results/A16_series_overlap_review.tsv',overlap_review)
    assert (out[out.scenario.eq('S0_submitted')].old_tier==out[out.scenario.eq('S0_submitted')].new_tier).all()
    out.groupby('scenario').apply(lambda x: pd.Series({'n_high':int(x.new_tier.eq('High-confidence').sum()),'n_provisional':int(x.new_tier.eq('Provisional').sum()),'n_changed':int(x.old_tier.ne(x.new_tier).sum()),'changed_TFs':';'.join(x.loc[x.old_tier.ne(x.new_tier),'TF'])}),include_groups=False).to_csv(A16/'results/A16_tier_sensitivity_summary.tsv',sep='\t',encoding='utf-8-sig')
    e=pd.read_csv(ROOT/'core_data/figure2/raw/figure2b__TCGA_results__TF_TF_regulatory_edge_evidence_score_table_with_TRRUST.tsv',sep='\t')
    e=e[e.regulator_TF.isin(TFS)&e.target_TF.isin(TFS)&e.regulator_TF.ne(e.target_TF)].copy();assert len(e)==45
    bits=['edge_CancerTarget_RP_evidence','edge_expression_correlation_evidence','edge_target_ATAC_peak_to_gene_evidence','edge_regulator_ATAC_motif_evidence','edge_regulator_ATAC_footprint_evidence','edge_target_H3K27ac_evidence','edge_activity_context_evidence','TRRUST_curated_edge_evidence']
    assert not e[bits].isna().any().any()
    assert np.array_equal(e[bits].astype(int).sum(axis=1),e.TF_TF_edge_evidence_score_with_TRRUST)
    historical=pd.read_csv(A14/'A14_occupancy_matrix_baseline.tsv',sep='\t').set_index('regulator_TF')
    network_rows=[];summ=[];motif_rows=[];source_members=[]
    def selected(reg,tgt,score,binding):
        if score>=7:return True
        if score<6 or not binding:return False
        for ctx in CTX:
            same=int(a.loc[reg,ctx+'_member']) and int(a.loc[tgt,ctx+'_member'])
            primary=(ctx in str(a.loc[reg,'primary_context']).split(';')) or (ctx in str(a.loc[tgt,'primary_context']).split(';'))
            if same and (primary or (ctx=='HER2' and 'TFAP2C' in [reg,tgt])):return True
        return False
    for name in scenarios:
        fr=network_frames.get(name)
        if fr is None:binary=historical.loc[TFS,TFS].gt(0)
        else:
            binary=fr.groupby('TF')[TFS].sum().reindex(TFS).fillna(0).gt(0)
            for r in fr.itertuples(index=False):
                source_members.append({'scenario':name,'TF':r.TF,'dataset_id':r.dataset_id,'source':r.source,'series':r.series,'path':r.path,'sha256':r.sha256})
        G=nx.DiGraph();G.add_nodes_from(TFS);ge6=0
        for r in e.itertuples(index=False):
            bind=bool(binary.loc[r.regulator_TF,r.target_TF]);score=int(r.TF_TF_edge_evidence_score_with_TRRUST)+int(bind)
            take=selected(r.regulator_TF,r.target_TF,score,bind);ge6+=score>=6
            network_rows.append({'scenario':name,'regulator_TF':r.regulator_TF,'target_TF':r.target_TF,'nonbinding_score_0_8':r.TF_TF_edge_evidence_score_with_TRRUST,'binding':bind,'score_0_9':score,'in_network':take,'context_caution':'historical network uses breast-context pooled binding; new MCF-7 evidence does not demonstrate HER2/TNBC-specific regulation'})
            if take:G.add_edge(r.regulator_TF,r.target_TF,weight=score)
        pairs=sorted({tuple(sorted((x,y))) for x,y in G.edges if G.has_edge(y,x)})
        cycles=sorted({tuple(c[c.index(min(c)):]+c[:c.index(min(c))]) for c in nx.simple_cycles(G) if len(c)==3})
        ffl=sorted({(u,v,w) for u,v in G.edges for w in G.successors(v) if w!=u and G.has_edge(u,w)})
        pr=nx.pagerank(G,weight='weight');summ.append({'scenario':name,'n_complete_pairs':len(e),'n_edges_ge6':ge6,'n_network_edges':G.number_of_edges(),'reciprocal_pairs':len(pairs),'three_node_cycles':len(cycles),'feedforward_loops':len(ffl),'top_pagerank':max(pr,key=pr.get)})
        for typ,objs in [('reciprocal',pairs),('cycle3',cycles),('feedforward',ffl)]:
            for obj in objs:motif_rows.append({'scenario':name,'type':typ,'nodes':';'.join(obj),'meaning':'inferred topology under fixed historical rule, not functional feedback'})
    pd.DataFrame(network_rows).to_csv(A16/'results/A16_network_edge_sensitivity.tsv',sep='\t',index=False,encoding='utf-8-sig')
    save_tsv(A16/'results/A16_network_summary.tsv',summ);save_tsv(A16/'results/A16_network_motifs.tsv',motif_rows)
    save_tsv(A16/'results/A16_scenario_source_membership.tsv',source_members)
    s0=next(x for x in summ if x['scenario']=='S0_submitted');assert s0['n_network_edges']==24 and s0['n_edges_ge6']==35,s0
    save_json(A16/'results/A16_sensitivity_parameters.json',{'breast_normalization_denominator':denom,'TNBC_baseline_normalization_denominator':tnbc_denom,'score_promoters':C16,'changing_layers':'binding only; supportive and other discovery layers fixed','inference_scope':'conditional replay on original 13 TFs only; no global selection-bias correction','endogenous_tag_decision':'SPDEF ENCSR042GSX/GSE170028 retained after review; tagging is not treated as functional perturbation','S0_network_reproduced':True,'checked_utc':stamp()})
    print(out.groupby('scenario').new_tier.value_counts().to_string());print(pd.DataFrame(summ).to_string(index=False))

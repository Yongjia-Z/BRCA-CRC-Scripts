import sys, re
from pathlib import Path
import numpy as np, pandas as pd, networkx as nx
W, C, A16, O = map(Path, sys.argv[1:5]); O.mkdir(parents=True, exist_ok=True)
pd.set_option('display.width', 230, 'display.max_columns', 40, 'display.max_rows', 200, 'display.max_colwidth', 80)
T = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']
A14 = C/'A14_occupancy_and_network_checks'
base = pd.read_csv(A14/'A14_occupancy_matrix_baseline.tsv', sep='\t').set_index('regulator_TF').loc[T, T].astype(float)
d = pd.read_csv(A16/'results/A16_occupancy_by_dataset.tsv', sep='\t', encoding='utf-8-sig').fillna('')
d['use_strict'] = d.use_strict.astype(str).str.lower().eq('true')
serials = lambda s: set(re.findall(r'GSE\d+|ENCSR[A-Z0-9]+', str(s)))
prior = d[d.use_strict & d.source.isin(['ReMap baseline', 'TNBC baseline author/ReMap', 'PGR curated ChIP-Atlas'])]
seen = {tf: set().union(*([serials(s) for s in prior.loc[prior.TF.eq(tf), 'series']] or [set()])) for tf in T}
ca = d[d.source.eq('ChIP-Atlas supplemental') & d.use_strict & d.decision.ne('main_condition_sorted') & (d.n_peaks.astype(float) > 0)].copy()
extra = ca[[not bool(serials(r.series) & seen[r.TF]) for r in ca.itertuples()]].copy()
pbx = d[d.source.eq('GEO author supplemental') & d.use_strict & d.TF.eq('PBX1')].copy()
new = pd.concat([extra, pbx], ignore_index=True); new.to_csv(O/'A19_S2_new_study_profiles_used.tsv', sep='\t', index=False)
print(new.groupby('TF').agg(profiles=('dataset_id', 'size'), series=('series', lambda s: len(set().union(*(serials(x) for x in s))))))
M = np.expm1(base).copy(); add = pd.DataFrame(0.0, index=T, columns=T)
for r in new.itertuples():
    w = 8.0 if r.context == 'Basal/TNBC' else 1.0
    for t in T: add.loc[r.TF, t] += w*float(getattr(r, t))
M2 = M + add; m2 = np.log1p(M2); m2.index.name = 'regulator_TF'; m2.to_csv(O/'A19_occupancy_matrix.tsv', sep='\t')
changed = [(a, b, round(base.loc[a, b], 3), round(m2.loc[a, b], 3)) for a in T for b in T if abs(base.loc[a, b] - m2.loc[a, b]) > 1e-9]
newcells = [(a, b) for a in T for b in T if a != b and base.loc[a, b] == 0 and m2.loc[a, b] > 0]
print('changed cells', len(changed), 'new off-diagonal relations', newcells)
print('outgoing occupancy counts (other 12 promoters):', {a: int((m2.loc[a].drop(a) > 0).sum()) for a in T})
un = pd.read_csv(A16/'results/A16_strict_evidence_union_binary.tsv', sep='\t', encoding='utf-8-sig').set_index('regulator_TF').loc[T, T]
alt = [(a, b) for a in T for b in T if a != b and float(un.loc[a, b] if pd.notna(un.loc[a, b]) else 0) > 0 and m2.loc[a, b] == 0]
print('relations only with alternative processing of the same studies', alt)
pd.DataFrame(alt, columns=['regulator', 'target']).to_csv(O/'A19_alternative_processing_only_relations.tsv', sep='\t', index=False)
s = pd.read_csv(A16/'results/A16_self_occupancy_summary.tsv', sep='\t', encoding='utf-8-sig')
su = pd.read_csv(A14/'A14_self_occupancy_support.tsv', sep='\t').set_index('TF')
rows = []
for t in T:
    rec = dict(TF=t, earlier_datasets=int(su.loc[t, 'n_datasets']) if pd.notna(su.loc[t, 'n_datasets']) else 0, earlier_with_self_peak=int(su.loc[t, 'datasets_with_self_peak']) if pd.notna(su.loc[t, 'datasets_with_self_peak']) else 0)
    for lab, src in [('atlas', 'ChIP-Atlas supplemental'), ('author', 'GEO author supplemental')]:
        x = s[(s.TF == t) & (s.source == src) & (s.stratum == 'main_unsorted')]
        rec[lab + '_profiles'] = int(x.n_profiles.sum()); rec[lab + '_with_self_peak'] = int(x.n_with_self_peak.sum()); rec[lab + '_max_peaks'] = int(x.max_peaks.max()) if len(x) else 0
    rows.append(rec)
S = pd.DataFrame(rows); print(S.to_string()); S.to_csv(O/'A19_self_locus_by_source.tsv', sep='\t', index=False)
e = pd.read_csv(A16/'results/A16_network_edge_sensitivity.tsv', sep='\t', encoding='utf-8-sig'); e2 = e[e.scenario == 'S2_new_studies_unsorted'].copy(); e1 = e[e.scenario == 'S1_curated_baseline']
ev = pd.read_csv(W/'core_data/figure2/raw/figure2b__TCGA_results__TF_TF_regulatory_edge_evidence_score_table_with_TRRUST.tsv', sep='\t')
ev = ev[ev.regulator_TF.isin(T) & ev.target_TF.isin(T) & ev.regulator_TF.ne(ev.target_TF)]
BITS = {'TRRUST': [c for c in ev.columns if 'TRRUST' in c and 'score' not in c.lower()][0], 'CancerTarget': 'edge_CancerTarget_RP_evidence', 'Expr corr.': 'edge_expression_correlation_evidence', 'Target ATAC': 'edge_target_ATAC_peak_to_gene_evidence', 'Reg. motif': 'edge_regulator_ATAC_motif_evidence', 'Reg. footprint': 'edge_regulator_ATAC_footprint_evidence', 'Target H3K27ac': 'edge_target_H3K27ac_evidence', 'Activity': [c for c in ev.columns if c.startswith('edge_activity')][0]}
print(BITS)
tab = e2.merge(ev[['regulator_TF', 'target_TF'] + list(BITS.values())], on=['regulator_TF', 'target_TF']).rename(columns={v: k for k, v in BITS.items()}); tab['Direct'] = tab.binding.astype(int)
tab = tab[['regulator_TF', 'target_TF', 'score_0_9', 'in_network', 'TRRUST', 'Direct', 'CancerTarget', 'Expr corr.', 'Target ATAC', 'Reg. motif', 'Reg. footprint', 'Target H3K27ac', 'Activity']]
assert (tab[['TRRUST', 'Direct', 'CancerTarget', 'Expr corr.', 'Target ATAC', 'Reg. motif', 'Reg. footprint', 'Target H3K27ac', 'Activity']].astype(int).sum(1) == tab.score_0_9).all()
tab.sort_values(['score_0_9', 'regulator_TF'], ascending=[False, True]).to_csv(O/'A19_edges.tsv', sep='\t', index=False)
net = tab[tab.in_network]; G = nx.DiGraph(); G.add_nodes_from(T)
for r in net.itertuples(): G.add_edge(r.regulator_TF, r.target_TF, weight=r.score_0_9)
pr = nx.pagerank(G, weight='weight'); prs = pd.Series(pr).reindex(T); print(prs.round(3).sort_values(ascending=False).to_dict(), 'edges', G.number_of_edges(), 'ge6', int((tab.score_0_9 >= 6).sum()))
print('new in network vs S1:', sorted(set(map(tuple, net[['regulator_TF', 'target_TF']].values)) - set(map(tuple, e1[e1.in_network][['regulator_TF', 'target_TF']].values))))
print('outdegree', dict(G.out_degree()))
LAY = ['TRRUST', 'Direct', 'CancerTarget', 'Expr corr.', 'Target ATAC', 'Reg. motif', 'Reg. footprint', 'Target H3K27ac', 'Activity']
tt = tab.set_index(['regulator_TF', 'target_TF'])
CTX = {'Luminal': {'ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C'}, 'HER2': {'FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1'}, 'Basal/TNBC': {'RUNX1','TP63','FOXC1'}}
def motifs(g, label):
    out = []
    for a, b in sorted({tuple(sorted(x)) for x in g.edges if g.has_edge(x[1], x[0])}): out.append(('REC', f'{a} ↔ {b}', [(a, b), (b, a)]))
    for c in sorted({tuple(c[c.index(min(c)):] + c[:c.index(min(c))]) for c in nx.simple_cycles(g) if len(c) == 3}): out.append(('CYC', ' → '.join(c) + f' → {c[0]}', [(c[0], c[1]), (c[1], c[2]), (c[2], c[0])]))
    for u, v, w_ in sorted({(u, v, w) for u, v in g.edges for w in g.successors(v) if w != u and g.has_edge(u, w)}): out.append(('FFL', f'{u} → {v} → {w_}; {u} → {w_}', [(u, v), (v, w_), (u, w_)]))
    rows = []
    for ty, name, es in out:
        nodes = sorted({x for p in es for x in p}); row = dict(network=label, type=ty, motif=name, motif_score=int(sum(tt.loc[p, 'score_0_9'] for p in es)), n_edges_direct=int(sum(tt.loc[p, 'Direct'] for p in es)), n_edges=len(es), contains_ESR1_or_FOXA1=bool(set(nodes) & {'ESR1', 'FOXA1'}), contexts_all_nodes=';'.join(k for k, v in CTX.items() if set(nodes) <= v) or 'crosses contexts')
        row['all_nodes_luminal_members'] = set(nodes) <= CTX['Luminal']
        for l in LAY: row[l] = float(np.mean([tt.loc[p, l] for p in es]))
        rows.append(row)
    return pd.DataFrame(rows)
mn = motifs(G, 'network'); print(mn[['type', 'motif', 'motif_score', 'n_edges_direct', 'contains_ESR1_or_FOXA1', 'all_nodes_luminal_members', 'contexts_all_nodes']].to_string())
G6 = nx.DiGraph(); [G6.add_edge(r.regulator_TF, r.target_TF) for r in tab[tab.score_0_9 >= 6].itertuples()]
m6 = motifs(G6, 'all_edges_ge6'); print('ge6: ', m6.type.value_counts().to_dict(), 'without ESR1/FOXA1:'); print(m6[~m6.contains_ESR1_or_FOXA1][['type', 'motif', 'motif_score', 'n_edges_direct', 'n_edges', 'contexts_all_nodes']].to_string())
pd.concat([mn, m6]).to_csv(O/'A19_motifs.tsv', sep='\t', index=False)
prs.rename('PageRank').to_csv(O/'A19_pagerank.tsv', sep='\t')
rep = []
for tf, g in extra.groupby('TF'):
    if tf not in T: continue
    oldocc = {t for t in T if t != tf and base.loc[tf, t] > 0}
    if not oldocc: continue
    anyp = {t for t in T if t != tf and float(g[t].astype(float).sum()) > 0}
    rep.append(dict(TF=tf, new_profiles=len(g), new_series=len(set().union(*(serials(x) for x in g.series))), earlier_occupied_promoters=len(oldocc), reproduced=len(oldocc & anyp), not_reproduced=';'.join(sorted(oldocc - anyp)), additional=';'.join(sorted(anyp - oldocc)), median_peaks=float(g.n_peaks.astype(float).median())))
R = pd.DataFrame(rep); print(R.to_string()); R.to_csv(O/'A19_replication_in_new_studies.tsv', sep='\t', index=False)
print('total', R.reproduced.sum(), '/', R.earlier_occupied_promoters.sum())

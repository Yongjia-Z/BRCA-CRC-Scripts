import sys, gzip, re
from pathlib import Path
import numpy as np, pandas as pd, networkx as nx
W, O = Path(sys.argv[1]), Path(sys.argv[2]); R = W/'core_data/remap_reanalysis'
g = pd.read_csv(W/'datasets/GENCODE/results/gencode_v36_gene_annotation.tsv', sep='\t').set_index('gene_name')
r = g.loc['PGR']; t = int(r.start if r.strand == '+' else r.end); CHR, A, B = r.chromosome, t-2000, t+2000
def hits(path):
    op = gzip.open if str(path).endswith('.gz') else open
    n = 0
    with op(path, 'rt') as h:
        for l in h:
            if l.startswith(('#', 'track', 'browser')): continue
            p = l.split('\t', 3)
            if p[0] == CHR and int(p[1]) < B and int(p[2]) > A: n += 1
    return n
st = pd.read_csv(R/'results/remap_breast_tier1_download_status.tsv', sep='\t'); st = st[st.local_peak_file_present].copy()
st['pgr_peaks'] = [hits(f) for f in st.local_peak_file]
st[['target_name','biotype','biotype_modification','experiment','n_peaks_reported','pgr_peaks']].to_csv(O/'A14_PGR_target_by_dataset.tsv', sep='\t', index=False)
b = st.groupby('target_name').agg(n_datasets=('pgr_peaks','size'), datasets_with_peak=('pgr_peaks', lambda x: int((x>0).sum())), peaks=('pgr_peaks','sum'))
tn = pd.read_csv(R/'results/TNBC_direct_binding/tnbc_direct_binding_dataset_manifest.tsv', sep='\t')
tn = tn[tn.usable_for_hg38_overlap.astype(str).str.lower().eq('true')] if 'usable_for_hg38_overlap' in tn else tn
tn['pgr_peaks'] = [hits(f) for f in tn.local_peak_file]
print(tn[['regulator_TF','cell_line','n_peaks_in_file','pgr_peaks']].to_string())
tt = tn.groupby('regulator_TF').pgr_peaks.sum()
b['tnbc_peaks'] = tt.reindex(b.index).fillna(0).astype(int)
for k in tt.index.difference(b.index): b.loc[k, ['n_datasets','datasets_with_peak','peaks','tnbc_peaks']] = [0,0,0,int(tt[k])]
b['matrix_value'] = np.log1p(b.peaks + 8*b.tnbc_peaks)
b.to_csv(O/'A14_PGR_target_summary.tsv', sep='\t'); print(b.to_string())
rx = re.compile(r'(^|[_-])si[A-Za-z]|^SI[A-Z]', re.I)
m = st.biotype_modification.fillna('').astype(str)
flag = st[m.map(lambda x: bool(rx.search(x)) and not re.search('siC(on|tr)', x, re.I))]
print('tier-1 knockdown-labelled:'); print(flag[['target_name','biotype','biotype_modification','experiment']].to_string())
F = W/'data/submitted_network'; a = pd.read_csv(F/'figure2a_all_candidate_edges_score_ge6.tsv', sep='\t'); e = pd.read_csv(F/'figure2a_edges.tsv', sep='\t')
def motifs(d):
    G = nx.DiGraph(); G.add_edges_from(zip(d.regulator_TF, d.target_TF))
    pairs = sorted({tuple(sorted(x)) for x in G.edges if G.has_edge(x[1], x[0])})
    cyc = ['>'.join(c) for c in nx.simple_cycles(G) if len(c) == 3]
    ffl = [f'{u}>{v}>{w}' for u, v in G.edges for w in G.successors(v) if w != u and G.has_edge(u, w)]
    return G, pairs, cyc, ffl
w = lambda d: dict(zip(zip(d.regulator_TF, d.target_TF), d.context_edge_score))
e2 = pd.concat([e, a[(a.regulator_TF=='ESR1')&(a.target_TF=='PGR')].assign(context_edge_score=7.0, context_direct_binding=True)], ignore_index=True)
a2 = a.copy(); a2.loc[(a2.regulator_TF=='ESR1')&(a2.target_TF=='PGR'), ['context_edge_score','context_direct_binding']] = [7.0, True]
out = []
for name, d in [('submitted rule, 23 edges', e), ('submitted rule + ESR1->PGR, 24 edges', e2), ('all score>=6, 35 edges', a2)]:
    G, p, c, f = motifs(d)
    H = nx.DiGraph(); [H.add_edge(x, y, weight=s) for (x, y), s in w(d).items()]
    pr = nx.pagerank(H, weight='weight')
    out.append(dict(network=name, edges=G.number_of_edges(), reciprocal_pairs=';'.join('<->'.join(x) for x in p), cycles=';'.join(c), n_ffl=len(f), ffl=';'.join(f),
                    pagerank_top3=';'.join(f'{k}={v:.3f}' for k, v in sorted(pr.items(), key=lambda x: -x[1])[:3]), FOXA1_out=G.out_degree('FOXA1')))
pd.DataFrame(out).to_csv(O/'A14_network_motifs.tsv', sep='\t', index=False); print(pd.DataFrame(out).drop(columns='ffl').T.to_string())
a2.assign(in_23=[(x, y) in w(e) for x, y in zip(a2.regulator_TF, a2.target_TF)])[['regulator_TF','target_TF','context_edge_score','context_direct_binding','in_23','figure2a_selection_reason']].to_csv(O/'A14_all_edges_ge6.tsv', sep='\t', index=False)
E = pd.read_csv(W/'core_data/figure2/raw/figure2b__TCGA_results__TF_TF_regulatory_edge_evidence_score_table_with_TRRUST.tsv', sep='\t')
T = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']
x = E[E.regulator_TF.isin(T) & E.target_TF.isin(T) & (E.regulator_TF != E.target_TF)]
print('scored ordered pairs', len(x), 'of 156; regulators', x.regulator_TF.value_counts().to_dict(), 'never regulators', sorted(set(T)-set(x.regulator_TF)))
print('all regulators in table', E.regulator_TF.nunique(), 'min RP among rows', E.max_BRCA_regular_potential.min())

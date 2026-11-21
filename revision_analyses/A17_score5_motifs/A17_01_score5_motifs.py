import sys, itertools
from pathlib import Path
import pandas as pd
A, OUT = Path(sys.argv[1]), Path(sys.argv[2])
E = pd.read_csv(A/'A16_network_edge_sensitivity.tsv', sep='\t', encoding='utf-8-sig')
CTX = {'Luminal': {'ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C'},
       'HER2': {'FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1'},
       'Basal/TNBC': {'RUNX1','TP63','FOXC1'}}
CORE = {'ESR1','FOXA1'}
rows = []; detail = []
for sc in ['S1_curated_baseline', 'S2_new_studies_unsorted']:
    e = E[E.scenario == sc]
    for thr in [6, 5]:
        ed = e[e.score_0_9 >= thr]; S = set(zip(ed.regulator_TF, ed.target_TF)); B = dict(zip(zip(e.regulator_TF, e.target_TF), e.binding))
        nodes = sorted({x for p in S for x in p})
        mot = []
        for a, b in itertools.combinations(nodes, 2):
            if (a, b) in S and (b, a) in S: mot.append(('reciprocal', (a, b), [(a, b), (b, a)]))
        for a, b, c in itertools.permutations(nodes, 3):
            if a == min(a, b, c) and (a, b) in S and (b, c) in S and (c, a) in S: mot.append(('cycle', (a, b, c), [(a, b), (b, c), (c, a)]))
        for t, n, es in mot:
            ctx = [k for k, v in CTX.items() if set(n) <= v]
            detail.append(dict(scenario=sc, min_score=thr, type=t, nodes='→'.join(n), n_edges_binding=sum(bool(B[x]) for x in es), n_edges=len(es),
                               all_edges_binding=all(bool(B[x]) for x in es), scores=';'.join(str(int(e[(e.regulator_TF == x[0]) & (e.target_TF == x[1])].score_0_9.iloc[0])) for x in es),
                               contains_ESR1_or_FOXA1=bool(set(n) & CORE), contexts_containing_all_nodes=';'.join(ctx) or 'none (crosses contexts)'))
        d = [m for m in mot]
        def cnt(c, onlybind=False):
            return sum(1 for t, n, es in mot if set(n) <= CTX[c] and (not onlybind or all(bool(B[x]) for x in es)))
        rows.append(dict(scenario=sc, min_score=thr, n_edges=len(S), reciprocal=sum(t == 'reciprocal' for t, _, _ in mot), cycles=sum(t == 'cycle' for t, _, _ in mot),
                         **{f'closed_in_{c}': cnt(c) for c in CTX}, **{f'closed_with_all_binding_in_{c}': cnt(c, True) for c in CTX}))
S = pd.DataFrame(rows); D = pd.DataFrame(detail)
S.to_csv(OUT/'A17_score5_summary.tsv', sep='\t', index=False); D.to_csv(OUT/'A17_score5_motifs.tsv', sep='\t', index=False)
print(S.to_string())
for sc, thr in [('S1_curated_baseline', 5), ('S2_new_studies_unsorted', 5), ('S2_new_studies_unsorted', 6)]:
    d = D[(D.scenario == sc) & (D.min_score == thr)]
    print('\n', sc, thr, len(d)); print(d[d.contexts_containing_all_nodes.str.contains('HER2|Basal')][['type', 'nodes', 'scores', 'n_edges_binding', 'all_edges_binding', 'contains_ESR1_or_FOXA1', 'contexts_containing_all_nodes']].to_string())
e = E[E.scenario == 'S2_new_studies_unsorted']; print(e.score_0_9.value_counts().sort_index()); print(e[e.score_0_9 == 5][['regulator_TF', 'target_TF', 'binding']].to_string())
print(e[e.regulator_TF.isin(['RUNX1', 'TP63', 'FOXC1']) | e.target_TF.isin(['TP63', 'FOXC1'])][['regulator_TF', 'target_TF', 'binding', 'score_0_9', 'in_network']].to_string())

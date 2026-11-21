import sys, re
from pathlib import Path
import pandas as pd, networkx as nx
A16, A19, O = map(Path, sys.argv[1:4]); O.mkdir(parents=True, exist_ok=True)
T = ['ESR1', 'GATA3', 'PGR', 'PBX1', 'FOXA1', 'SPDEF', 'ELF3', 'EHF', 'AR', 'TFAP2C', 'RUNX1', 'TP63', 'FOXC1']
LINE_CTX = {'MCF-7': 'Luminal', 'T-47D': 'Luminal', 'ZR-75-1': 'Luminal', 'BT-474': 'HER2', 'SK-BR-3': 'HER2', 'HCC1954': 'HER2',
            'MDA-MB-231': 'Basal/TNBC', 'MDA-MB-468': 'Basal/TNBC', 'BT-549': 'Basal/TNBC', 'BT549': 'Basal/TNBC', 'Hs578T': 'Basal/TNBC',
            'HCC70': 'Basal/TNBC', 'HCC38': 'Basal/TNBC', 'SUM159PT': 'Basal/TNBC', 'HCC1806': 'Basal/TNBC', 'HCC1395': 'Basal/TNBC'}
MEMBERS = {'HER2': {'FOXA1', 'SPDEF', 'ELF3', 'EHF', 'AR', 'TFAP2C', 'RUNX1'}, 'Basal/TNBC': {'RUNX1', 'TP63', 'FOXC1'}}
d = pd.read_csv(A16/'results/A16_occupancy_by_dataset.tsv', sep='\t', encoding='utf-8-sig').fillna('')
d['use_strict'] = d.use_strict.astype(str).str.lower().eq('true')
prior = d[d.use_strict & d.source.isin(['ReMap baseline', 'TNBC baseline author/ReMap', 'PGR curated ChIP-Atlas'])]
new = pd.read_csv(A19/'A19_S2_new_study_profiles_used.tsv', sep='\t').fillna('')
fr = pd.concat([prior, new[prior.columns.intersection(new.columns)]], ignore_index=True)
fr['line'] = fr.cell_line.astype(str).str.strip()
unknown = sorted(set(fr.line) - set(LINE_CTX)); assert not unknown, unknown
fr['line_context'] = fr.line.map(LINE_CTX)
print('profiles', len(fr), fr.groupby(['TF', 'line']).size().unstack(fill_value=0).to_string())
e = pd.read_csv(A19/'A19_edges.tsv', sep='\t')
net = e[e.in_network.astype(str).str.lower().eq('true')]
rows = []
for r in net.itertuples():
    sub = fr[fr.TF.eq(r.regulator_TF) & (fr[r.target_TF].astype(float) > 0)]
    lines = sorted(set(sub.line)); ctxs = sorted(set(sub.line_context))
    rows.append(dict(regulator=r.regulator_TF, target=r.target_TF, score=int(r.score_0_9), binding_in_score=int(r.Direct),
                     n_profiles_with_peak=len(sub), n_profiles_assayed=int(fr.TF.eq(r.regulator_TF).sum()),
                     cell_lines_with_peak=';'.join(lines), contexts_with_peak=';'.join(ctxs),
                     HER2_line_peak=int('HER2' in ctxs), TNBC_line_peak=int('Basal/TNBC' in ctxs)))
E = pd.DataFrame(rows).sort_values(['score', 'regulator'], ascending=[False, True]); E.to_csv(O/'A22_1_network_edges_by_cell_line.tsv', sep='\t', index=False)
print(E.to_string())
def motifs(edges):
    g = nx.DiGraph(); g.add_edges_from(edges)
    pairs = sorted({tuple(sorted(x)) for x in g.edges if g.has_edge(x[1], x[0])})
    cyc = sorted({tuple(c[c.index(min(c)):] + c[:c.index(min(c))]) for c in nx.simple_cycles(g) if len(c) == 3})
    return pairs, cyc
summ = []
for line, ctx in sorted({(l, LINE_CTX[l]) for l in fr.line}):
    eds = [(x.regulator, x.target) for x in E.itertuples() if line in x.cell_lines_with_peak.split(';')]
    pairs, cyc = motifs(eds)
    regs = sorted(set(fr.TF[fr.line.eq(line)]))
    inside = {k: [p for p in pairs if set(p) <= v] for k, v in MEMBERS.items()}
    summ.append(dict(cell_line=line, context=ctx, regulators_profiled=';'.join(regs), network_edges_with_peak=len(eds),
                     reciprocal_pairs=';'.join('↔'.join(p) for p in pairs), three_node_cycles=';'.join('→'.join(c) for c in cyc),
                     HER2_member_pairs=';'.join('↔'.join(p) for p in inside['HER2']), BasalTNBC_member_pairs=';'.join('↔'.join(p) for p in inside['Basal/TNBC'])))
Sm = pd.DataFrame(summ).sort_values(['context', 'network_edges_with_peak'], ascending=[True, False]); Sm.to_csv(O/'A22_2_motifs_within_cell_line.tsv', sep='\t', index=False)
print(Sm.to_string())
print('edges with binding:', int(E.binding_in_score.sum()), 'of', len(E), '| with a peak in a HER2 line:', int(E.HER2_line_peak.sum()), '| in a TNBC line:', int(E.TNBC_line_peak.sum()))

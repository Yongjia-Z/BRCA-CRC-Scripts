import sys, re
from pathlib import Path
import pandas as pd
mot = pd.read_csv(sys.argv[1], encoding='utf-8-sig'); OUT = Path(sys.argv[3])
calls = pd.read_csv(sys.argv[2], sep='\t', header=[0, 1], index_col=0)['corrected']
mem = {tf: set(str(calls.loc[tf, 'active']).split(';')) for tf in calls.index}
rows = []
for _, m in mot.iterrows():
    nodes = sorted({n for e in m[['edge_1', 'edge_2', 'edge_3']].dropna() for n in e.split('->')})
    ctx = {n: '/'.join(sorted(mem[n])) for n in nodes}
    all_lum = all('Luminal' in mem[n] for n in nodes)
    all_her2 = all('HER2' in mem[n] for n in nodes)
    non_lum = [n for n in nodes if 'Luminal' not in mem[n]]
    rows.append({'motif': m.motif, 'type': m.type, 'closed_feedback': m.type in ('reciprocal_pair', 'three_node_cycle'),
                 'motif_score': m.motif_score, 'nodes': ';'.join(nodes), 'node_contexts': '; '.join(f'{n}:{c}' for n, c in ctx.items()),
                 'contains_ESR1_or_FOXA1': bool({'ESR1', 'FOXA1'} & set(nodes)), 'all_nodes_Luminal_member': all_lum,
                 'all_nodes_HER2_member': all_her2, 'nodes_without_Luminal_membership': ';'.join(non_lum),
                 'class': ('cross-context (includes non-Luminal node)' if non_lum else
                           'Luminal/HER2 shared members only' if all_her2 else 'Luminal-anchored')})
R = pd.DataFrame(rows); R.to_csv(OUT/'A13_motif_node_membership.tsv', sep='\t', index=False)
S = R.groupby(['type', 'class']).size().rename('n').reset_index(); S.to_csv(OUT/'A13_motif_class_summary.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(R[['motif', 'type', 'node_contexts', 'class']].to_string(index=False)); print(S.to_string(index=False))
print('all contain ESR1/FOXA1:', R.contains_ESR1_or_FOXA1.all(), '| closed feedback motifs all Luminal-member:', R[R.closed_feedback].all_nodes_Luminal_member.all())

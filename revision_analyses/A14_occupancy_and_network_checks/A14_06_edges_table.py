import sys
from pathlib import Path
import pandas as pd
O = Path(sys.argv[1]); a = pd.read_csv(O/'A14_all_edges_ge6.tsv', sep='\t')
a['in_24'] = a.in_23 | ((a.regulator_TF == 'ESR1') & (a.target_TF == 'PGR'))
a = a[['regulator_TF', 'target_TF', 'context_edge_score', 'context_direct_binding', 'in_23', 'in_24']]
a.columns = ['regulator_TF', 'target_TF', 'edge_score', 'direct_binding', 'in_original_23_edges', 'in_network_24_edges']
a.to_csv(O/'A14_edges_ge6_for_S15b.tsv', sep='\t', index=False); print(int(a.in_network_24_edges.sum()), len(a))

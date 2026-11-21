import sys, re
from pathlib import Path
import numpy as np, pandas as pd
W, O = Path(sys.argv[1]), Path(sys.argv[2]); R = W/'core_data/remap_reanalysis/results'
d = pd.read_csv(R/'remap_breast_tier1_promoter_overlap_by_dataset.tsv', sep='\t')
st = pd.read_csv(R/'remap_breast_tier1_download_status.tsv', sep='\t')
kd = st[st.biotype_modification.fillna('').astype(str).isin(['SIGATA_E2', 'siFEN1', 'siFOXA1-IL6', 'siFOXA1-Veh'])]
print(kd[['target_name','biotype','biotype_modification','experiment','local_peak_file_present']].to_string())
d2 = d[~d.dataset_file.isin(kd.dataset_file)]
def tab(x):
    return x.groupby(['regulator_TF','target_TF']).agg(ds=('dataset_file','nunique'), peaks=('n_peak_overlaps','sum'))
a, b = tab(d), tab(d2)
c = a.join(b, rsuffix='_excl', how='outer').fillna(0).astype(int); c['changed'] = (c.ds != c.ds_excl) | (c.peaks != c.peaks_excl)
c.to_csv(O/'A14_knockdown_label_effect.tsv', sep='\t')
core = ['ESR1','FOXA1','GATA3']
print(c.loc[[(x, y) for x in core for y in core if (x, y) in c.index]].to_string())
print('pairs lost entirely', c[(c.ds > 0) & (c.ds_excl == 0)].index.tolist())
n0 = a.groupby('regulator_TF').ds.sum(); n1 = b.groupby('regulator_TF').ds.sum()
s = pd.DataFrame({'n': n0, 'n_excl': n1}); s['score'] = np.log1p(s.n)/np.log1p(s.n.max()); s['score_excl'] = np.log1p(s.n_excl)/np.log1p(s.n_excl.max())
print(s.round(3).to_string())

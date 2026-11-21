import sys
from pathlib import Path
import numpy as np, pandas as pd
f = pd.read_csv(sys.argv[1]).fillna({'active_contexts': ''}); OUT = Path(sys.argv[2])
SEED = 20261026; rng = np.random.default_rng(SEED)
cols = ['RNA','TF_activity_regulon','ATAC_motif','ATAC_footprint','ATAC_peak_to_gene_accessibility','context_matched_H3K27ac_SE_MAX','context_matched_direct_binding_MAX']
mat = f[cols].fillna(0).to_numpy(float); s = f['max_state_score_0_5'].to_numpy(); spread = f['state_spread_max_minus_second'].to_numpy(); sup = f['Supportive_score_0_3'].to_numpy()
def select(core, offset=0):
    high = (core >= 5.5+offset) & (s >= 2.6)
    prov = ((core >= 4.75+offset) & (s >= 3) & (sup >= .5)) | ((core >= 4.75+offset) & (s >= 3.5) & (spread >= 1)) | ((core+sup >= 6+offset) & (s >= 2.4))
    active = (s >= 2.8) | ((core >= 5.5+offset) & (s >= 2.6)) | ((core >= 4.5+offset) & (s >= 2.35)) | (s >= 2.35)
    return (high | prov) & active
base = select(mat.sum(1)); formal = f.loc[base, 'TF'].tolist()
assert set(formal) == {'ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1'}, formal
sens = []
for j, col in enumerate(cols):
    hit = select((mat.sum(1)-mat[:, j])*7/6)
    sens.append({'analysis':'leave_one_discovery_layer_out_rescaled','omitted_layer':col,'n_selected':int(hit.sum()),'original_retained':int((hit&base).sum()),'lost':';'.join(f.loc[base&~hit,'TF']),'new':';'.join(f.loc[~base&hit,'TF'])})
for off in [-.5, -.25, 0, .25, .5]:
    hit = select(mat.sum(1), off)
    sens.append({'analysis':'core_total_threshold_shift','omitted_layer':str(off),'n_selected':int(hit.sum()),'original_retained':int((hit&base).sum()),'lost':';'.join(f.loc[base&~hit,'TF']),'new':';'.join(f.loc[~base&hit,'TF'])})
pd.DataFrame(sens).to_csv(OUT/'A8_sensitivity_summary.tsv', sep='\t', index=False)
hits = []
for i in range(2000):
    w = rng.uniform(.8, 1.2, 7); w = w*7/w.sum(); hits.append(select(mat@w))
freq = np.mean(hits, axis=0)
ft = f[['TF','confidence_tier','Core_discovery_score_v2_0_7','max_state_score_0_5','active_contexts']].copy(); ft['selection_frequency_2000_weight_perturbations'] = freq
ft.to_csv(OUT/'A8_weight_sensitivity_all_TFs.tsv', sep='\t', index=False)
null = []
for i in range(2000):
    ix = rng.choice(len(f), 13, replace=False)
    null.append({'replicate':i+1,'mean_core_score':float(mat[ix].sum(1).mean()),'n_passing_original_rules':int(base[ix].sum())})
null = pd.DataFrame(null); null.to_csv(OUT/'A8_selection_stage_random_sets.tsv', sep='\t', index=False)
print(pd.DataFrame(sens).to_string())
print('weight freq formal:', dict(zip(f.loc[base,'TF'], np.round(freq[base],3))))
print('nonformal >0.05:', ft[~base & (freq > .05)][['TF','selection_frequency_2000_weight_perturbations']].values.tolist())
print('random: observed mean core', round(mat[base].sum(1).mean(),3), 'median', round(null.mean_core_score.median(),3), 'max', round(null.mean_core_score.max(),3), 'max passing', null.n_passing_original_rules.max())
print(f[f.TF.isin(['XBP1','ZNF703'])][['TF','Core_discovery_score_v2_0_7','max_state_score_0_5']])

import sys, glob
from pathlib import Path
import numpy as np, pandas as pd
R, G, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); OUT.mkdir(parents=True, exist_ok=True)
D = R/'datasets'
a = pd.read_csv(R/'source_tables/figure1/figure1a_all_tf_scores.csv').set_index('TF')
U = a.index
g = pd.read_csv(D/'GENCODE/results/gencode_v36_gene_annotation.tsv', sep='\t').rename(columns={'gene_name':'TF','chromosome':'chrom'})
g = g[g.TF.isin(U)].drop_duplicates('TF')
g['rs'] = (g.start-100000).clip(lower=0); g['re'] = g.end+100000
def counts(files):
    res = {}
    for f in files:
        p = pd.read_csv(f, sep='\t', header=None, usecols=[0,1,2], comment='#')
        n = {}
        for chrom, sub in g.groupby('chrom'):
            pk = p[p[0]==chrom]
            for tf, rs, re_ in zip(sub.TF, sub.rs, sub.re): n[tf] = int(((pk[1]<re_)&(pk[2]>rs)).sum())
        res[Path(f).name] = pd.Series(n)
    return pd.DataFrame(res).reindex(U).fillna(0)
def ls01(x):
    x = np.log1p(pd.Series(x, index=U).astype(float).fillna(0)); return (x-x.min())/(x.max()-x.min())
h = pd.read_csv(D/'GSE85158_H3K27ac/results/GSE85158_TF_window_H3K27ac_subtype_summary_hg19_to_hg38.tsv', sep='\t').drop_duplicates('TF').set_index('TF')
db = pd.read_csv(D/'SEdb/results/dbSUPER_breast_CRC_TF_100kb_overlap_all.tsv', sep='\t')
dbh = db[db.cell_line.eq('HCC1954')].groupby('TF').n_dbSUPER_SE_100kb.max().reindex(U).fillna(0)
her2c = counts(sorted(glob.glob(str(D/'GSE251868_HER2_H3K27ac/data/*.bed.gz'))))
her2 = np.maximum.reduce([ls01(her2c.max(1)), ls01(dbh), ls01(h['GSE85158_HER2_minus_normal_mean_signal'].reindex(U).fillna(0))])
tn = counts(sorted(glob.glob(str(G/'*.narrowPeak.gz'))))
tn.loc[['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']].to_csv(OUT/'A7_GSE245830_peak_counts_CRC13.tsv', sep='\t')
gb = ls01(h['GSE85158_Basal_TNBC_minus_normal_mean_signal'].reindex(U).fillna(0))
EV = [c for c in tn.columns if '_EV_' in c]; assert len(EV) == 3
basal_corr = np.maximum(ls01(tn[EV].max(1)), gb)
val = pd.DataFrame({'check': ['HER2 layer reproduced (Pearson r)', 'stored Basal layer vs GSE85158-only (r)', 'stored Basal layer vs GSE85158+GSE245830 (r)',
                              'max abs diff HER2 reproduced', 'GSE245830 files found by figure1a.py glob'],
                    'value': [np.corrcoef(a.HER2_H3K27ac_SE, her2)[0,1], np.corrcoef(a.Basal_TNBC_H3K27ac_SE, gb)[0,1], np.corrcoef(a.Basal_TNBC_H3K27ac_SE, basal_corr)[0,1],
                              float(np.abs(a.HER2_H3K27ac_SE-her2).max()), len(glob.glob(str(D/'GSE245830_TNBC_H3K27ac/data/*.bed*')))]})
val.to_csv(OUT/'A7_validation.tsv', sep='\t', index=False); print(val)
b = a.copy()
b['Basal_TNBC_H3K27ac_SE'] = basal_corr
newmax = b[['Luminal_H3K27ac_SE','HER2_H3K27ac_SE','Basal_TNBC_H3K27ac_SE']].max(axis=1)
b['Core_discovery_score_v2_0_7'] = a.Core_discovery_score_v2_0_7 - a.context_matched_H3K27ac_SE_MAX + newmax
b['context_matched_H3K27ac_SE_MAX'] = newmax
b['Audited_total_score_0_10'] = b.Core_discovery_score_v2_0_7 + a.Supportive_score_0_3
b['state_score_Basal/TNBC'] = a['state_score_Basal/TNBC'] + (basal_corr - a.Basal_TNBC_H3K27ac_SE)
CTX = ['Luminal','HER2','Basal/TNBC']
def call(r):
    s = {c: float(r['state_score_'+c]) for c in CTX}; mx = max(s.values()); sec = sorted(s.values())[-2]
    core, tot = r.Core_discovery_score_v2_0_7, r.Audited_total_score_0_10
    P = [c for c, v in s.items() if v >= mx-0.20 and v >= 2.80]
    if not P and core >= 5.50 and mx >= 2.60: P = [max(s, key=s.get)]
    if not P and core >= 4.50 and mx >= 2.35: P = [max(s, key=s.get)]
    S = [c for c, v in s.items() if c not in P and v >= 2.35 and v >= 0.82*mx]
    if core >= 5.50 and mx >= 2.60: t = 'High-confidence'
    elif (core >= 4.75 and mx >= 3.00 and tot-core >= 0.50) or (core >= 4.75 and mx >= 3.50 and mx-sec >= 1.00) or (tot >= 6.00 and mx >= 2.40): t = 'Provisional'
    else: t = 'Supporting'
    st = {c: ('P' if c in P else ('S' if c in S else '')) for c in CTX}
    return pd.Series({'tier': t, 'primary': ';'.join(P), 'active': ';'.join(P+S), **{c+'_state': st[c] for c in CTX}, 'formal': t != 'Supporting' and bool(P+S)})
orig = a.apply(call, axis=1); corr = b.apply(call, axis=1)
assert (orig.tier == a.confidence_tier).all() and (orig.active.fillna('') == a.active_contexts.fillna('')).all(), 'rule re-implementation does not reproduce stored calls'
fo, fc = set(orig.index[orig.formal]), set(corr.index[corr.formal])
print('formal stored', len(fo), 'corrected', len(fc), 'added', sorted(fc-fo), 'lost', sorted(fo-fc))
cols = ['Basal_TNBC_H3K27ac_SE','Core_discovery_score_v2_0_7','state_score_Luminal','state_score_HER2','state_score_Basal/TNBC']
rows = sorted(fo | fc, key=lambda t: list(a.index).index(t))
tab = pd.concat({'stored': a.loc[rows, cols].join(orig.loc[rows]), 'corrected': b.loc[rows, cols].join(corr.loc[rows])}, axis=1).round(3)
tab.to_csv(OUT/'A7_CRC_TF_calls_stored_vs_corrected.tsv', sep='\t'); print(tab.to_string())
chg = [t for t in rows if (orig.loc[t, ['tier','active','Luminal_state','HER2_state','Basal/TNBC_state']] != corr.loc[t, ['tier','active','Luminal_state','HER2_state','Basal/TNBC_state']]).any()]
print('TFs with any change in tier/context:', chg)
sc = b[['state_score_Luminal','state_score_HER2','state_score_Basal/TNBC']].to_numpy()
b['max_state_score_0_5'] = sc.max(1); b['second_state_score_0_5'] = np.sort(sc, 1)[:, -2]
b['state_spread_max_minus_second'] = b.max_state_score_0_5 - b.second_state_score_0_5
for c in CTX: b['state_strength_'+c] = b['state_score_'+c]/5.0
b['Core_discovery_score_0_7'] = b.Core_discovery_score_v2_0_7; b['H3K27ac_SE_MAX'] = b.context_matched_H3K27ac_SE_MAX
b['confidence_tier'] = corr.tier; b['primary_context'] = corr.primary; b['active_contexts'] = corr.active
for c in CTX: b[c+'_state'] = corr[c+'_state']; b[c+'_member'] = corr[c+'_state'].isin(['P','S']).astype(int)
b.to_csv(OUT/'A7_all_tf_scores_corrected_basal_H3K27ac.csv')

ball = np.maximum(ls01(tn.max(1)), gb)
sens = []
for lab, layer in [('EV control files only (main analysis)', basal_corr), ('all 6 GSE245830 files incl. shRNA (sensitivity)', ball)]:
    bb = a.copy(); bb['Basal_TNBC_H3K27ac_SE'] = layer
    nm = bb[['Luminal_H3K27ac_SE','HER2_H3K27ac_SE','Basal_TNBC_H3K27ac_SE']].max(axis=1)
    bb['Core_discovery_score_v2_0_7'] = a.Core_discovery_score_v2_0_7 - a.context_matched_H3K27ac_SE_MAX + nm
    bb['Audited_total_score_0_10'] = bb.Core_discovery_score_v2_0_7 + a.Supportive_score_0_3
    bb['state_score_Basal/TNBC'] = a['state_score_Basal/TNBC'] + (layer - a.Basal_TNBC_H3K27ac_SE)
    cc = bb.apply(call, axis=1); f = set(cc.index[cc.formal])
    sens.append({'variant': lab, 'n_formal': len(f), 'same_13': f == fo, 'same_as_main': f == fc, 'RUNX1_Basal_layer': round(float(layer['RUNX1']), 4),
                 'RUNX1_state_Luminal': round(float(bb.loc['RUNX1','state_score_Luminal']), 4), 'RUNX1_state_HER2': round(float(bb.loc['RUNX1','state_score_HER2']), 4),
                 'RUNX1_state_Basal': round(float(bb.loc['RUNX1','state_score_Basal/TNBC']), 4), 'RUNX1_active': cc.loc['RUNX1','active'],
                 'formal_TFs_with_tier_or_context_differing_from_main': ';'.join(t for t in fc if (cc.loc[t,['tier','active']] != corr.loc[t,['tier','active']]).any()),
                 'max_abs_core_score_diff_vs_main_13TF': round(float((bb.loc[sorted(fc),'Core_discovery_score_v2_0_7'] - b.loc[sorted(fc),'Core_discovery_score_v2_0_7']).abs().max()), 3)})
pd.DataFrame(sens).to_csv(OUT/'A7_sensitivity_GSE245830_samples.tsv', sep='\t', index=False); print(pd.DataFrame(sens).to_string())

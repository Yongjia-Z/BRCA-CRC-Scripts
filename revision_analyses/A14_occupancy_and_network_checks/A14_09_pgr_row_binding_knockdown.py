import sys
from pathlib import Path
import numpy as np, pandas as pd
W, X, O = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]).resolve(); R = W/'core_data/remap_reanalysis/results'
T = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']
C16 = ['ESR1','FOXA1','SPDEF','GATA3','PBX1','AR','CREB3L4','TFAP2C','EHF','MYBL2','NR2F6','STAT1','ELF3','RUNX1','TP63','FOXC1']
g = pd.read_csv(W/'datasets/GENCODE/results/gencode_v36_gene_annotation.tsv', sep='\t').set_index('gene_name')
win = {n: (g.loc[n, 'chromosome'], int(g.loc[n, 'start'] if g.loc[n, 'strand'] == '+' else g.loc[n, 'end'])) for n in dict.fromkeys(T + C16)}
byc = {}
for n, (c, t) in win.items(): byc.setdefault(c, []).append((n, t - 2000, t + 2000))
man = pd.read_csv(W/'data/dataset_lists/A14_chipatlas_PGR_manifest.tsv', sep='\t')
man = man[man.srx.isin(pd.read_csv(W/'data/dataset_lists/A14_PGR_curation.tsv', sep='\t').query('keep').srx)]
rows = []
for r in man.itertuples(index=False):
    cnt = dict.fromkeys(win, 0); tot = 0
    for l in open(O/'raw_chipatlas_PGR'/f'{r.srx}.05.bed'):
        if l.startswith(('track', '#', 'browser')): continue
        p = l.split('\t', 3); tot += 1
        for n, a, b in byc.get(p[0], []):
            if int(p[1]) < b and int(p[2]) > a: cnt[n] += 1
    rows.append(dict(srx=r.srx, cell_type=r.cell_type, title=r.title, n_peaks=tot, **cnt))
P = pd.DataFrame(rows); P.to_csv(O/'A14_chipatlas_PGR_by_dataset.tsv', sep='\t', index=False)
print('PGR profiles', len(P), P.cell_type.value_counts().to_dict(), 'median peaks', P.n_peaks.median(), 'empty', int((P.n_peaks == 0).sum()))
d = pd.read_csv(R/'remap_breast_tier1_promoter_overlap_by_dataset.tsv', sep='\t')
st = pd.read_csv(R/'remap_breast_tier1_download_status.tsv', sep='\t')
cu = pd.read_csv(W/'data/dataset_lists/A14_remap_curation.tsv', sep='\t'); kd = cu[~cu.keep].dataset_file
d2 = d[~d.dataset_file.isin(kd)]
m0 = pd.read_excel(X, sheet_name='S5_Promoter_occupancy').set_index('regulator_TF').loc[T, T].astype(float)
pk = d.groupby(['regulator_TF', 'target_TF']).n_peak_overlaps.sum(); pk2 = d2.groupby(['regulator_TF', 'target_TF']).n_peak_overlaps.sum()
m = m0.copy()
for (a, b), v in pk.items():
    if a in T and b in T:
        assert abs(np.log1p(v) - m0.loc[a, b]) < 1e-6 or m0.loc[a, b] > np.log1p(v), (a, b)
        tn = (np.expm1(m0.loc[a, b]) - v)/8
        m.loc[a, b] = np.log1p(pk2.get((a, b), 0) + 8*tn)
tg = pd.read_csv(O/'A14_PGR_target_by_dataset.tsv', sep='\t')
tg = tg.reset_index(drop=True); stl = st[st.local_peak_file_present].reset_index(drop=True); assert len(tg) == len(stl); tg = tg[~stl.dataset_file.isin(kd).values]
s = tg.groupby('target_name').agg(n_datasets=('pgr_peaks', 'size'), datasets_with_peak=('pgr_peaks', lambda x: int((x > 0).sum())), peaks=('pgr_peaks', 'sum'))
s.to_csv(O/'A14_PGR_target_summary_final.tsv', sep='\t'); print(s.to_string())
for tf in T:
    m.loc[tf, 'PGR'] = np.log1p(s.loc[tf, 'peaks']) if tf in s.index else 0.0
for tf in T: m.loc['PGR', tf] = np.log1p(P[tf].sum())
m.reset_index().to_csv(O/'A14_occupancy_matrix_baseline.tsv', sep='\t', index=False); print(m.round(2).to_string())
su = pd.read_csv(O/'A14_self_occupancy_support.tsv', sep='\t').set_index('TF')
u = st[st.local_peak_file_present & ~st.dataset_file.isin(kd)]
for tf in ['ESR1', 'FOXA1', 'GATA3', 'AR', 'TFAP2C', 'SPDEF']:
    x = d2[(d2.regulator_TF == tf) & (d2.target_TF == tf)]
    su.loc[tf, ['n_datasets', 'datasets_with_self_peak', 'self_peaks']] = [(u.target_name == tf).sum(), x.dataset_file.nunique(), x.n_peak_overlaps.sum()]
su.loc['PGR', 'source'] = 'ChIP-Atlas (curated MCF-7, T-47D and ZR-75-1 profiles); ReMap PGR datasets not available'
su.loc['PGR', ['n_datasets', 'datasets_with_self_peak', 'self_peaks', 'largest_dataset_peaks']] = [len(P), int((P.PGR > 0).sum()), int(P.PGR.sum()), int(P.n_peaks.max())]
su = su.drop(columns=['ChIPAtlas_crosscheck_datasets', 'ChIPAtlas_with_self_peak', 'ChIPAtlas_largest_peaks']).join(
    pd.read_csv(O/'A14_chipatlas_promoter_occupancy.tsv', sep='\t').query("regulator_TF != 'PGR'").groupby('regulator_TF').agg(ChIPAtlas_crosscheck_datasets=('srx', 'size'), ChIPAtlas_with_self_peak=('self_peaks', lambda x: int((x > 0).sum())), ChIPAtlas_largest_peaks=('n_peaks', 'max')))
su.reset_index().to_csv(O/'A14_self_occupancy_support.tsv', sep='\t', index=False); print(su.drop(columns='source').to_string())
n_pgr = int((P[C16] > 0).sum().sum())
n0 = d2.groupby('regulator_TF').apply(lambda x: x.groupby('target_TF').dataset_file.nunique().sum(), include_groups=False)
mx = max(n0.max(), n_pgr); sc = np.log1p(n_pgr)/np.log1p(mx)
A7 = pd.read_csv(O.parent/'A7_Basal_TNBC_H3K27ac_layer/A7_all_tf_scores_corrected_basal_H3K27ac.csv').set_index(pd.read_csv(O.parent/'A7_Basal_TNBC_H3K27ac_layer/A7_all_tf_scores_corrected_basal_H3K27ac.csv').columns[0])
core = A7.loc['PGR', 'Core_discovery_score_v2_0_7']; cols = [c for c in A7.columns if 'tate_score' in c or 'context' in c.lower()][:12]
print('n ESR1', n0.max(), 'n PGR', n_pgr, 'binding score', round(sc, 3), 'core', core, '->', round(core + sc, 3)); print(A7.loc['PGR', cols].to_string())
pd.DataFrame([dict(n_dataset_promoter_overlaps=n_pgr, max_n=int(mx), binding_score=round(sc, 3), core_submitted=core, core_with_binding=round(core + sc, 3), high_confidence_threshold=5.50, tier_changes=bool(core + sc >= 5.50))]).to_csv(O/'A14_PGR_binding_sensitivity.tsv', sep='\t', index=False)

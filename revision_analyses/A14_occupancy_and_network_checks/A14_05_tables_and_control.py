import sys
from pathlib import Path
import numpy as np, pandas as pd
W, X, O = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
T = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']
m = pd.read_excel(X, sheet_name='S5_Promoter_occupancy').set_index('regulator_TF').loc[T, T].astype(float)
s = pd.read_csv(O/'A14_PGR_target_summary.tsv', sep='\t', index_col=0)
for tf in T:
    if tf in s.index: m.loc[tf, 'PGR'] = s.loc[tf, 'matrix_value']
ca = pd.read_csv(O/'A14_chipatlas_promoter_occupancy.tsv', sep='\t')
pgr = ca[ca.srx.eq('SRX185709')].iloc[0]
for tf in T: m.loc['PGR', tf] = np.log1p(pgr[tf])
m.reset_index().to_csv(O/'A14_occupancy_matrix_baseline.tsv', sep='\t', index=False)
R = W/'core_data/remap_reanalysis/results'
st = pd.read_csv(R/'remap_breast_tier1_download_status.tsv', sep='\t'); st = st[st.local_peak_file_present]
e = pd.read_csv(R/'remap_breast_tier1_promoter_overlap_edges.tsv', sep='\t'); e = e[e.regulator_TF.eq(e.target_TF)].set_index('regulator_TF')
d = pd.read_csv(R/'TNBC_direct_binding/tnbc_direct_binding_overlap_by_dataset.tsv', sep='\t')
rows = []
for tf in T:
    a = st[st.target_name.eq(tf)]; c = ca[ca.regulator_TF.eq(tf)]
    n, hit, pk = len(a), int(e.loc[tf, 'n_breast_remap_dataset_overlaps']) if tf in e.index else 0, int(e.loc[tf, 'n_breast_remap_peak_overlaps']) if tf in e.index else 0
    src, mx = 'ReMap 2022 (breast tier 1)', int(a.n_peaks_reported.max()) if n else 0
    if tf == 'FOXC1':
        f = d[d.regulator_TF.eq('FOXC1') & d.target_TF.eq('FOXC1') & d.window_type.eq('TSS_pm_2kb_promoter')]
        n, hit, pk, src, mx = 4, f.dataset.size, int(f.n_peak_overlaps.sum()), 'GSE213840 (TNBC ChIP-seq)', 14311
    if tf == 'PGR':
        n, hit, pk, src, mx = 1, int(pgr.self_peaks > 0), int(pgr.self_peaks), 'ChIP-Atlas SRX185709 (T-47D); ReMap PGR datasets not acquired', int(pgr.n_peaks)
    best = max(mx, int(c.n_peaks.max()) if len(c) else 0)
    status = 'self-locus peak' if hit else ('no eligible dataset' if best == 0 else 'unresolved: low-yield data only' if best < 3000 else 'no self-locus peak')
    rows.append(dict(TF=tf, source=src, n_datasets=n, datasets_with_self_peak=hit, self_peaks=pk, largest_dataset_peaks=mx,
                     ChIPAtlas_crosscheck_datasets=len(c), ChIPAtlas_with_self_peak=int((c.self_peaks > 0).sum()), ChIPAtlas_largest_peaks=int(c.n_peaks.max()) if len(c) else 0, status=status))
S = pd.DataFrame(rows); S.to_csv(O/'A14_self_occupancy_support.tsv', sep='\t', index=False); print(S.drop(columns='source').to_string())

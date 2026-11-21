import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
W, C, O = map(Path, sys.argv[1:4]); O.mkdir(parents=True, exist_ok=True)
SEED = 20261102; REPS = 200
F = ['ESR1', 'GATA3', 'PGR', 'PBX1', 'FOXA1', 'SPDEF', 'ELF3', 'EHF', 'AR', 'TFAP2C', 'RUNX1', 'TP63', 'FOXC1']
INTER = W/'source_tables/intermediate/figure1a'
X = pd.read_csv(INTER/'figure1a_tcga_tf_expression_matrix.csv.gz', compression='gzip', index_col=0)
arch = pd.read_csv(INTER/'figure1a_tcga_expression_evidence.csv').set_index('TF')
cli = pd.read_csv(W/'datasets/TCGA/data/cbioportal_brca_patient_clinical_wide.tsv', sep='\t').set_index('patientId').SUBTYPE
prim = [c for c in X.columns if c[13:15] == '01']
lab = pd.Series({c: cli.get(c[:12]) for c in prim}).dropna()
ORIG = {'BRCA_LumA': 'LumA', 'BRCA_LumB': 'LumB', 'BRCA_Her2': 'HER2', 'BRCA_Basal': 'Basal'}
g_orig = lab[lab.isin(ORIG)].map(ORIG)
g_five = lab[lab.isin(list(ORIG) + ['BRCA_Normal'])].map({**ORIG, 'BRCA_Normal': 'Normal-like'})
print('original groups', g_orig.value_counts().to_dict(), 'n', len(g_orig))

def kw_fdr(samples, groups):
    M = X[samples].to_numpy(float); gs = np.asarray(groups); P = np.ones(len(M))
    idx = {g: np.where(gs == g)[0] for g in pd.unique(gs)}
    for i, row in enumerate(M):
        arrs = [row[k] for k in idx.values() if len(k) >= 3 and row[k].std(ddof=1) > 0]
        if len(arrs) >= 2: P[i] = stats.kruskal(*arrs).pvalue
    return pd.Series(multipletests(P, method='fdr_bh')[1], index=X.index)

q0 = kw_fdr(list(g_orig.index), list(g_orig.values))
diff = (q0 - arch.subtype_kruskal_FDR.reindex(q0.index)).abs()
rna_filter = (arch.mean_primary_tpm_log > 1) & (arch.detected_fraction_expr_gt1 > 0.20)
rna0 = (rna_filter & (q0 <= 0.05)).astype(float)
a = pd.read_csv(C/'A7_Basal_TNBC_H3K27ac_layer/A7_all_tf_scores_corrected_basal_H3K27ac.csv').set_index('TF')
agree = int((rna0.reindex(a.index).fillna(0) == a.RNA).sum())
print('max |FDR - archived|', diff.max(), '; RNA layer agreement', agree, 'of', len(a))
fig1e = pd.read_csv(W/'source_tables/figure1/figure1e_all_tf_subtype_association.csv').set_index('TF')
T1 = pd.DataFrame({'selection_FDR_4_original_groups': q0.reindex(F), 'RNA_layer': a.RNA.reindex(F),
                   'Fig1E_FDR_collapsed_contexts_plus_Normal_like': fig1e.subtype_kruskal_FDR.reindex(F)}).round(6)
T1.index.name = 'TF'; T1.to_csv(O/'A21_1_grouping_by_test.tsv', sep='\t'); print(T1.to_string())
pd.DataFrame([dict(test='RNA and activity layers (selection)', groups='LumA, LumB, HER2-enriched, Basal-like (Normal-like excluded)',
                   n_tumors=len(g_orig), fdr_family=f'all {len(X)} TFs in the TCGA TF expression matrix', script='figure1a.py, kruskal_by_tf'),
              dict(test='Fig. 1E ranking', groups='Luminal (A+B), HER2, Basal/TNBC, Normal-like', n_tumors=int(len(g_five)),
                   fdr_family='1,631 TFs tested in figure1e.py', script='figure1e.py')]).to_csv(O/'A21_1_test_definitions.tsv', sep='\t', index=False)

core, CTX = 'Core_discovery_score_v2_0_7', ['Luminal', 'HER2', 'Basal/TNBC']
SUP = ['scRNA_malignant', 'scATAC_cis_accessibility', 'DepMap_breast_dependency']
ST = a[['state_score_Luminal', 'state_score_HER2', 'state_score_Basal/TNBC']].to_numpy()
def call(st, c, t):
    s = dict(zip(CTX, st)); mx = max(st); sec = sorted(st)[-2]
    P = [k for k, v in s.items() if v >= mx - 0.20 and v >= 2.80]
    if not P and c >= 5.50 and mx >= 2.60: P = [max(s, key=s.get)]
    if not P and c >= 4.50 and mx >= 2.35: P = [max(s, key=s.get)]
    S = [k for k, v in s.items() if k not in P and v >= 2.35 and v >= 0.82*mx]
    if c >= 5.50 and mx >= 2.60: tier = 'High-confidence'
    elif (c >= 4.75 and mx >= 3.00 and t - c >= 0.50) or (c >= 4.75 and mx >= 3.50 and mx - sec >= 1.00) or (t >= 6.00 and mx >= 2.40): tier = 'Provisional'
    else: tier = 'Supporting'
    return tier, (tier != 'Supporting') and bool(P + S)
def formal(rna):
    c = (a[core] - a.RNA + rna.reindex(a.index).fillna(0)).to_numpy(); s = a[SUP].sum(1).to_numpy()
    return {t for i, t in enumerate(a.index) if call(list(ST[i]), c[i], c[i] + s[i])[1]}
assert formal(a.RNA) == set(F)

rng = np.random.default_rng(SEED)
counts = g_orig.value_counts(); prop = (counts/counts.sum()*312).round().astype(int)
designs = [('a balanced: 78 per original group', g_orig, {g: 78 for g in counts.index}),
           ('b proportional, same total (312)', g_orig, prop.to_dict()),
           ('c five classes incl. Normal-like, 36 each', g_five, {g: 36 for g in g_five.unique()})]
rows, freq = [], {}
for name, lab_, k in designs:
    kept, rates = [], []
    for _ in range(REPS):
        cols, gr = [], []
        for g, n in k.items():
            ids = lab_.index[lab_ == g]; s_ = rng.choice(ids, n, replace=False); cols += list(s_); gr += [g]*n
        q = kw_fdr(cols, gr); rna = (rna_filter & (q <= 0.05)).astype(float)
        kept.append(formal(rna)); rates.append(rna.reindex(F))
    R = pd.concat(rates, axis=1).mean(1); fq = pd.Series({t: np.mean([t in m for m in kept]) for t in a.index})
    new = fq[~fq.index.isin(F)].sort_values(ascending=False)
    freq[name] = fq.reindex(F); freq[name + ' | RNA layer rate'] = R
    rows.append(dict(design=name, group_sizes='; '.join(f'{g} {n}' for g, n in k.items()), draws=REPS,
                     all_13_retained=np.mean([set(F) <= m for m in kept]), median_formal=np.median([len(m) for m in kept]),
                     median_new=np.median([len(m - set(F)) for m in kept]), max_new_TF_rate=round(float(new.iloc[0]), 3) if len(new) else 0,
                     ELF3=fq['ELF3'], EHF=fq['EHF'], TFAP2C_RNA=R['TFAP2C']))
T3 = pd.DataFrame(rows); T3.to_csv(O/'A21_3_subsampling_summary.tsv', sep='\t', index=False); print(T3.to_string())
pd.DataFrame(freq).round(3).rename_axis('TF').to_csv(O/'A21_3_subsampling_TF_frequency.tsv', sep='\t'); print(pd.DataFrame(freq).round(3).to_string())

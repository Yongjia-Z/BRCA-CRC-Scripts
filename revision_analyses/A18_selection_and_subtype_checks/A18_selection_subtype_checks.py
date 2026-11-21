import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
W, C, A16, O = map(Path, sys.argv[1:5]); O.mkdir(parents=True, exist_ok=True)
pd.set_option('display.width', 230, 'display.max_columns', 40, 'display.max_rows', 120, 'display.max_colwidth', 90)
SEED = 20261101
a = pd.read_csv(C/'A7_Basal_TNBC_H3K27ac_layer/A7_all_tf_scores_corrected_basal_H3K27ac.csv').set_index('TF')
F = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']
core, tot = 'Core_discovery_score_v2_0_7', 'Audited_total_score_0_10'
CORE = ['RNA', 'TF_activity_regulon', 'ATAC_motif', 'ATAC_footprint', 'ATAC_peak_to_gene_accessibility', 'context_matched_H3K27ac_SE_MAX', 'Direct_binding_MAX']
SUP = ['scRNA_malignant', 'scATAC_cis_accessibility', 'DepMap_breast_dependency']
O5 = ['RNA', 'TF_activity_regulon', 'ATAC_motif', 'ATAC_peak_to_gene_accessibility', 'context_matched_H3K27ac_SE_MAX']
CTX = ['Luminal', 'HER2', 'Basal/TNBC']
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
def replay(c, s):
    res = [call(list(ST[i]), c[i], c[i] + s[i]) for i in range(len(a))]
    return pd.DataFrame(res, index=a.index, columns=['tier', 'formal'])
d0 = replay(a[core].to_numpy(), a[SUP].sum(1).to_numpy()); assert set(d0.index[d0.formal]) == set(F)
e = pd.read_csv(W/'source_tables/figure1/figure1e_all_tf_subtype_association.csv').sort_values('subtype_kruskal_FDR').reset_index(drop=True)
e['rank_fig1E'] = np.arange(1, len(e) + 1); e = e.set_index('TF')
print('TFs tested', len(e), 'FDR<=0.05', int((e.subtype_kruskal_FDR <= 0.05).sum()), 'RNA layer = 1', int((a.RNA == 1).sum()))
def why(t):
    if t not in a.index: return 'not in score table'
    c, s_, m = a.loc[t, core], a.loc[t, 'Supportive_score_0_3'], a.loc[t, 'max_state_score_0_5']
    if t in F: return f"{a.loc[t, 'confidence_tier']}: " + ('core >= 5.50, max context >= 2.60' if c >= 5.5 and m >= 2.6 else 'provisional rule met')
    miss = [k for k in ['ATAC_motif', 'ATAC_footprint', 'Direct_binding_MAX', 'TF_activity_regulon', 'ATAC_peak_to_gene_accessibility'] if a.loc[t, k] == 0]
    return f'core {c:.2f}, core + supportive {c + s_:.2f}, max context {m:.2f}: below thresholds; zero layers: ' + (', '.join(miss) or 'none')
top = list(e.index[:26]); rows = list(dict.fromkeys(top + F))
T1 = pd.DataFrame([dict(TF=t, rank_fig1E=int(e.loc[t, 'rank_fig1E']), subtype_FDR=e.loc[t, 'subtype_kruskal_FDR'], formal=t in F, **{k: round(a.loc[t, k], 3) if t in a.index else np.nan for k in CORE + [core, 'Supportive_score_0_3', 'max_state_score_0_5']}, explanation=why(t)) for t in rows])
T1.to_csv(O/'A18_1_fig1E_selection_reasons.tsv', sep='\t', index=False); print(T1[['TF', 'rank_fig1E', 'formal', core, 'Supportive_score_0_3', 'max_state_score_0_5', 'explanation']].to_string())
a['open5'] = a[O5].sum(1)
lines = {'MCF-7', 'T-47D', 'ZR-75-1', 'BT-474', 'SK-BR-3', 'HCC1954', 'MDA-MB-231', 'MDA-MB-468', 'BT-549', 'Hs 578T', 'HCC70', 'HCC38', 'SUM 159PT', 'HCC1806', 'HCC1395'}
n = {}
for b in open(A16/'raw/metadata/experimentList.tab', 'rb'):
    p = b.decode('utf-8', 'replace').split('\t')
    if len(p) < 9 or p[1] != 'hg38' or p[2] != 'TFs and others' or p[4] != 'Breast' or p[5] not in lines: continue
    n[p[3]] = n.get(p[3], 0) + 1
assert n
a['chip_profiles'] = pd.Series(n).reindex(a.index).fillna(0).astype(int)
a['chip_stratum'] = pd.cut(a.chip_profiles, [-1, 0, 4, 1e9], labels=['0', '1-4', '>=5'])
print(a.groupby('chip_stratum', observed=True).agg(n=('open5', 'size'), mean_open5=('open5', 'mean')), a.loc[F, 'chip_stratum'].value_counts())
rng = np.random.default_rng(SEED); obs = a.loc[F, 'open5'].mean(); N = 2000
un = np.array([a.open5.values[rng.choice(len(a), 13, replace=False)].mean() for _ in range(N)])
mt = []
for _ in range(N):
    v = []
    for st, k in a.loc[F, 'chip_stratum'].value_counts().items():
        if k: v += list(rng.choice(a.open5[a.chip_stratum == st].values, k, replace=False))
    mt.append(np.mean(v))
mt = np.array(mt)
T2 = pd.DataFrame([dict(reference='random 13-TF sets', n_sets=N, observed_mean_5layer=obs, median=np.median(un), p95=np.percentile(un, 95), max=un.max(), empirical_p=(1 + (un >= obs).sum())/(N + 1)),
                   dict(reference='13-TF sets matched on breast ChIP-Atlas profile strata (0, 1-4, >=5)', n_sets=N, observed_mean_5layer=obs, median=np.median(mt), p95=np.percentile(mt, 95), max=mt.max(), empirical_p=(1 + (mt >= obs).sum())/(N + 1))])
T2.to_csv(O/'A18_2_five_layer_reference.tsv', sep='\t', index=False); print(T2.round(4).to_string())
a[['open5', 'chip_profiles', 'chip_stratum']].to_csv(O/'A18_2_five_layer_scores.tsv', sep='\t')
print('TFs with >=1 breast ChIP profile', int((a.chip_profiles > 0).sum()), '; mean core by stratum', a.groupby('chip_stratum', observed=True)[core].mean().round(3).to_dict(), '; formal by stratum', a.assign(f=a.index.isin(F)).groupby('chip_stratum', observed=True).f.sum().to_dict())
print('spearman chip_profiles vs core', stats.spearmanr(a.chip_profiles, a[core]), 'vs open5', stats.spearmanr(a.chip_profiles, a.open5))
tpm = pd.read_csv(W/'datasets/TCGA/data/expression/TCGA-BRCA.star_tpm.tsv.gz', sep='\t')
tpm['g'] = tpm.Ensembl_ID.astype(str).str.replace(r'\..*$', '', regex=True)
ann = pd.read_csv(W/'datasets/GENCODE/results/gencode_v36_gene_annotation.tsv', sep='\t')[['ensembl_gene_id', 'gene_name']].drop_duplicates('ensembl_gene_id')
sc = [c for c in tpm.columns if c.startswith('TCGA')]
ex = tpm.merge(ann, left_on='g', right_on='ensembl_gene_id').groupby('gene_name')[sc].mean(); ex = ex.loc[[t for t in a.index if t in ex.index]]
prim = [c for c in sc if c[13:15] == '01']
cli = pd.read_csv(W/'datasets/TCGA/data/cbioportal_brca_patient_clinical_wide.tsv', sep='\t').set_index('patientId').SUBTYPE
lab = pd.Series({c: cli.get(c[:12]) for c in prim}).dropna()
g4 = {'BRCA_LumA': 'LumA', 'BRCA_LumB': 'LumB', 'BRCA_Her2': 'HER2', 'BRCA_Basal': 'Basal', 'BRCA_Normal': 'Normal-like'}
lab = lab.map(g4).dropna(); print(lab.value_counts().to_dict(), 'primary', len(prim))
pass0 = (ex[prim].mean(1) > 1) & ((ex[prim] > 1).mean(1) > 0.20)
X = ex.loc[pass0, lab.index]
def kw_fdr(cols, groups):
    M = X[cols].to_numpy(); R = stats.rankdata(M, axis=1); N_ = M.shape[1]; gs = np.array(groups); H = np.zeros(M.shape[0])
    for g in np.unique(gs):
        k = gs == g; H += R[:, k].sum(1)**2/k.sum()
    H = 12/(N_*(N_ + 1))*H - 3*(N_ + 1); p = stats.chi2.sf(H, len(np.unique(gs)) - 1)
    o = np.argsort(p); q = np.empty_like(p); q[o] = np.minimum.accumulate((p[o]*len(p)/np.arange(1, len(p) + 1))[::-1])[::-1]
    return pd.Series(q, index=X.index)
rows = []; freq = {}
rng = np.random.default_rng(SEED)
q_full = kw_fdr(list(lab.index), list(lab.values)); print('full cohort five-class FDR, formal TFs', q_full.reindex(F).round(4).to_dict(), 'layer agreement', int(((q_full <= 0.05).reindex(a.index).fillna(False) == (a.RNA == 1)).sum()), 'of', len(a))
for name, groups, k in [('five PAM50 classes, n = 36 each', ['LumA', 'LumB', 'HER2', 'Basal', 'Normal-like'], 36), ('four classes without Normal-like, n = 78 each', ['LumA', 'LumB', 'HER2', 'Basal'], 78)]:
    kept = []; rna = []
    for rep in range(200):
        cols, gr = [], []
        for g in groups:
            ids = lab.index[lab == g]; s_ = rng.choice(ids, min(k, len(ids)), replace=False); cols += list(s_); gr += [g]*len(s_)
        q = kw_fdr(cols, gr); new = (q <= 0.05).reindex(a.index, fill_value=False).astype(float)
        c = (a[core] - a.RNA + new).to_numpy(); d = replay(c, a[SUP].sum(1).to_numpy()); fm = set(d.index[d.formal])
        kept.append(fm); rna.append(new)
    R_ = pd.concat(rna, axis=1).mean(1); fq = pd.Series([sum(t in m for m in kept)/len(kept) for t in a.index], index=a.index)
    freq[name] = fq; freq[name + ' (RNA layer rate)'] = R_
    rows.append(dict(design=name, reps=200, median_formal_set=np.median([len(m) for m in kept]), all_13_retained_fraction=np.mean([set(F) <= m for m in kept]), median_13_retained=np.median([len(set(F) & m) for m in kept]), median_new_TFs=np.median([len(m - set(F)) for m in kept]), rna_layer_rate_13_min=R_[F].min(), tf_lowest_rna=R_[F].idxmin(), most_frequent_new=';'.join(f'{i}:{v:.2f}' for i, v in fq[~fq.index.isin(F)].sort_values(ascending=False).head(5).items())))
    print(name, 'selection', fq[F].round(2).to_dict(), 'RNA', R_[F].round(2).to_dict())
T3 = pd.DataFrame(rows); T3.to_csv(O/'A18_3_balanced_subtype_replay.tsv', sep='\t', index=False); print(T3.to_string())
pd.DataFrame(freq).loc[F].round(3).to_csv(O/'A18_3_balanced_subtype_TF_frequency.tsv', sep='\t')
ms = pd.read_csv(W/'datasets/METABRIC/data/metabric_sample_clinical_wide.tsv', sep='\t', low_memory=False); mp = pd.read_csv(W/'datasets/METABRIC/data/metabric_patient_clinical_wide.tsv', sep='\t', low_memory=False)
if 'patientId' not in ms.columns: ms['patientId'] = ms.sampleId
mm = ms.merge(mp, on='patientId')
mm = pd.DataFrame(dict(sample=mm.sampleId, pam50=mm.CLAUDIN_SUBTYPE, ER=mm.ER_STATUS.map({'Positive': 1, 'Negative': 0}), PR=mm.PR_STATUS.map({'Positive': 1, 'Negative': 0}), HER2=mm.HER2_STATUS.map({'Positive': 1, 'Negative': 0})))
gs = pd.read_csv(W/'datasets/GSE96058/results/GSE96058_CRC_module_scores.tsv', sep='\t')
ss = pd.DataFrame(dict(sample=gs.geoAcc, title=gs.title, repl=gs.isRepl, pam50=gs.pam50_subtype, ER=gs.er_Status, PR=gs.pgr_Status, HER2=gs.her2_Status)); ss = ss[ss.repl.astype(str) == 'False']
def rec(r):
    if pd.isna(r.ER) or pd.isna(r.HER2): return np.nan
    if r.HER2 == 1: return 'HER2+'
    if r.ER == 1 or r.PR == 1: return 'HR+/HER2-'
    if r.ER == 0 and r.HER2 == 0 and pd.notna(r.PR) and r.PR == 0: return 'TNBC'
    return np.nan
cr = []; store = {}
for co, x in [('METABRIC', mm), ('GSE96058', ss)]:
    x = x.copy(); x['pam'] = x.pam50.replace({'LumA': 'Luminal', 'LumB': 'Luminal', 'Her2': 'HER2-enriched', 'Basal': 'Basal-like', 'Normal': 'Normal-like'}); x['receptor'] = x.apply(rec, axis=1)
    store[co] = x; y = x.dropna(subset=['receptor']); y = y[y.pam.isin(['Luminal', 'HER2-enriched', 'Basal-like', 'Normal-like'])]
    tab = pd.crosstab(y.pam, y.receptor); print(co, len(y)); print(tab)
    for p_, r_ in [('Luminal', 'HR+/HER2-'), ('HER2-enriched', 'HER2+'), ('Basal-like', 'TNBC')]:
        cr.append(dict(cohort=co, n_tumors=len(y), pam50=p_, receptor=r_, n_pam50=int(tab.loc[p_].sum()), n_receptor=int(tab[r_].sum()), n_both=int(tab.loc[p_, r_]), pct_pam50_with_receptor=round(100*tab.loc[p_, r_]/tab.loc[p_].sum(), 1), pct_receptor_with_pam50=round(100*tab.loc[p_, r_]/tab[r_].sum(), 1)))
    tab.assign(cohort=co).to_csv(O/f'A18_4_crosstab_{co}.tsv', sep='\t')
T4 = pd.DataFrame(cr); T4.to_csv(O/'A18_4_pam50_receptor_overlap.tsv', sep='\t', index=False); print(T4.to_string())
ext = pd.read_csv(W/'revision_tables/R1_7_external_state_scores.tsv', sep='\t')
print(ext.head(3).to_string())
out = []
for co, x in store.items():
    e_ = ext[ext.cohort == co]; ids = set(e_.sample_id); key = 'sample'
    if not ids & set(x['sample']) and 'title' in x: key = 'title'
    print(co, 'id match', len(ids & set(x[key])))
    e2 = e_.merge(x[[key, 'receptor']], left_on='sample_id', right_on=key).dropna(subset=['receptor'])
    for stt, g in e2.groupby('state_id'):
        med = g.groupby('receptor').state_score.median(); H, p = stats.kruskal(*[v.state_score.values for _, v in g.groupby('receptor')])
        out.append(dict(cohort=co, state_id=stt, n=len(g), **{f'median_{k}': round(v, 3) for k, v in med.items()}, top_receptor_group=med.idxmax(), kw_p=p))
T5 = pd.DataFrame(out); T5.to_csv(O/'A18_4_states_by_receptor_group.tsv', sep='\t', index=False); print(T5.to_string())
sc_ = {'all layers (manuscript)': (a[core], a[SUP].sum(1)),
       'without DepMap': (a[core], a[SUP].sum(1) - a.DepMap_breast_dependency),
       'without GSE176078 malignant enrichment': (a[core], a[SUP].sum(1) - a.scRNA_malignant),
       'without scATAC': (a[core], a[SUP].sum(1) - a.scATAC_cis_accessibility),
       'without all three supportive layers': (a[core], a[SUP].sum(1)*0)}
rows = []
for k, (c, s_) in sc_.items():
    d = replay(c.to_numpy(), s_.to_numpy()); fm = set(d.index[d.formal]); hi = set(d.index[d.formal & (d.tier == 'High-confidence')])
    rows.append(dict(scenario=k, n_formal=len(fm), n_high=len(hi), lost=';'.join(sorted(set(F) - fm)), added=';'.join(sorted(fm - set(F)))))
T6 = pd.DataFrame(rows); T6.to_csv(O/'A18_5_reuse_replay.tsv', sep='\t', index=False); print(T6.to_string())

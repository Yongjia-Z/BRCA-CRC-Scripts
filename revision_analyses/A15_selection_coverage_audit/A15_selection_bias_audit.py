import sys, re
from pathlib import Path
import numpy as np, pandas as pd
W, C, O, EXP = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]); O.mkdir(parents=True, exist_ok=True)
pd.set_option('display.width', 200, 'display.max_columns', 30, 'display.max_rows', 100, 'display.max_colwidth', 120)
a = pd.read_csv(C/'A7_Basal_TNBC_H3K27ac_layer/A7_all_tf_scores_corrected_basal_H3K27ac.csv'); a = a.set_index(a.columns[0])
F = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']
C16 = ['ESR1','FOXA1','SPDEF','GATA3','PBX1','AR','CREB3L4','TFAP2C','EHF','MYBL2','NR2F6','STAT1','ELF3','RUNX1','TP63','FOXC1']
C18 = ['ESR1','FOXA1','GATA3','XBP1','TFAP2A','TFAP2C','TP63','FOXC1','SOX9','RUNX1','RUNX2','JUN','FOS','NFKB1','RELA','STAT3','TEAD1','MYC']
core, tot = 'Core_discovery_score_v2_0_7', 'Audited_total_score_0_10'
CORE = ['RNA', 'TF_activity_regulon', 'ATAC_motif', 'ATAC_footprint', 'ATAC_peak_to_gene_accessibility', 'context_matched_H3K27ac_SE_MAX', 'Direct_binding_MAX']
SUP = ['scRNA_malignant', 'scATAC_cis_accessibility', 'DepMap_breast_dependency']
assert (a[CORE].sum(1) - a[core]).abs().max() < 1e-9
cov = pd.DataFrame({'layer_role': ['core']*7 + ['supportive']*3, 'n_nonzero': [(a[c] != 0).sum() for c in CORE + SUP], 'pct_nonzero': [round((a[c] != 0).mean()*100, 1) for c in CORE + SUP],
                    'n_at_max_1': [(a[c] == 1).sum() for c in CORE + SUP], 'n_formal_13_nonzero': [(a.loc[F, c] != 0).sum() for c in CORE + SUP],
                    'nonzero_outside_formal_13': [';'.join(a.index[(a[c] != 0) & ~a.index.isin(F)][:12]) for c in CORE + SUP]}, index=CORE + SUP)
cov.to_csv(O/'A15_1_layer_coverage.tsv', sep='\t'); print(cov.to_string())
print('RNA=activity=motif=1:', int(((a.RNA == 1) & (a.TF_activity_regulon == 1) & (a.ATAC_motif == 1)).sum()), '; core >= 5.5:', a.index[a[core] >= 5.5].tolist(), '; core>=5.5 without footprint or binding:', int(((a[core] >= 5.5) & (a.ATAC_footprint == 0) & (a.Direct_binding_MAX == 0)).sum()))
print('formal in C16:', len(set(F) & set(C16)), sorted(set(F) - set(C16)), 'C16 not formal:', sorted(set(C16) - set(F)), 'formal in C18:', sorted(set(F) & set(C18)))
print('scATAC nonzero outside C16:', sorted(set(a.index[a.scATAC_cis_accessibility != 0]) - set(C16)), 'footprint nonzero outside C18:', sorted(set(a.index[a.ATAC_footprint != 0]) - set(C18)), 'binding nonzero outside C16:', sorted(set(a.index[a.Direct_binding_MAX != 0]) - set(C16)))
O5 = ['RNA', 'TF_activity_regulon', 'ATAC_motif', 'ATAC_peak_to_gene_accessibility', 'context_matched_H3K27ac_SE_MAX']
a['open5'] = a[O5].sum(1); a['rank_open5'] = a.open5.rank(ascending=False, method='min').astype(int)
r = a.sort_values(['open5', core], ascending=False)[O5 + ['open5', 'rank_open5', core, 'max_state_score_0_5']].round(3); r['formal_13'] = r.index.isin(F); r['preselected_16'] = r.index.isin(C16)
r.head(120).to_csv(O/'A15_2_rank_open5_top120.tsv', sep='\t'); print(a.loc[F, ['open5', 'rank_open5', core]].round(3).sort_values('rank_open5').to_string())
print('formal in top 13/20/50/100:', [int(r.index[:k].isin(F).sum()) for k in (13, 20, 50, 100)], 'ties at open5>=4.75:', int((a.open5 >= 4.75).sum()))
fp = pd.read_csv(W/'datasets/TCGA/results/ATAC_footprint_expression_correlation.tsv', sep='\t')
for c in ['depth_FDR', 'flank_FDR', 'flank_correlation']: fp[c] = pd.to_numeric(fp[c], errors='coerce')
key = {t.replace('-', '').upper(): t for t in a.index}
fp['TF'] = fp.TF_Name.astype(str).str.upper().map(lambda x: key.get(x.replace('-', '')))
best = fp.dropna(subset=['TF']).sort_values('flank_correlation', ascending=False).drop_duplicates('TF').set_index('TF')
uf = best.flank_FDR.lt(0.05).astype(float).reindex(a.index).fillna(0.0); assert (uf[a.ATAC_footprint == 1] == 1).all()
print('motifs in S6 correlation sheet', len(fp), 'mapped to TF universe', best.shape[0], 'uniform footprint TFs', int(uf.sum()), 'previously', int((a.ATAC_footprint == 1).sum()))
pd.DataFrame({'flank_correlation': best.flank_correlation, 'flank_FDR': best.flank_FDR, 'footprint_in_manuscript': a.ATAC_footprint.reindex(best.index), 'footprint_uniform': uf.reindex(best.index)}).query('footprint_uniform == 1').sort_values('flank_FDR').round(4).to_csv(O/'A15_3_uniform_footprint.tsv', sep='\t')
CTX = ['Luminal', 'HER2', 'Basal/TNBC']
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
    st = a[['state_score_Luminal', 'state_score_HER2', 'state_score_Basal/TNBC']].to_numpy()
    res = [call(list(st[i]), c.iloc[i], c.iloc[i] + s.iloc[i]) for i in range(len(a))]
    d = pd.DataFrame(res, index=a.index, columns=['tier', 'formal']); return d
s0 = a[SUP].sum(1); s_no = a.scRNA_malignant + a.DepMap_breast_dependency
LINES = {'MCF-7', 'T-47D', 'ZR-75-1', 'BT-474', 'SK-BR-3', 'HCC1954', 'MDA-MB-231', 'MDA-MB-468', 'BT-549', 'Hs 578T', 'HCC70', 'HCC38', 'SUM 159PT', 'HCC1806', 'HCC1395'}
n, n3 = {}, {}
for b in open(EXP, 'rb'):
    p = b.decode('utf-8', 'replace').split('\t')
    if len(p) < 9 or p[1] != 'hg38' or p[2] != 'TFs and others' or p[4] != 'Breast' or p[5] not in LINES: continue
    m = re.search(r',(\d+)\s*$', p[7]); k = int(m.group(1)) if m else 0; n[p[3]] = n.get(p[3], 0) + 1
    if k >= 3000: n3[p[3]] = n3.get(p[3], 0) + 1
av = pd.DataFrame({'chipatlas_breast_profiles': pd.Series(n), 'profiles_ge3000_peaks': pd.Series(n3)}).fillna(0).astype(int)
a = a.join(av).fillna({'chipatlas_breast_profiles': 0, 'profiles_ge3000_peaks': 0}); assert n
z = a[(a.Direct_binding_MAX == 0) & (a.profiles_ge3000_peaks > 0)]
z.sort_values(core, ascending=False)[[core, 'max_state_score_0_5', 'chipatlas_breast_profiles', 'profiles_ge3000_peaks']].round(3).to_csv(O/'A15_4_chip_available_binding_zero.tsv', sep='\t')
print('TFs with any breast-line profile', int((a.chipatlas_breast_profiles > 0).sum()), '>=3000 peaks', int((a.profiles_ge3000_peaks > 0).sum()), 'binding 0 with >=3000-peak profile', len(z))
bind = lambda v: a.Direct_binding_MAX.where(a.Direct_binding_MAX > 0, (a.profiles_ge3000_peaks > 0).astype(float)*v)
SC = {'S0 manuscript': (a[core], s0),
      'S1 uniform footprint': (a[core] - a.ATAC_footprint + uf, s0),
      'S2 no scATAC in supportive': (a[core], s_no),
      'S3 uniform footprint, no scATAC': (a[core] - a.ATAC_footprint + uf, s_no),
      'S4 no footprint, no binding, no scATAC; core rescaled 7/5': (a.open5*7/5, s_no),
      'S5 S3 + binding 0.50 where >=3000-peak breast ChIP exists': (a[core] - a.ATAC_footprint + uf - a.Direct_binding_MAX + bind(0.50), s_no),
      'S6 S3 + binding 0.75 where >=3000-peak breast ChIP exists': (a[core] - a.ATAC_footprint + uf - a.Direct_binding_MAX + bind(0.75), s_no)}
rows = []; member = pd.DataFrame(index=a.index)
for k, (c, s) in SC.items():
    d = replay(c, s); fm = d.index[d.formal].tolist(); hi = d.index[d.formal & d.tier.eq('High-confidence')].tolist(); member[k] = d.tier.where(d.formal, '')
    rows.append(dict(scenario=k, n_formal=len(fm), n_high_confidence=len(hi), formal13_kept=len(set(F) & set(fm)), formal13_lost=';'.join(sorted(set(F) - set(fm))), n_added=len(set(fm) - set(F)), high_confidence=';'.join(sorted(hi)), added=';'.join(sorted(set(fm) - set(F)))))
R = pd.DataFrame(rows); R.to_csv(O/'A15_5_tier_replay.tsv', sep='\t', index=False); print(R.drop(columns=['added']).to_string())
assert set(member.index[member['S0 manuscript'] != '']) == set(F), 'tier rule does not reproduce the 13 TFs'
member = member[(member != '').any(axis=1)]; member['chipatlas_breast_profiles'] = a.chipatlas_breast_profiles; member.to_csv(O/'A15_5_tier_replay_membership.tsv', sep='\t')
print(R[['scenario', 'added']].to_string())
print('near thresholds:', a.loc[F, [core, 'max_state_score_0_5', tot]].assign(support=lambda d: d[tot] - d[core]).round(3).to_string())

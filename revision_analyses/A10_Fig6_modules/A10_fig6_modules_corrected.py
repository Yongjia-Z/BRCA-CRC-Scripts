import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT, A7, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
F6 = ROOT/'core_data/remap_reanalysis/figure5-8/figure6'
TC = ROOT/'datasets/TCGA/results'
CRC = ["ESR1","GATA3","PGR","PBX1","FOXA1","SPDEF","ELF3","EHF","AR","TFAP2C","RUNX1","TP63","FOXC1"]
SC = ['pan_CRC13_score','luminal_CRC13_score','her2_CRC13_score','basal_TNBC_CRC13_score']
LAB = {'pan_CRC13_score':'Pan CRC TFs','luminal_CRC13_score':'Luminal CRC TFs','her2_CRC13_score':'HER2 CRC TFs','basal_TNBC_CRC13_score':'Basal/TNBC CRC TFs'}
def zrows(df):
    v = df.astype(float); sd = v.std(axis=1).replace(0, np.nan)
    return v.sub(v.mean(axis=1), axis=0).div(sd, axis=0).fillna(0)
expr = pd.read_csv(TC/'TCGA_BRCA_TF_expression_matrix.tsv.gz', sep='\t', index_col=0)
act = pd.read_csv(TC/'TCGA_BRCA_DoRothEA_like_TF_activity_matrix.tsv.gz', sep='\t', index_col=0)
common = [s for s in expr.columns if s in act.columns]
ez = zrows(expr.loc[[g for g in CRC if g in expr.index], common]); az = zrows(act.loc[[g for g in CRC if g in act.index], common])
def modscore(genes):
    e = ez.loc[[g for g in genes if g in ez.index]].mean(0); a = az.loc[[g for g in genes if g in az.index]].mean(0)
    return (e + a) / 2
old_mem = json.loads((F6/'figure6_crc13_recompute_summary.json').read_text(encoding='utf8'))['module_membership']
t = pd.read_csv(A7); t = t[t.TF.isin(CRC)].set_index('TF').loc[CRC]
new_mem = {'luminal': [g for g in CRC if int(t.loc[g, 'Luminal_member']) == 1],
           'her2': [g for g in CRC if int(t.loc[g, 'HER2_member']) == 1],
           'basal_TNBC': [g for g in CRC if int(t.loc[g, 'Basal/TNBC_member']) == 1], 'pan': CRC}
key = {'pan':'pan_CRC13_score','luminal':'luminal_CRC13_score','her2':'her2_CRC13_score','basal_TNBC':'basal_TNBC_CRC13_score'}
stored = pd.read_csv(F6/'figure6_crc13_scores_with_tumor_state_mcp.tsv', sep='\t').set_index('sample_id')
chk = {}
for m, g in old_mem.items():
    s = modscore(g); chk[key[m]] = float((s.reindex(stored.index) - stored[key[m]]).abs().max())
print('reproduction of stored scores with old membership, max |diff|:', chk)
assert max(chk.values()) < 1e-9
d = stored.copy()
for m, g in new_mem.items(): d[key[m]] = modscore(g).reindex(d.index)
d = d.reset_index()
changed = {key[m]: (old_mem[m], new_mem[m]) for m in new_mem if old_mem[m] != new_mem[m]}
print('changed modules:', changed)
d['context'] = d.SUBTYPE.map({'BRCA_LumA':'Luminal','BRCA_LumB':'Luminal','BRCA_Her2':'HER2','BRCA_Basal':'Basal/TNBC','BRCA_Normal':'Normal-like'})
d.to_csv(OUT/'A10_fig6_scores_corrected_membership.tsv', sep='\t', index=False)
r_old_new = stats.spearmanr(stored['luminal_CRC13_score'].reindex(d.sample_id).to_numpy(), d['luminal_CRC13_score'])[0]
ORD = ['Luminal','HER2','Basal/TNBC','Normal-like']
rows = []
for s in SC:
    sub = d.dropna(subset=['context', s]); groups = [sub.loc[sub.context.eq(c), s] for c in ORD]
    H, p = stats.kruskal(*groups)
    for c, gv in zip(ORD, groups):
        rows.append({'module': LAB[s], 'context': c, 'n': len(gv), 'median': gv.median(), 'q1': gv.quantile(.25), 'q3': gv.quantile(.75), 'kruskal_H': H, 'kruskal_p': p})
med = pd.DataFrame(rows); med.to_csv(OUT/'A10_fig6a_medians_corrected.tsv', sep='\t', index=False)
so = stored.reset_index(); so['context'] = so.SUBTYPE.map({'BRCA_LumA':'Luminal','BRCA_LumB':'Luminal','BRCA_Her2':'HER2','BRCA_Basal':'Basal/TNBC','BRCA_Normal':'Normal-like'})
mo = so.groupby('context')['luminal_CRC13_score'].median()
SIG = ["ER signaling","HER2 amplicon","EMT","Cell cycle","MYC targets","DNA repair","Hypoxia","Inflammatory response","T cell signature","B cell signature","Macrophage signature","Stromal signature","Immune checkpoint"]
MCP = ['T cells','CD8 T cells','Cytotoxic lymphocytes','NK cells','B lineage','Monocytic lineage','Myeloid dendritic cells','Neutrophils','Endothelial cells','Fibroblasts']
def sp(df, ys, ylab):
    rr = []
    for x in SC:
        for y in ys:
            sub = df[[x, y]].dropna(); r, p = stats.spearmanr(sub[x], sub[y]); rr.append({'CRC_score': x, ylab: y, 'spearman_rho': r, 'pvalue': p, 'n': len(sub)})
    o = pd.DataFrame(rr); o['fdr'] = multipletests(o.pvalue, method='fdr_bh')[1]; return o
b = sp(d, SIG, 'signature'); b.to_csv(OUT/'A10_fig6b_signature_correlation_corrected.tsv', sep='\t', index=False)
c = sp(d, MCP, 'MCPcounter_population'); c.to_csv(OUT/'A10_fig6c_mcp_correlation_corrected.tsv', sep='\t', index=False)
def resid(df, cols, grp):
    X = np.column_stack([np.ones(len(df)), pd.get_dummies(grp.fillna('missing'), dtype=float).to_numpy()])
    o = pd.DataFrame(index=df.index)
    for col in cols:
        y = df[col].to_numpy(float); b_ = np.linalg.lstsq(X, y, rcond=None)[0]; o[col] = y - X @ b_
    return o
final6d = pd.read_csv(ROOT/'source_tables/figure6/figure6d.csv')
conv = {}
for name, grp in [('PAM50', so.SUBTYPE), ('collapsed_context', so.context)]:
    rs = sp(resid(so, SC + MCP, grp), MCP, 'MCPcounter_population')
    m = rs.merge(final6d, on=['CRC_score', 'MCPcounter_population'], suffixes=('', '_final'))
    conv[name] = float((m.spearman_rho - m.spearman_rho_final).abs().max())
best = min(conv, key=conv.get); print('6D convention reproduction max|diff|:', conv, '->', best)
grp = d.SUBTYPE if best == 'PAM50' else d.context
dd = sp(resid(d, SC + MCP, grp), MCP, 'MCPcounter_population'); dd.to_csv(OUT/'A10_fig6d_adjusted_mcp_correlation_corrected.tsv', sep='\t', index=False)
def oldnew(new, old_path, ylab):
    old = pd.read_csv(old_path, sep='\t') if str(old_path).endswith('.tsv') else old_path
    m = new[new.CRC_score.eq('luminal_CRC13_score')].merge(old[old.CRC_score.eq('luminal_CRC13_score')][['CRC_score', ylab, 'spearman_rho']], on=['CRC_score', ylab], suffixes=('_corrected', '_stored'))
    m['delta'] = m.spearman_rho_corrected - m.spearman_rho_stored; return m
ob = oldnew(b, sp(so, SIG, 'signature'), 'signature')
oc_old = sp(so, MCP, 'MCPcounter_population')
storedc = pd.read_csv(F6/'figure6c_crc13_mcp_counter_correlation.tsv', sep='\t').merge(oc_old, on=['CRC_score','MCPcounter_population'])
conv['6C_stored_reproduction'] = float((storedc.spearman_rho_x - storedc.spearman_rho_y).abs().max())
oc = oldnew(c, oc_old, 'MCPcounter_population')
od_old = sp(resid(so, SC + MCP, so.SUBTYPE if best == 'PAM50' else so.context), MCP, 'MCPcounter_population')
od = oldnew(dd, od_old, 'MCPcounter_population')
cmp = pd.concat([ob.assign(panel='6B').rename(columns={'signature': 'feature'}), oc.assign(panel='6C').rename(columns={'MCPcounter_population': 'feature'}), od.assign(panel='6D').rename(columns={'MCPcounter_population': 'feature'})])
cmp.to_csv(OUT/'A10_luminal_module_stored_vs_corrected.tsv', sep='\t', index=False)
summ = {'old_membership': old_mem, 'corrected_membership': new_mem, 'reproduction_check_max_abs_diff': chk, 'fig6d_adjustment_convention': best, 'fig6d_reproduction_max_abs_diff': conv,
        'spearman_old_vs_new_luminal_score': r_old_new, 'luminal_median_old': mo.to_dict(), 'max_abs_delta_rho': {p: float(g.delta.abs().max()) for p, g in cmp.groupby('panel')}}
(OUT/'A10_summary.json').write_text(json.dumps(summ, indent=2, ensure_ascii=False), encoding='utf8')
print(json.dumps(summ, indent=1, ensure_ascii=False)); print(med.round(3).to_string(index=False)); print(cmp.round(3).to_string(index=False))

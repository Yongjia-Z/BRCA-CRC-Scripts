import sys, gzip
from pathlib import Path
import numpy as np, pandas as pd
W, A16, O = map(Path, sys.argv[1:4])
T = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']; F = 100_000
d = pd.read_csv(A16/'results/A16_occupancy_by_dataset.tsv', sep='\t', encoding='utf-8-sig').fillna('')
d['use_strict'] = d.use_strict.astype(str).str.lower().eq('true'); d = d[d.use_strict & d.TF.isin(T)].copy()
print(d.groupby(['TF', 'source']).size().to_string())
g = pd.read_csv(W/'datasets/GENCODE/results/gencode_v36_gene_annotation.tsv', sep='\t'); print(g.columns.tolist()[:12])
tc, schr, sst, sen, sstr = 'gene_type', 'chromosome', 'start', 'end', 'strand'
g = g[g[tc].eq('protein_coding') & g[schr].astype(str).str.match(r'^chr([0-9]+|X|Y)$')].drop_duplicates('gene_name').copy()
g['tss'] = np.where(g[sstr].eq('+'), g[sst], g[sen]); print('genes', len(g))
W_ = {c: x.sort_values('tss') for c, x in g.groupby(schr)}
rows = []
for r in d.itertuples():
    p = W/r.path if not Path(r.path).is_absolute() else Path(r.path)
    if not p.exists(): p = Path(str(r.path).replace('\\', '/')); p = p if p.is_absolute() else W/p
    if not p.exists(): rows.append(dict(TF=r.TF, dataset_id=r.dataset_id, source=r.source, series=r.series, cell_line=r.cell_line, n_peaks=0, status='file missing')); continue
    op = gzip.open if str(p).endswith('.gz') else open; starts = {}
    with op(p, 'rt', errors='replace') as f:
        for l in f:
            if l.startswith(('#', 'track', 'browser')): continue
            a = l.split('\t')
            if len(a) < 3: continue
            try: starts.setdefault(a[0], []).append((int(a[1]) + int(a[2]))//2)
            except ValueError: continue
    n = sum(len(v) for v in starts.values())
    if n == 0: rows.append(dict(TF=r.TF, dataset_id=r.dataset_id, source=r.source, series=r.series, cell_line=r.cell_line, n_peaks=0, status='no peaks')); continue
    cnt = []
    for c, x in W_.items():
        m = np.sort(np.array(starts.get(c, []), dtype=np.int64)); t = x.tss.to_numpy()
        cnt.append(pd.Series(np.searchsorted(m, t + F, 'right') - np.searchsorted(m, t - F, 'left'), index=x.gene_name.values))
    cnt = pd.concat(cnt); self_n = int(cnt.get(r.TF, 0)); pct = float((cnt < self_n).mean() + 0.5*(cnt == self_n).mean())
    rows.append(dict(TF=r.TF, dataset_id=r.dataset_id, source=r.source, series=r.series, cell_line=r.cell_line, n_peaks=n, self_peaks_100kb=self_n, self_promoter_peaks=r.self_peaks, median_gene_window=float(cnt.median()), p95_gene_window=float(cnt.quantile(.95)), percentile=round(pct, 4), enriched=pct >= 0.95, status='ok'))
D = pd.DataFrame(rows); D.to_csv(O/'A19_self_locus_100kb_by_dataset.tsv', sep='\t', index=False); print(D.status.value_counts().to_dict())
ok = D[D.status == 'ok']
S = ok.groupby('TF').agg(profiles=('dataset_id', 'size'), with_peak_100kb=('self_peaks_100kb', lambda x: int((x > 0).sum())), enriched_top5pct=('enriched', 'sum'), median_percentile=('percentile', 'median'), with_promoter_peak=('self_promoter_peaks', lambda x: int((pd.to_numeric(x, errors='coerce') > 0).sum())), max_peaks=('n_peaks', 'max')).reindex(T); print(S.to_string())
S.to_csv(O/'A19_self_locus_100kb_summary.tsv', sep='\t')

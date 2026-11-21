import sys
from pathlib import Path
import pandas as pd
W, O = Path(sys.argv[1]), Path(sys.argv[2])
T = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']
g = pd.read_csv(W/'datasets/GENCODE/results/gencode_v36_gene_annotation.tsv', sep='\t').set_index('gene_name').loc[T]
win = {n: (r.chromosome, int(r.start if r.strand == '+' else r.end)-2000, int(r.start if r.strand == '+' else r.end)+2000) for n, r in g.iterrows()}
rows = []
for f in sorted((W/'datasets/CRC13_Breast_ChIPAtlas_expansion/data').glob('*.bed')):
    srx, tf, line, ctx = f.name.split('.')[:4]
    c = dict.fromkeys(T, 0); tot = 0
    for l in open(f):
        if l.startswith(('track', '#', 'browser')): continue
        p = l.split('\t', 3); tot += 1
        for n, (ch, a, b) in win.items():
            if p[0] == ch and int(p[1]) < b and int(p[2]) > a: c[n] += 1
    rows.append(dict(srx=srx, regulator_TF=tf, cell_line=line, context=ctx, n_peaks=tot, self_peaks=c[tf], **c))
d = pd.DataFrame(rows); d.to_csv(O/'A14_chipatlas_promoter_occupancy.tsv', sep='\t', index=False)
pd.set_option('display.width', 230, 'display.max_columns', 30)
print(d.drop(columns=['context']).to_string())

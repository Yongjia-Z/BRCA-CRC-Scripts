from __future__ import annotations
import csv, gzip, hashlib, json, os, time
from datetime import datetime, timezone
from pathlib import Path
A16 = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('CRC_PROJECT_ROOT', Path(__file__).resolve().parents[3]))
PRIOR = A16.parent
REMAP = ROOT / 'core_data/remap_reanalysis'
A14 = PRIOR / 'A14_occupancy_and_network_checks'
TFS = ['ESR1','GATA3','PGR','PBX1','FOXA1','SPDEF','ELF3','EHF','AR','TFAP2C','RUNX1','TP63','FOXC1']
C16 = ['ESR1','FOXA1','SPDEF','GATA3','PBX1','AR','CREB3L4','TFAP2C','EHF','MYBL2','NR2F6','STAT1','ELF3','RUNX1','TP63','FOXC1']
LINES = ['MCF-7','T-47D','ZR-75-1','BT-474','SK-BR-3','HCC1954','MDA-MB-231','MDA-MB-468','BT-549','Hs578T','HCC70','HCC38','SUM159PT','HCC1806','HCC1395']
CONTEXT = {'MCF-7':'Luminal','T-47D':'Luminal','ZR-75-1':'Luminal','BT-474':'Luminal/HER2','SK-BR-3':'HER2','HCC1954':'HER2', **{x:'Basal/TNBC' for x in LINES[6:]}}
def stamp(): return datetime.now(timezone.utc).isoformat()
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(2**20),b''): h.update(b)
    return h.hexdigest()
def save_json(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def save_tsv(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);rows=list(rows)
    fields=fields or (list(rows[0]) if rows else [])
    with open(path,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def open_text(path):
    return gzip.open(path,'rt',encoding='utf-8-sig',errors='replace') if str(path).endswith('.gz') else open(path,encoding='utf-8-sig',errors='replace')
def promoter_table():
    import pandas as pd
    g=pd.read_csv(ROOT/'datasets/GENCODE/results/gencode_v36_gene_annotation.tsv',sep='\t').set_index('gene_name');out=[]
    for tf in dict.fromkeys(TFS+C16):
        r=g.loc[tf];t=int(r.start if r.strand=='+' else r.end)
        out.append({'TF':tf,'chrom':r.chromosome,'tss_legacy':t,'start':t-2000,'end':t+2000,'strand':r.strand,'coordinate_note':'legacy gene-level TSS +/-2000, as A14; GTF 1-based TSS to BED0 minus-one sensitivity audited separately'})
    return out
def peak_hits(path,promoters=None):
    p=promoters or promoter_table();byc={}
    for r in p:byc.setdefault(r['chrom'],[]).append(r)
    hits=dict.fromkeys([r['TF'] for r in p],0);n=0;bad=0
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#','track','browser')):continue
            x=line.split()
            try:ch=x[0];a=int(x[1]);b=int(x[2]);assert 0<=a<b
            except (ValueError,IndexError,AssertionError):bad+=1;continue
            n+=1
            for r in byc.get(ch,[]):
                if a<r['end'] and b>r['start']:hits[r['TF']]+=1
    return n,bad,hits

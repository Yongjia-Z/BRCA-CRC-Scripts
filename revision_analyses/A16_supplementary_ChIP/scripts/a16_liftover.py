from a16_common import *
from bisect import bisect_right
import numpy as np

class ChainMap:
    def __init__(self,path):
        self.by_chrom={}
        with gzip.open(path,'rt') as f:
            for line in f:
                if line.startswith('chain '):
                    z=line.split();_,score,tn,ts,tstrand,ta,tb,qn,qs,qstrand,qa,qb,cid=z
                    assert tstrand=='+'
                    rec={'id':cid,'score':int(score),'src_start':int(ta),'src_end':int(tb),'qname':qn,'qsize':int(qs),'strand':qstrand,'blocks':[]}
                    self.by_chrom.setdefault(tn,[]).append(rec);t=int(ta);q=int(qa)
                elif line.strip():
                    z=list(map(int,line.split()));n=z[0]
                    rec['blocks'].append((t,t+n,q,q+n));t+=n;q+=n
                    if len(z)==3:t+=z[1];q+=z[2]
        for chrom,recs in self.by_chrom.items():
            for r in recs:r['starts']=[b[0] for b in r['blocks']]
    def interval(self,ch,a,b,min_match=.95):
        candidates=[]
        for r in self.by_chrom.get(ch,[]):
            if r['src_start']>=b or r['src_end']<=a:continue
            j=max(0,bisect_right(r['starts'],a)-1);covered=0;dest=[]
            for s,e,q,qe in r['blocks'][j:]:
                if s>=b:break
                lo=max(a,s);hi=min(b,e)
                if lo>=hi:continue
                covered+=hi-lo;x=q+lo-s;y=q+hi-s
                if r['strand']=='-':x,y=r['qsize']-y,r['qsize']-x
                dest.extend((x,y))
            if covered/(b-a)>=min_match:
                candidates.append({'chrom':r['qname'],'start':min(dest),'end':max(dest),'strand':r['strand'],'match_fraction':covered/(b-a),'chain_id':r['id'],'chain_score':r['score']})
        return candidates

def convert(name,build):
    p=A16/'raw/peaks'/name;chain=A16/'raw/metadata'/f'{build}ToHg38.over.chain.gz';mapper=ChainMap(chain)
    from pyliftover import LiftOver
    independent=LiftOver(str(chain))
    audit=[];passed=[];endpoint_mismatch=0
    with open_text(p) as f:
        for i,line in enumerate(f,1):
            if not line.strip() or line.startswith(('track','#','browser')):continue
            x=line.split()
            try:a=int(x[1]);b=int(x[2]);assert b>a>=0
            except (ValueError,IndexError,AssertionError):
                audit.append({'source_line':i,'status':'invalid_interval','source':line.strip()[:100]});continue
            cand=mapper.interval(x[0],a,b)
            r={'source_line':i,'source_chrom':x[0],'source_start':a,'source_end':b,'n_qualifying_chains':len(cand)}
            if len(cand)!=1:
                r['status']='unmapped_or_insufficient_match' if not cand else 'multiple_mapping';audit.append(r);continue
            z=cand[0];r.update(z);r['status']='mapped_unique';r['length_ratio']=(z['end']-z['start'])/(b-a)
            ends=[]
            for pos in [a,b-1]:
                matches=[v for v in independent.convert_coordinate(x[0],pos) or [] if v[0]==z['chrom'] and v[2]==z['strand'] and v[3]==z['chain_score']]
                ends.append(matches)
            if all(len(v)==1 for v in ends):
                expect=(min(ends[0][0][1],ends[1][0][1]),max(ends[0][0][1],ends[1][0][1])+1)
                if expect!=(z['start'],z['end']):endpoint_mismatch+=1
            audit.append(r);passed.append((z['chrom'],z['start'],z['end'],f'peak_{i}'))
    dest=A16/'derived/hg38'/name.replace(f'_{build}.bed.gz','_hg38.bed')
    with open(dest,'w',encoding='utf-8') as f:
        for row in sorted(set(passed)):f.write('\t'.join(map(str,row))+'\n')
    fields=sorted(set().union(*(r.keys() for r in audit)))
    audit_path=A16/'results'/f'A16_liftover_{name.removesuffix(".bed.gz")}.tsv'
    save_tsv(audit_path,audit,fields)
    assert endpoint_mismatch==0
    return {'source_path':str(p.relative_to(ROOT)),'source_sha256':sha(p),'chain_path':str(chain.relative_to(ROOT)),'chain_sha256':sha(chain),'dest_path':str(dest.relative_to(ROOT)),'dest_sha256':sha(dest),'audit_path':str(audit_path.relative_to(ROOT)),'input_intervals':len(audit),'mapped_unique':len(passed),'output_unique_intervals':len(set(passed)),'multiple_mapping':sum(r['status']=='multiple_mapping' for r in audit),'unmapped':sum(r['status']=='unmapped_or_insufficient_match' for r in audit),'independent_endpoint_mismatches':endpoint_mismatch,'min_match':.95,'source_build':build,'target_build':'hg38'}

if __name__=='__main__':
    rows=[convert('PBX1_GSM692743_hg18.bed.gz','hg18'),convert('RUNX1_GSE75070_hg19.bed.gz','hg19')]
    save_json(A16/'results/A16_liftover_verification.json',rows);print(json.dumps(rows,ensure_ascii=False),flush=True)

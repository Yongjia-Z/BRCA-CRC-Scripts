from a16_common import *
import pandas as pd, re, struct

if __name__=='__main__':
    p=A16/'raw/metadata/hg38.chrom.sizes'
    if not p.exists():raise RuntimeError('Cannot verify hg38 chromosome bounds without chromosome lengths')
    sizes={x[0]:int(x[1]) for s in p.read_text().splitlines() if len(x:=s.split())==2}
    assert sizes['chr1']==248956422
    d=pd.read_csv(A16/'results/A16_occupancy_by_dataset.tsv',sep='\t').fillna('')
    quality=[];cache={}
    for r in d.itertuples(index=False):
        f=ROOT/r.path
        if r.path not in cache:
            seen=set();counts={'n_peaks':0,'malformed':0,'unknown_chromosome':0,'beyond_chromosome_end':0,'exact_duplicate_intervals':0}
            unknown={}
            with open_text(f) as h:
                for s in h:
                    if not s.strip() or s.startswith(('#','track','browser')):continue
                    x=s.split()
                    try:
                        ch=x[0];start=int(x[1]);end=int(x[2]);assert 0<=start<end
                    except (ValueError,IndexError,AssertionError):counts['malformed']+=1;continue
                    counts['n_peaks']+=1
                    if ch not in sizes:counts['unknown_chromosome']+=1;unknown[ch]=unknown.get(ch,0)+1
                    elif end>sizes[ch]:counts['beyond_chromosome_end']+=1
                    key=(ch,start,end)
                    if key in seen:counts['exact_duplicate_intervals']+=1
                    seen.add(key)
            cache[r.path]={**counts,'unknown_chromosome_names':json.dumps(unknown),'sha256_verified':sha(f)==r.sha256}
        quality.append({'TF':r.TF,'dataset_id':r.dataset_id,'source':r.source,'path':r.path,**cache[r.path]})
    save_tsv(A16/'results/A16_interval_quality.tsv',quality)

    group=[]
    for tf,sub in d.groupby('TF'):
        for r in sub.itertuples(index=False):
            own=set(re.findall(r'GSE\d+|ENCSR[A-Z0-9]+',str(r.series)))
            peers=[]
            for s in sub.itertuples(index=False):
                if s.source==r.source:continue
                gs=set(re.findall(r'GSE\d+|ENCSR[A-Z0-9]+',str(s.series)))
                if own & gs:peers.append(f'{s.source}:{s.dataset_id}')
            group.append({'TF':tf,'source':r.source,'dataset_id':r.dataset_id,'GSM':r.gsm,'SRX':r.srx,
                'series':r.series,'same_TF_series_other_processing_rows':';'.join(peers),
                'relationship':'same TF and GEO series; review sample/condition before asserting identical assay' if peers else 'no matched study identifier in other processing source',
                'identity_caution':'pooled author tables and per-replicate database tables are not independent replication'})
    save_tsv(A16/'results/A16_processing_correspondence.tsv',group)

    co=pd.read_csv(A16/'results/A16_coordinate_1bp_sensitivity.tsv',sep='\t').fillna('')
    details=[]
    for r in co[co.n_changed_promoters_1bp_conversion.gt(0)].itertuples(index=False):
        old=json.loads(r.legacy_counts);new=json.loads(r.BED0_counts)
        details.append({'dataset_id':r.dataset_id,'TF':r.TF,'legacy_counts':r.legacy_counts,'BED0_counts':r.BED0_counts,
            'binary_promoter_presence_changed':any(bool(old[k])!=bool(new[k]) for k in old),
            'delta_old16_dataset_promoter_count':sum(int(new[k]>0)-int(old[k]>0) for k in old if k in C16),
            'note':'one-base coordinate correction is audited separately; historical windows are preserved for the main replay'})
    save_tsv(A16/'results/A16_coordinate_1bp_conclusions.tsv',details)

    bw=json.loads((A16/'results/A16_bigwig_manifest.json').read_text())
    bw_checks=[]
    for r in bw['files']:
        f=ROOT/r['path']
        with open(f,'rb') as h:magic=struct.unpack('<I',h.read(4))[0]
        bw_checks.append({'gsm':r['gsm'],'path':r['path'],'sha256_verified':sha(f)==r['sha256'],
            'valid_bigwig_magic':magic==0x888FFC26,'bytes_match':f.stat().st_size==r['bytes']})
    sig=pd.read_csv(A16/'results/A16_bigwig_promoter_signals.tsv',sep='\t')
    contrast=[]
    for (tf,region),sub in sig.groupby(['TF','region']):
        means=sub.set_index('label').mean_zero_filled
        for chip,inp in [('ChIP_rep1','Input_rep1'),('ChIP_rep2','Input_rep23'),('ChIP_rep3','Input_rep23')]:
            contrast.append({'TF':tf,'region':region,'ChIP':chip,'input':inp,'ChIP_mean':means[chip],
                'input_mean':means[inp],'descriptive_difference':means[chip]-means[inp],
                'interpretation':'coverage support only; no significance, peak, binding score, or causal inference'})
    save_tsv(A16/'results/A16_FOXC1_signal_input_comparison.tsv',contrast)
    q=pd.DataFrame(quality)
    save_json(A16/'results/A16_quality_summary.json',{'dataset_source_rows':len(d),'unique_paths':len(cache),
        'malformed_rows':int(q.drop_duplicates('path').malformed.sum()),
        'unknown_chromosomes':int(q.drop_duplicates('path').unknown_chromosome.sum()),
        'out_of_bounds':int(q.drop_duplicates('path').beyond_chromosome_end.sum()),
        'files_with_duplicate_intervals':int(q.drop_duplicates('path').exact_duplicate_intervals.gt(0).sum()),
        'all_peak_hashes_match':bool(q.sha256_verified.all()),'bigwig_checks':bw_checks,
        'signal_rows':len(sig),'all_signal_regions_quantified':len(sig)==5*17*2,
        'coverage_fraction_valid':bool(sig.coverage_fraction.between(0,1+1e-12).all()),
        'QC_unavailable':'FRiP/IDR/alignment metrics not recomputed from BED or coverage tracks; low peak count is only a flag',
        'checked_utc':stamp()})
    print(q.groupby('source')[['malformed','unknown_chromosome','beyond_chromosome_end','exact_duplicate_intervals']].sum().to_string())
    print('bigWig integrity:',all(all(x[k] for k in ['sha256_verified','valid_bigwig_magic','bytes_match']) for x in bw_checks))

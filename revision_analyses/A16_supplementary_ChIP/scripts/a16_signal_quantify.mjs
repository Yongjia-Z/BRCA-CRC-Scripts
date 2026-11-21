import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { createRequire } from 'node:module';
const root=path.resolve(process.argv[2]);
const a16=path.join(root,'revision_analyses/A16_supplementary_ChIP');
const require=createRequire(path.join(a16,'runtime/package.json'));
const {BigWig}=require('@gmod/bbi');
const lines=fs.readFileSync(path.join(a16,'results/A16_promoter_windows.tsv'),'utf8').replace(/^\uFEFF/,'').trim().split(/\r?\n/);
const keys=lines.shift().split('\t');
const promoters=lines.map(s=>Object.fromEntries(s.split('\t').map((v,i)=>[keys[i],v])));
const regions=promoters.flatMap(p=>[
 {...p,region:'promoter_TSS_pm2kb'},
 {...p,region:'neighborhood_TSS_pm100kb',start:Math.max(0,+p.tss_legacy-100000),end:+p.tss_legacy+100000},
]);
const samples=[
 ['GSM6595097_MDA231-ChIP_S2.bigWig','GSM6595097','ChIP_rep1'],
 ['GSM6595098_MDA231-ChIP-rep2_S1.bigWig','GSM6595098','ChIP_rep2'],
 ['GSM6595099_MDA231-ChIP-rep3_S5.bigWig','GSM6595099','ChIP_rep3'],
 ['GSM6595100_MDA231-Input_S1.bigWig','GSM6595100','Input_rep1'],
 ['GSM6595101_MDA231-Input-rep23_S9.bigWig','GSM6595101','Input_rep23'],
];
async function hash(file){const h=crypto.createHash('sha256');for await(const b of fs.createReadStream(file)){h.update(b)}return h.digest('hex')}
const rows=[],manifest=[];
for(const [name,gsm,label] of samples){
 const file=path.join(root,'datasets/Figure4_browser_tracks/data',name);
 if(!fs.existsSync(file)){manifest.push({gsm,label,status:'missing_local_file',path:file});continue;}
 const bw=new BigWig({path:file});const header=await bw.getHeader();
 manifest.push({gsm,label,status:'reused_read_only',path:path.relative(root,file),bytes:fs.statSync(file).size,sha256:await hash(file),genome:'hg38',source_url:`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm}`,type:'library-size-scaled coverage; not a significance track',bigwig_header_version:header.version});
 for(const p of regions){
  const start=+p.start,end=+p.end;const feats=await bw.getFeatures(p.chrom,start,end);
  let covered=0,integral=0,max=0;
  for(const f of feats){const n=Math.max(0,Math.min(end,f.end)-Math.max(start,f.start));covered+=n;integral+=n*f.score;max=Math.max(max,f.score);}
  rows.push({gsm,label,TF:p.TF,region:p.region,chrom:p.chrom,start,end,n_features:feats.length,covered_bases:covered,coverage_fraction:covered/(end-start),mean_zero_filled:integral/(end-start),mean_covered:covered?integral/covered:null,max_signal:max});
 }
 console.log(gsm,label,'quantified');
}
const fields=Object.keys(rows[0]||{});fs.writeFileSync(path.join(a16,'results/A16_bigwig_promoter_signals.tsv'),'\uFEFF'+fields.join('\t')+'\n'+rows.map(r=>fields.map(k=>r[k]??'').join('\t')).join('\n')+'\n');
fs.writeFileSync(path.join(a16,'results/A16_bigwig_manifest.json'),JSON.stringify({files:manifest,software:'@gmod/bbi',version:JSON.parse(fs.readFileSync(path.join(a16,'runtime/node_modules/@gmod/bbi/package.json'),'utf8')).version,note:'Signal-only supporting characterization from existing GSE213840 files; not new independent binding evidence and not used for scores.',checked_utc:new Date().toISOString()},null,2));

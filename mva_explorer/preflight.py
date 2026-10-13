"""Blind patient VCF preflight and deterministic real-data smoke subset."""
from collections import Counter, defaultdict
import argparse, csv, hashlib, json, os
from pathlib import Path
from datetime import datetime, timezone
import pysam

PRIMARY = tuple(map(str,range(1,23))) + ('X','Y','MT')
GRCH38_LENGTHS={'1':248956422,'2':242193529,'3':198295559,'4':190214555,'5':181538259,
'6':170805979,'7':159345973,'8':145138636,'9':138394717,'10':133797422,
'11':135086622,'12':133275309,'13':114364328,'14':107043718,'15':101991189,
'16':90338345,'17':83257441,'18':80373285,'19':58617616,'20':64444167,
'21':46709983,'22':50818468,'X':156040895,'Y':57227415,'MT':16569}

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def atomic_json(path,value):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2));tmp.replace(path)

def classify(ref,alts):
    if not alts:return 'NO_ALT'
    if any(a.startswith('<') or '[' in a or ']' in a or a=='*' for a in alts):return 'SV_OR_SYMBOLIC'
    if all(len(ref)==len(a)==1 for a in alts):return 'SNV'
    if any(len(ref)!=len(a) for a in alts):return 'INDEL'
    return 'MNV'

def index_check(path,vf):
    candidates=[Path(str(path)+'.tbi'),Path(str(path)+'.csi')]
    present=[str(p) for p in candidates if p.is_file()]
    if not present:return {'present':False,'files':[],'fetch_test':'NOT_TESTED_NO_INDEX'}
    try:
        contig=next(iter(vf.header.contigs),None)
        if contig: next(vf.fetch(contig,0,min(1000000,vf.header.contigs[contig].length or 1000000)),None)
        return {'present':True,'files':present,'fetch_test':'PASS'}
    except Exception as e:return {'present':True,'files':present,'fetch_test':'FAIL','error':str(e)}

def run_preflight(vcf_path,out,sample=None,per_contig=40):
    vcf_path=Path(vcf_path).resolve();out=Path(out).resolve();out.mkdir(parents=True,exist_ok=True)
    report_path=out/'01_preflight.json';qc_path=out/'02_variant_qc_summary.csv';subset_path=out/'patient_smoke_subset.vcf.gz'
    started=datetime.now(timezone.utc).isoformat();digest=sha256(vcf_path)
    counts=Counter();filters=Counter();genotypes=Counter();formats=Counter();subset_counts=Counter()
    chosen=defaultdict(list)
    with pysam.VariantFile(str(vcf_path)) as vf:
        samples=list(vf.header.samples)
        if sample is None:
            if len(samples)!=1:raise ValueError('Choose one sample explicitly.')
            sample=samples[0]
        if sample not in samples:raise ValueError('Sample not found in VCF.')
        index=index_check(vcf_path,vf)
        header_text=str(vf.header)
        references=[line for line in header_text.splitlines() if line.startswith('##reference=')]
        contig_lengths={c:h.length for c,h in vf.header.contigs.items()}
        canonical={c.removeprefix('chr'):n for c,n in contig_lengths.items() if c.removeprefix('chr') in GRCH38_LENGTHS}
        matches={c:(canonical.get(c)==n) for c,n in GRCH38_LENGTHS.items() if c in canonical}
        build='GRCh38' if len(matches)>=24 and all(matches.values()) else 'UNVERIFIED'
        annotation_fields=[x for x in ('CSQ','ANN','CLNSIG','gnomAD_AF','gnomADg_AF','gnomADe_AF') if x in vf.header.info]
        for i,rec in enumerate(vf,1):
            counts['records']+=1;typ=classify(rec.ref,rec.alts or ());counts['type_'+typ]+=1
            if len(rec.alts or ())>1:counts['multiallelic_records']+=1
            filt=';'.join(rec.filter) or '.';filters[filt]+=1
            call=rec.samples[sample];gt=call.get('GT') or ()
            carried=any(a is not None and a>0 for a in gt)
            if carried:counts['carried_records']+=1
            if not gt or any(a is None for a in gt):genotypes['missing_or_partial']+=1
            elif all(a==0 for a in gt):genotypes['hom_ref']+=1
            elif len(set(gt))==1:genotypes['hom_alt_or_haploid_alt']+=1
            else:genotypes['heterozygous_or_multiallelic']+=1
            for key in ('GT','DP','GQ','AD','PGT','PID','PS'):
                if key in call and call.get(key) is not None:formats[key+'_present']+=1
                else:formats[key+'_missing']+=1
            canon=rec.contig.removeprefix('chr');canon='MT' if canon=='M' else canon
            # Deterministic, genome-distributed smoke subset. No gene/disease selection.
            if carried and canon in PRIMARY and len(chosen[canon])<per_contig:
                chosen[canon].append(rec.copy());subset_counts[canon]+=1
        header=vf.header.copy()
    with pysam.VariantFile(str(subset_path),'wz',header=header) as dest:
        for contig in PRIMARY:
            for rec in chosen[contig]:dest.write(rec)
    try:pysam.tabix_index(str(subset_path),preset='vcf',force=True)
    except Exception:pass
    report={'status':'PREFLIGHT_PASS' if build=='GRCh38' and counts['records'] and counts['carried_records'] else 'PREFLIGHT_FAIL',
      'started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),'input':str(vcf_path),
      'input_bytes':vcf_path.stat().st_size,'input_sha256':digest,'samples':samples,'selected_sample':sample,
      'build_assessment':build,'build_basis':{'reference_headers':references,'canonical_contig_length_matches':matches},
      'index':index,'header':{'format_fields':sorted(header.formats),
      'info_fields':sorted(header.info),'preexisting_annotation_fields':annotation_fields},
      'counts':dict(counts),'filters':dict(filters),'genotypes':dict(genotypes),'format_availability':dict(formats),
      'smoke_subset':{'path':str(subset_path),'sha256':sha256(subset_path),'selection':'first carried records per primary contig in VCF order; no gene or disease selection','per_contig_limit':per_contig,'records':sum(subset_counts.values()),'by_contig':dict(subset_counts)},
      'phenotype_read':False,'gene_whitelist_used':False}
    atomic_json(report_path,report)
    with qc_path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['category','metric','value'])
        for category,data in [('counts',counts),('filters',filters),('genotypes',genotypes),('format',formats),('subset_contig',subset_counts)]:
            for key,value in sorted(data.items()):w.writerow([category,key,value])
    if report['status']!='PREFLIGHT_PASS':raise RuntimeError('Preflight failed; inspect 01_preflight.json.')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('vcf');p.add_argument('--out',required=True);p.add_argument('--sample');p.add_argument('--per-contig',type=int,default=40)
    a=p.parse_args();print(json.dumps(run_preflight(a.vcf,a.out,a.sample,a.per_contig),indent=2))

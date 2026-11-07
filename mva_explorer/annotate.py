"""Resumable reference validation, normalization and offline VEP annotation.

No patient phenotype access, gene whitelist, causal score or ACMG classification.
The reference is an annotation reference, NOT asserted to be the exact calling FASTA.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import pysam
from .core import sha256, allele_balance

VERSION = 'annotation-0.2.0'

def write_json(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2))
    temp.replace(path)

def timestamp():
    return datetime.now(timezone.utc).isoformat()

def variant_fingerprint(path):
    """Detect dropped, duplicated or changed records across annotation."""
    h=hashlib.sha256(); count=0
    with pysam.VariantFile(str(path)) as vf:
        for r in vf:
            count+=1
            h.update(json.dumps([r.contig,r.pos,r.ref,r.alts],separators=(',',':')).encode()+b'\n')
    return {'records':count,'ordered_locus_sha256':h.hexdigest()}

def partition(vcf, fasta, out, sample, progress=lambda **kw: None):
    counts = Counter()
    with pysam.VariantFile(str(vcf)) as src, pysam.FastaFile(str(fasta)) as ref:
        if sample not in src.header.samples:
            raise ValueError('Selected sample not in VCF.')
        header = src.header.copy()
        if 'END' not in header.info:
            header.info.add('END',1,'Integer','End position of the variant')
        for name, n, typ, description in [
            ('MVA_GT',1,'String','Original unsplit genotype for the selected sample'),
            ('MVA_SITE',1,'String','Original CHROM:POS before normalization'),
            ('MVA_AB','A','Float','Original alternate depth / sum of all original allelic depths'),
            ('MVA_REASON',1,'String','Reason retained outside current annotation scope')]:
            if name in header.info: raise ValueError('MVA tags already present; choose original input.')
            header.info.add(name,n,typ,description)
        targets = dict(zip(ref.references,ref.lengths))
        for name,length in targets.items():
            if name not in header.contigs: header.contigs.add(name,length=length)
        with pysam.VariantFile(str(out/'reference_compatible.vcf.gz'),'wz',header=header) as good, pysam.VariantFile(str(out/'unresolved.vcf.gz'),'wz',header=header) as unresolved:
            for i, rec in enumerate(src,1):
                counts['records_seen'] += 1
                if i % 50000 == 0: progress(records_scanned=i, counts=dict(counts))
                call = rec.samples[sample]
                gt = call.get('GT') or ()
                if not any(a is not None and a>0 for a in gt):
                    counts['no_called_alt_in_selected_sample'] += 1
                    continue
                counts['carried_records'] += 1
                original_contig = rec.contig
                normalized_contig = original_contig.removeprefix('chr')
                if normalized_contig == 'M': normalized_contig='MT'
                target = original_contig if original_contig in targets else normalized_contig
                reason = None
                if target not in targets:
                    reason='CONTIG_NOT_IN_ANNOTATION_REFERENCE'
                elif src.header.contigs[original_contig].length not in (None,targets[target]):
                    reason='CONTIG_LENGTH_MISMATCH'
                elif not rec.alts or any(not re.fullmatch('[ACGT]+',a,re.I) for a in rec.alts) or not re.fullmatch('[ACGT]+',rec.ref,re.I):
                    reason='SYMBOLIC_SV_OR_AMBIGUOUS_ALLELE'
                else:
                    observed = ref.fetch(target,rec.start,rec.start+len(rec.ref)).upper()
                    if observed != rec.ref.upper(): reason='REFERENCE_MISMATCH'
                copy=rec.copy(); copy.translate(header)
                copy.info['MVA_GT']=('|' if call.phased else '/').join('.' if a is None else str(a) for a in gt)
                copy.info['MVA_SITE']=original_contig+':'+str(rec.pos)
                copy.info['MVA_AB']=tuple(allele_balance(call.get('AD'),j) for j in range(1,len(rec.alts or ())+1))
                if reason:
                    copy.info['MVA_REASON']=reason
                    unresolved.write(copy); counts[reason]+=1
                else:
                    copy.contig=target
                    good.write(copy); counts['reference_compatible_records']+=1
    if not counts['reference_compatible_records']:
        raise ValueError('No carried alleles have a compatible annotation reference.')
    write_json(out/'reference_audit.json',dict(counts))
    return dict(counts)

def population_max(annotations):
    values=[]
    for a in annotations:
        for k,v in a.items():
            if k.startswith(('gnomADe_','gnomADg_')) and k.endswith('_AF'):
                for token in re.split('[&/]',v):
                    try:
                        n=float(token)
                        if 0<=n<=1: values.append(n)
                    except (ValueError,TypeError): pass
    return max(values) if values else None

def extract(annotated, out, sample, progress=lambda **kw:None):
    counts=Counter(); genes=defaultdict(Counter)
    fields=['variant_key','chrom','pos','ref','alt','gene_id','symbol','consequence','impact',
            'transcripts','mane_transcripts','mane_plus_clinical_transcripts','HGVSc','HGVSp','gnomAD_max_observed_AF',
            'population_status','cached_CLIN_SIG','normalized_GT','original_GT','original_site',
            'DP','GQ','original_allele_balance','FILTER','PID','PGT','PS','review_bucket',
            'ACMG_classification','case_causality']
    handles={}; writers={}
    try:
        for name in ('review_candidates.csv.gz','other_rare_or_unknown.csv.gz'):
            handles[name]=gzip.open(out/name,'wt',newline='')
            writers[name]=csv.DictWriter(handles[name],fieldnames=fields);writers[name].writeheader()
        with pysam.VariantFile(str(annotated)) as vf:
            if 'CSQ' not in vf.header.info: raise ValueError('VEP did not produce CSQ annotations.')
            keys=vf.header.info['CSQ'].description.split('Format: ',1)[-1].split('|')
            if 'Gene' not in keys or 'Consequence' not in keys or 'gnomADg_AF' not in keys:
                raise ValueError('Missing required VEP fields; inspect VEP options and cache.')
            for i,rec in enumerate(vf,1):
                if i%50000==0:progress(annotated_records_exported=i,counts=dict(counts))
                if len(rec.alts or ())!=1: raise ValueError('Expected split biallelic VEP records.')
                call=rec.samples[sample];gt=call.get('GT') or ()
                def info_value(key):return rec.info.get(key) if key in rec.header.info else None
                if 1 not in gt: continue
                counts['carried_annotated_alleles']+=1
                annotations=[dict(zip(keys,item.split('|'))) for item in rec.info.get('CSQ',())]
                af=population_max(annotations)
                if af is None: counts['population_frequency_unknown']+=1
                # This is a broad export convenience threshold, not PM2/BA1 or a diagnosis filter.
                if af is not None and af>.01:
                    counts['observed_AF_above_export_threshold_retained_in_VCF']+=1
                    continue
                groups=defaultdict(list)
                for a in annotations:groups[a.get('Gene') or 'UNASSIGNED'].append(a)
                if not groups:groups['UNASSIGNED']=[];counts['missing_CSQ_retained']+=1
                for gene, entries in groups.items():
                    def union(key):return '&'.join(sorted({a.get(key,'') for a in entries}-{''}))
                    impact=union('IMPACT'); consequence=union('Consequence')
                    clinical=union('CLIN_SIG')
                    is_review=bool({'HIGH','MODERATE'} & set(impact.split('&'))) or 'splice' in consequence or 'pathogenic' in clinical.lower()
                    bucket='coding_splice_or_cached_clinical_label' if is_review else 'other_non_coding_or_unassigned'
                    row=dict(variant_key=f'GRCh38:{rec.contig}:{rec.pos}:{rec.ref}:{rec.alts[0]}',chrom=rec.contig,pos=rec.pos,ref=rec.ref,alt=rec.alts[0],gene_id=gene,symbol=union('SYMBOL'),consequence=consequence,impact=impact,transcripts=union('Feature'),mane_transcripts=union('MANE_SELECT'),mane_plus_clinical_transcripts=union('MANE_PLUS_CLINICAL'),HGVSc=union('HGVSc'),HGVSp=union('HGVSp'),gnomAD_max_observed_AF=af,population_status='UNKNOWN_CACHE_NOT_COMPLETE_GNOMAD' if af is None else 'OBSERVED_IN_VEP_CACHE',cached_CLIN_SIG=clinical,normalized_GT=('|' if call.phased else '/').join('.' if a is None else str(a) for a in gt),original_GT=info_value('MVA_GT'),original_site=info_value('MVA_SITE'),DP=call.get('DP'),GQ=call.get('GQ'),original_allele_balance=(info_value('MVA_AB') or (None,))[0],FILTER=';'.join(rec.filter) or '.',PID=call.get('PID'),PGT=call.get('PGT'),PS=call.get('PS'),review_bucket=bucket,ACMG_classification='NOT_ASSESSED',case_causality='NOT_ASSESSED')
                    name='review_candidates.csv.gz' if is_review else 'other_rare_or_unknown.csv.gz'
                    writers[name].writerow(row); counts[name]+=1
                    if is_review:
                        genes[(gene,row['symbol'])]['allele_gene_rows']+=1
                        if tuple(gt)==(1,1):genes[(gene,row['symbol'])]['homozygous_alt_rows']+=1
                        if any(x is None for x in gt):genes[(gene,row['symbol'])]['partial_genotype_rows']+=1
        with open(out/'gene_inventory.csv','w',newline='') as f:
            w=csv.writer(f);w.writerow(['gene_id','symbol','allele_gene_rows','homozygous_alt_rows','partial_genotype_rows'])
            for (gene,symbol),c in sorted(genes.items()):w.writerow([gene,symbol,c['allele_gene_rows'],c['homozygous_alt_rows'],c['partial_genotype_rows']])
    finally:
        for f in handles.values():f.close()
    write_json(out/'annotation_summary.json',dict(counts))
    return dict(counts)

def run(args):
    base=Path(__file__).resolve().parents[1]
    out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    lock=open(out/'.lock','w')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise ValueError('This output folder already has an active annotation run.')
    fasta=base/'resources/reference/GRCh38_primary.fa'
    cache=base/'resources/vep115'
    info=cache/'homo_sapiens/115_GRCh38/info.txt'
    vcf=Path(args.vcf).resolve()
    if not all(p.is_file() for p in (vcf,fasta,info)):
        raise ValueError('Input/resources missing. Run bash setup_annotation.sh first.')
    prefix=base/'.tools/vep115'
    conda=shutil.which('conda')
    if not conda:raise ValueError('Conda is not on PATH.')
    state_path=out/'annotation_status.json'
    old=json.loads(state_path.read_text()) if state_path.exists() else {}
    state=dict(version=VERSION,status='STARTING',phase='HASHING',updated=timestamp(),phenotype_read=False,
               final_ranking='NOT_PERFORMED',warnings=[
               'Annotation reference is Ensembl GRCh38 primary assembly; original masked/decoy calling FASTA not supplied.',
               'Contig and REF-incompatible records and structural alleles retained in unresolved.vcf.gz.',
               'Missing gnomAD cache frequencies are UNKNOWN, not zero or proven absence.',
               'Export cutoff AF<=0.01 is only a review convenience; all annotated calls remain in annotated.vcf.gz.',
               'Cached clinical labels lack full assertion review; no ACMG criteria assigned.'])
    def progress(**kwargs):
        state.update(kwargs);state['updated']=timestamp();write_json(state_path,state)
        print(state['phase'],kwargs,flush=True)
    def command(argv):
        state.setdefault('commands',[]).append(argv)
        progress(last_command=argv)
        with open(out/'commands.log','a') as log:
            log.write(json.dumps(argv)+'\n');log.flush()
            subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,check=True)
    def exe(name):
        # VEP lives in its pinned Perl environment. The host bcftools/samtools
        # are deliberately kept outside it to avoid incompatible zlib pins.
        if name in ('bcftools','samtools'):
            path=shutil.which(name)
            if not path:raise ValueError(name+' is not on PATH.')
            return [path]
        return [conda,'run','--no-capture-output','-p',str(prefix),name]
    try:
        identity={'input_sha256':sha256(vcf),'reference_sha256':sha256(fasta),'cache_info_sha256':sha256(info),'sample':args.sample,'version':VERSION}
        if old.get('identity') and old['identity']!=identity:raise ValueError('Inputs/resources changed; choose a new output folder.')
        state['identity']=identity
        if (base/'resources/annotation_resources.sha256').exists():state['resource_checksums']=(base/'resources/annotation_resources.sha256').read_text()
        state['cache_metadata']=info.read_text()
        stages=[('partition',[out/'reference_compatible.vcf.gz',out/'unresolved.vcf.gz',out/'reference_audit.json']),
                ('normalize',[out/'normalized.vcf.gz']),('vep',[out/'annotated.vcf.gz',out/'record_conservation.json']),
                ('export',[out/'review_candidates.csv.gz',out/'other_rare_or_unknown.csv.gz',out/'gene_inventory.csv',out/'annotation_summary.json'])]
        for name,outputs in stages:
            marker=out/(name+'.done.json')
            # Completed-stage artifacts are hashed; a stale/partial output is not reused.
            if marker.exists():
                saved=json.loads(marker.read_text())
                if saved.get('identity')==identity and all(p.exists() and saved['outputs'].get(p.name)==sha256(p) for p in outputs):
                    progress(status='RUNNING',phase=name.upper()+'_REUSED');continue
            progress(status='RUNNING',phase=name.upper())
            if name=='partition':
                partition(vcf,fasta,out,args.sample,progress)
            elif name=='normalize':
                command(exe('bcftools')+['norm','-f',str(fasta),'-c','e','-m','-any','--old-rec-tag','MVA_NORM_ORIG','--keep-sum','AD','-Oz','-o',str(out/'normalized.unsorted.vcf.gz'),str(out/'reference_compatible.vcf.gz')])
                command(exe('bcftools')+['sort','-Oz','-o',str(outputs[0]),str(out/'normalized.unsorted.vcf.gz')])
            elif name=='vep':
                command(exe('vep')+['--offline','--cache','--cache_version','115','--dir_cache',str(cache),'--assembly','GRCh38','--species','homo_sapiens','--fasta',str(fasta),'--format','vcf','--vcf','--input_file',str(out/'normalized.vcf.gz'),'--output_file',str(outputs[0]),'--compress_output','bgzip','--force_overwrite','--symbol','--canonical','--mane','--biotype','--hgvs','--numbers','--variant_class','--allele_number','--check_existing','--af_gnomade','--af_gnomadg','--regulatory','--fork',str(args.threads),'--buffer_size','1000','--stats_file',str(out/'vep_summary.html')])
                before=variant_fingerprint(out/'normalized.vcf.gz');after=variant_fingerprint(outputs[0])
                write_json(out/'record_conservation.json',{'before':before,'after':after,'equal':before==after})
                if before!=after:raise ValueError('VEP changed/dropped/reordered records; inspect record_conservation.json before exporting.')
                command(exe('bcftools')+['index','-f',str(outputs[0])])
            elif name=='export':extract(out/'annotated.vcf.gz',out,args.sample,progress)
            write_json(marker,{'identity':identity,'outputs':{p.name:sha256(p) for p in outputs}})
            # Recompute downstream stages after any stage was rerun.
            for later,_ in stages[[n for n,_ in stages].index(name)+1:]:
                (out/(later+'.done.json')).unlink(missing_ok=True)
        progress(status='ANNOTATION_COMPLETE',phase='COMPLETE',reference_audit=json.loads((out/'reference_audit.json').read_text()),annotation_summary=json.loads((out/'annotation_summary.json').read_text()))
    except Exception as exc:
        progress(status='FAILED',error=str(exc));raise
    finally:
        fcntl.flock(lock,fcntl.LOCK_UN);lock.close()

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('vcf');parser.add_argument('--sample',required=True)
    parser.add_argument('--out',default='data/analysis/genome_only')
    parser.add_argument('--threads',type=int,default=2)
    a=parser.parse_args()
    if not 1<=a.threads<=16:parser.error('threads must be 1–16')
    run(a)

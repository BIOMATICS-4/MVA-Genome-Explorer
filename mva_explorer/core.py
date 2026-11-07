"""Genome-only QC; no clinical interpretation or phenotype access."""
import csv
import hashlib
import json
import sqlite3
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone

VERSION = '0.1.0'

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def allele_balance(ad, index):
    if not ad or any(x is None or x < 0 for x in ad) or index >= len(ad):
        return None
    total = sum(ad)
    return ad[index] / total if total else None

def genotype_type(gt, alt_index):
    if not gt or any(x is None for x in gt):
        return 'missing_or_partial'
    n = gt.count(alt_index)
    if n == 0:
        return 'not_carried'
    if len(gt) == 1:
        return 'hemizygous_or_haploid'
    if n == len(gt):
        return 'homozygous_alt' if len(gt) == 2 else 'all_alt_polyploid'
    return 'heterozygous'

# Reference lengths provide build evidence without assuming that a file can
# be safely annotated against the matching FASTA.
_BUILD_LENGTHS = {
    'GRCh37': {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276,
        '5': 180915260, '6': 171115067, '7': 159138663, '8': 146364022,
        '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895,
        '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753,
        '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
        '21': 48129895, '22': 51304566, 'X': 155270560, 'Y': 59373566},
    'GRCh38': {'1': 248956422, '2': 242193529, '3': 198295559, '4': 190214555,
        '5': 181538259, '6': 170805979, '7': 159345973, '8': 145138636,
        '9': 138394717, '10': 133797422, '11': 135086622, '12': 133275309,
        '13': 114364328, '14': 107043718, '15': 101991189, '16': 90338345,
        '17': 83257441, '18': 80373285, '19': 58617616, '20': 64444167,
        '21': 46709983, '22': 50818468, 'X': 156040895, 'Y': 57227415}
}

def detect_genome_build(header):
    """Conservative QC-only assembly inference; never authorizes annotation."""
    import re
    references = [str(record) for record in header.records if record.key == 'reference']
    hint_text = ' '.join(references)
    hints = set()
    if re.search(r'(?i)grch37|hg19', hint_text): hints.add('GRCh37')
    if re.search(r'(?i)grch38|hg38', hint_text): hints.add('GRCh38')
    observed = {}
    for name, contig in header.contigs.items():
        key = name[3:] if name.startswith('chr') else name
        if key in _BUILD_LENGTHS['GRCh37'] and contig.length is not None:
            observed[key] = contig.length
    matches = {b: sum(observed[k] == length for k, length in table.items() if k in observed)
               for b, table in _BUILD_LENGTHS.items()}
    conflicts = {b: sum(observed[k] != length for k, length in table.items() if k in observed)
                 for b, table in _BUILD_LENGTHS.items()}
    basis = {'reference_headers': references, 'informative_contigs': len(observed),
             'matching_contig_lengths': matches, 'conflicting_contig_lengths': conflicts}
    plausible = [b for b in _BUILD_LENGTHS if matches[b] >= 4 and conflicts[b] == 0]
    if len(hints) > 1:
        basis['method'] = 'CONFLICTING_REFERENCE_HEADERS'
        return 'UNVERIFIED', basis
    if len(plausible) == 1:
        b = plausible[0]
        if hints and b not in hints:
            basis['method'] = 'REFERENCE_CONTIG_CONFLICT'
            return 'UNVERIFIED', basis
        basis['method'] = 'REFERENCE_AND_CONTIG_LENGTHS' if hints else 'CONTIG_LENGTHS'
        return b, basis
    if len(hints) == 1:
        b = next(iter(hints))
        if conflicts[b] == 0:
            basis['method'] = 'REFERENCE_HEADER_ONLY' if not observed else 'REFERENCE_WITH_PARTIAL_CONTIG_LENGTHS'
            return b, basis
    basis['method'] = 'INSUFFICIENT_OR_CONFLICTING_EVIDENCE'
    return 'UNVERIFIED', basis

def qc(vcf_path, out_dir, sample=None):
    import pysam
    vcf_path, out_dir = Path(vcf_path), Path(out_dir)
    if not vcf_path.is_file() or vcf_path.name.endswith(('.tbi', '.zip')):
        raise ValueError('A VCF/BCF file is required; an index is insufficient.')
    out_dir.mkdir(parents=True, exist_ok=True)
    if (out_dir/'run_manifest.json').exists():
        raise ValueError('Use a new run directory; existing run is preserved.')

    # Write a manifest immediately so the UI can show a real stage/progress state
    # while the potentially large input is being checksummed.
    manifest_path = out_dir/'run_manifest.json'
    input_size = vcf_path.stat().st_size
    manifest = dict(software_version=VERSION, pysam_version=pysam.__version__,
                    input_file=vcf_path.name, input_size_bytes=input_size,
                    input_sha256=None, sample=sample, samples=[],
                    started_utc=datetime.now(timezone.utc).isoformat(),
                    status='QC_HASHING_INPUT', phenotype_read=False,
                    genome_build='UNVERIFIED', normalization='NOT_PERFORMED',
                    ranking='BLOCKED', resources={},
                    progress={'stage':'Hashing input', 'percent':0.0,
                              'records_scanned':0, 'carried_alleles':0},
                    blockers=['Reference/build verification', 'Normalization',
                              'Genome-wide annotation', 'Reviewed evidence and inheritance engine'])
    def save():
        temp = out_dir/'run_manifest.tmp'
        temp.write_text(json.dumps(manifest, indent=2))
        temp.replace(manifest_path)
    save()

    digest = sha256(vcf_path)
    vf = pysam.VariantFile(str(vcf_path))
    samples = list(vf.header.samples)
    build, build_basis = detect_genome_build(vf.header)
    if sample is None:
        if len(samples) != 1:
            raise ValueError('Choose a sample explicitly for a multisample/site-only VCF.')
        sample = samples[0]
    if sample not in samples:
        raise ValueError('Sample not found in VCF.')
    manifest.update(input_sha256=digest, sample=sample, samples=samples,
                    genome_build=build, genome_build_basis=build_basis,
                    status='QC_RUNNING',
                    progress={'stage':'Scanning variants', 'percent':0.0,
                              'records_scanned':0, 'carried_alleles':0})
    save()
    counts = Counter()
    db = sqlite3.connect(out_dir/'variants.sqlite')
    db.execute('CREATE TABLE alleles (record_id INTEGER, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, alt_index INTEGER, gt TEXT, genotype TEXT, dp INTEGER, gq REAL, ab REAL, qual REAL, filter TEXT, carried INTEGER)')
    try:
        for record_id, rec in enumerate(vf, 1):
            counts['records'] += 1
            filters = list(rec.filter)
            filt = ';'.join(filters) or '.'
            counts['filter_'+filt] += 1
            if 'CSQ' in rec.info: counts['records_with_CSQ'] += 1
            if 'ANN' in rec.info: counts['records_with_ANN'] += 1
            call = rec.samples[sample]
            gt = call.get('GT') or ()
            sep = '|' if call.phased else '/'
            gt_text = sep.join('.' if x is None else str(x) for x in gt)
            dp, gq, ad = call.get('DP'), call.get('GQ'), call.get('AD')
            if dp is None: counts['missing_DP'] += 1
            if gq is None: counts['missing_GQ'] += 1
            if ad is None: counts['missing_AD'] += 1
            alts = rec.alts or ()
            if len(alts) > 1: counts['multiallelic_records'] += 1
            for i, alt in enumerate(alts, 1):
                kind = genotype_type(gt, i)
                counts['genotype_'+kind] += 1
                sv = alt.startswith('<') or '[' in alt or ']' in alt or alt == '*'
                typ = 'symbolic_or_SV' if sv else ('SNV' if len(rec.ref)==len(alt)==1 else 'indel_or_MNV')
                counts[typ] += 1
                carried = int(i in gt)
                counts['carried_alleles'] += carried
                db.execute('INSERT INTO alleles VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                           (record_id, rec.contig, rec.pos, rec.ref, alt, i, gt_text,
                            kind, dp, gq, allele_balance(ad, i), rec.qual, filt, carried))
            if record_id % 10000 == 0:
                db.commit()
            if record_id % 50000 == 0:
                # htslib virtual offsets encode the compressed byte position in the
                # high bits for bgzip/BCF. For plain VCF, tell() is a byte offset.
                try:
                    raw_pos = int(vf.tell())
                    compressed = vcf_path.name.endswith(('.gz', '.bcf'))
                    byte_pos = (raw_pos >> 16) if compressed else raw_pos
                    pct = min(99.5, max(0.0, (byte_pos / input_size * 100.0) if input_size else 0.0))
                except Exception:
                    pct = manifest.get('progress', {}).get('percent', 0.0)
                manifest['progress'] = {
                    'stage':'Scanning variants', 'percent':round(pct, 1),
                    'records_scanned':record_id,
                    'carried_alleles':int(counts.get('carried_alleles', 0)),
                }
                save()
        db.commit()
        db.execute('CREATE INDEX allele_locus ON alleles(chrom,pos)')
        db.commit()
        with open(out_dir/'02_variant_qc_summary.csv', 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['metric','value']); w.writerows(sorted(counts.items()))
        manifest.update(status='QC_COMPLETE_ANALYSIS_BLOCKED', counts=dict(counts),
                        progress={'stage':'Technical QC complete', 'percent':100.0,
                                  'records_scanned':int(counts.get('records', 0)),
                                  'carried_alleles':int(counts.get('carried_alleles', 0))},
                        header_reference=[str(x) for x in vf.header.records if x.key=='reference'],
                        annotations_in_header=[x for x in ('CSQ','ANN') if x in vf.header.info],
                        index_present=Path(str(vcf_path)+'.tbi').exists() or Path(str(vcf_path)+'.csi').exists())
        save()
    except Exception as e:
        db.rollback(); manifest.update(status='QC_FAILED', error=str(e),
                                      progress={**manifest.get('progress', {}), 'stage':'QC failed'}); save(); raise
    finally:
        db.close(); vf.close()
    return manifest

def normalize(vcf, fasta, output):
    """HTSlib-backed bcftools normalization with strict REF checking."""
    import pysam.bcftools
    if Path(output).exists(): raise ValueError('Output exists.')
    if not Path(str(fasta)+'.fai').exists(): raise ValueError('Indexed reference FASTA required.')
    pysam.bcftools.norm('-f', str(fasta), '-c', 'e', '-m', '-any', '-Oz',
                        '-o', str(output), str(vcf), catch_stdout=False)
    return {'input_sha256':sha256(vcf), 'reference_sha256':sha256(fasta),
            'output_sha256':sha256(output), 'method':'bcftools norm -c e -m -any',
            'note':'Normalization is not annotation or clinical classification.'}

def audit_evidence(rows):
    """Conservative review queue. Does not assign ACMG classifications."""
    seen = set(); result = []
    for raw in rows:
        row = dict(raw); reasons = []
        key = (row.get('variant_key'), row.get('disease_id'), row.get('dependency_group'))
        if not all(key): reasons.append('missing_variant_disease_or_dependency')
        if key in seen: reasons.append('dependent_evidence_requires_joint_review')
        seen.add(key)
        if row.get('criterion') in ('PP5','BP6'): reasons.append('deprecated_criterion')
        if not row.get('source') or not row.get('observation_id'): reasons.append('missing_provenance')
        if row.get('criterion') == 'PM2' and row.get('population_status') != 'queried_callable_absent':
            reasons.append('absence_not_established')
        if row.get('criterion') == 'PVS1':
            for field in ('lof_mechanism','transcript_relevant','nmd_review','rescue_review'):
                if row.get(field) != 'reviewed_supported': reasons.append('PVS1_'+field+'_unresolved')
        if row.get('criterion') == 'PM3' and row.get('phase') != 'confirmed_trans':
            reasons.append('PM3_requires_phase_strength_review')
        row['audit_flags'] = ';'.join(reasons)
        row['status'] = 'REVIEW_REQUIRED' if reasons else 'PROVENANCE_CHECKS_PASSED_NOT_CLASSIFIED'
        result.append(row)
    return result

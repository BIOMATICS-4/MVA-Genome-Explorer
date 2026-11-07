"""Frozen inputs, G0 exports, and immutable layer boundaries."""
from __future__ import annotations

import argparse
from collections import defaultdict
from contextlib import contextmanager
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys

from .inheritance import ScoreConfig, Variant, evaluate_gene, number

ROOT = Path(__file__).resolve().parents[1]
PROTECTED = ROOT / 'data/analysis/genome_only'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def write_json(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)


def write_csv(path, rows, fields=None):
    with Path(path).open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields or (list(rows[0]) if rows else ['gene']))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list)) else v
                             for k, v in row.items()})


def annotation_processes(proc_root=Path('/proc')):
    """Process metadata only; never opens VCFs or sends signals."""
    if not proc_root.is_dir():
        raise RuntimeError('Cannot verify annotation processes on this platform')
    found = []
    for entry in proc_root.iterdir():
        if not entry.name.isdigit() or int(entry.name) == os.getpid():
            continue
        try:
            args = (entry / 'cmdline').read_bytes().split(b'\0')
        except FileNotFoundError:
            continue
        except PermissionError as exc:
            raise RuntimeError('Cannot verify process completion') from exc
        names = [Path(a.decode(errors='replace')).name.lower() for a in args if a]
        if any(n in ('vep', 'vep.pl', 'bgzip') for n in names):
            found.append(int(entry.name))
    return found


def assert_snapshot_path(path):
    resolved = Path(path).resolve()
    # Deny source paths BEFORE stat, hashing, parsing, decompression or indexing.
    if resolved.parent == PROTECTED.resolve() and resolved.name in ('annotated.vcf.gz', 'normalized.vcf.gz'):
        raise ValueError('Active annotation/normalized VCF is forbidden as a snapshot input')
    if resolved.name == 'normalized.vcf.gz':
        raise ValueError('Ranking does not consume normalized VCFs')
    return resolved


@contextmanager
def completed_annotation(path):
    """Fail closed until all writers exited and annotation completion is recorded."""
    path = Path(path).resolve()
    if path.name != 'annotated.vcf.gz':
        raise ValueError('Completed input must be an annotated.vcf.gz from the annotation stage')
    active = annotation_processes()
    if active:
        raise RuntimeError(f'Annotation writers still active (PIDs {active}); VCF not opened')
    folder = path.parent
    with (folder / '.lock').open('rb') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Annotation stage still holds its lock; VCF not opened') from exc
        state = json.loads((folder / 'annotation_status.json').read_text())
        conservation = json.loads((folder / 'record_conservation.json').read_text())
        marker = json.loads((folder / 'vep.done.json').read_text())
        if state.get('status') != 'ANNOTATION_COMPLETE' or conservation.get('equal') is not True:
            raise RuntimeError('Completed annotation and record-conservation evidence required; VCF not opened')
        if annotation_processes():
            raise RuntimeError('Annotation writers started; VCF not opened')
        # Only here is opening the completed VCF permitted.
        with gzip.open(path, 'rb') as f:
            for _ in iter(lambda: f.read(4 * 1024 * 1024), b''):
                pass
        checksum = digest(path)
        if marker.get('outputs', {}).get(path.name) != checksum:
            raise ValueError('Completed VCF differs from its annotation-stage checksum')
        import pysam
        import pysam.bcftools
        with pysam.VariantFile(str(path)) as checked:
            checked.check_truncation()
        if not Path(str(path) + '.csi').exists() and not Path(str(path) + '.tbi').exists():
            pysam.bcftools.index(str(path), catch_stdout=False)
        yield checksum, state


def from_row(r):
    def get(*keys, default=''):
        return next((r[k] for k in keys if r.get(k) not in (None, '', '.')), default)
    def obj(key):
        value = r.get(key, {})
        return json.loads(value) if isinstance(value, str) and value else value or {}
    gene = get('gene', 'gene_id', 'symbol')
    if gene in ('', 'UNASSIGNED'):
        # Unmapped loci cannot accidentally form pairs with other unmapped loci.
        gene = 'UNASSIGNED:' + get('variant_key', default=f"{r.get('chrom')}:{r.get('pos')}:{r.get('ref')}:{r.get('alt')}")
    relevance = get('transcript_relevance', default='UNKNOWN')
    # Old union-only exports cannot prove the severe consequence affects MANE.
    return Variant(gene=gene, symbol=get('symbol'), chrom=get('chrom'), pos=int(r['pos']),
                   ref=r['ref'], alt=r['alt'], gt=get('gt', 'GT', 'normalized_GT'),
                   build=get('build', default='GRCh38'), consequence=get('consequence'),
                   transcript_relevance=relevance, af=number(get('af', 'gnomAD_max_observed_AF')),
                   population_status=get('population_status', default='UNKNOWN'),
                   population_source=get('population_source'),
                   complete_callable_resource=str(r.get('complete_callable_resource', '')).lower() == 'true',
                   dp=number(get('dp', 'DP')), gq=number(get('gq', 'GQ')),
                   ab=number(get('ab', 'AB', 'allele_balance', 'original_allele_balance')),
                   filter=get('filter', 'FILTER', default='.'), ps=get('ps', 'PS'),
                   pid=get('pid', 'PID'), pgt=get('pgt', 'PGT'),
                   variant_class=get('variant_class', default='SMALL_VARIANT'),
                   end=int(r['end']) if r.get('end') not in (None, '', '.') else None,
                   copy_number=number(get('copy_number')), expected_copy_number=number(get('expected_copy_number')),
                   structural_qc=get('structural_qc', default='UNKNOWN'),
                   transcripts=get('transcripts'), selected_transcript=get('selected_transcript'),
                   transcript_biotype=get('transcript_biotype'),
                   transcript_selection=get('transcript_selection'),
                   alternate_transcript_effects=obj('alternate_transcript_effects'),
                   cached_clinical=get('cached_CLIN_SIG', 'clinical'),
                   predictors=obj('predictors'), computational_review=obj('computational_review'), reviewed_acmg=obj('reviewed_acmg'),
                   parents=obj('parents'))


def csv_variants(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as f:
        yield from (from_row(r) for r in csv.DictReader(f))


_IMPACT_ORDER = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0, '': 0}


def _consequence_priority(value):
    terms = set(re.split('[&,]', value or ''))
    if terms & {'transcript_ablation', 'stop_gained', 'frameshift_variant',
                'splice_acceptor_variant', 'splice_donor_variant'}:
        return 4
    if terms & {'missense_variant', 'inframe_insertion', 'inframe_deletion',
                'protein_altering_variant', 'start_lost', 'stop_lost'}:
        return 3
    if terms & {'splice_region_variant', 'splice_donor_region_variant',
                'splice_polypyrimidine_tract_variant', 'splice_donor_5th_base_variant'}:
        return 2
    if terms & {'regulatory_region_ablation', 'regulatory_region_amplification',
                'TFBS_ablation'}:
        return 1
    return 0


def select_transcript_entry(entries):
    """Choose one transcript-consistent VEP annotation for a gene/variant.

    MANE Select and MANE Plus Clinical are the highest-priority transcript set.
    Within that set, molecular impact determines the selected entry. Canonical
    transcript is next, then the best remaining transcript. Consequence,
    biotype and transcript relevance always come from the same entry.
    """
    if not entries:
        return {}, 'NO_TRANSCRIPT'

    ranked = []
    for e in entries:
        mane_select = bool(e.get('MANE_SELECT'))
        mane_plus = bool(e.get('MANE_PLUS_CLINICAL'))
        canonical = str(e.get('CANONICAL', '')).upper() in ('YES', '1')
        transcript_tier = 3 if (mane_select or mane_plus) else 2 if canonical else 1
        label = ('MANE_SELECT' if mane_select else
                 'MANE_PLUS_CLINICAL' if mane_plus else
                 'VEP_CANONICAL' if canonical else 'OTHER_TRANSCRIPT')
        ranked.append((
            transcript_tier,
            _IMPACT_ORDER.get(e.get('IMPACT', ''), 0),
            _consequence_priority(e.get('Consequence', '')),
            1 if mane_select else 0,
            e.get('Feature', ''),
            label,
            e,
        ))
    ranked.sort(key=lambda x: (-x[0], -x[1], -x[2], -x[3], x[4]))
    return ranked[0][-1], ranked[0][-2]


def vcf_variants(path, sample=None, build='GRCh38', counts=None, threads=1, fast=True, progress=None):
    """Read ONLY after the caller's completion/frozen-input gate.

    Full-genome stream, normalized biallelic records; retain all coding/splice,
    explicit regulatory effects and structural candidates regardless of AF.
    Plain intergenic/intronic annotations are counted rather than paired.
    """
    import pysam
    from .annotate import population_max
    counts = counts if counts is not None else defaultdict(int)
    with pysam.VariantFile(str(path), threads=threads) as vf:
        samples = list(vf.header.samples)
        if sample is None and len(samples) == 1:
            sample = samples[0]
        if sample not in samples:
            raise ValueError('Explicit valid sample required')
        if 'CSQ' not in vf.header.info:
            raise ValueError('VEP CSQ annotation required')
        keys = vf.header.info['CSQ'].description.split('Format: ', 1)[-1].split('|')
        if not {'Gene', 'Consequence', 'IMPACT'} <= set(keys):
            raise ValueError('Required VEP CSQ fields missing')
        gene_index = keys.index('Gene')
        for rec in vf:
            counts['records_scanned'] += 1
            if progress and counts['records_scanned'] % 100000 == 0:
                progress(dict(counts))
            if len(rec.alts or ()) != 1:
                raise ValueError('Normalize and split multiallelic records before annotation; ranking never changes inputs')
            call = rec.samples[sample]
            gt = call.get('GT') or ()
            if 1 not in gt:
                continue
            raw_annotations = rec.info.get('CSQ', ())
            alt = rec.alts[0]
            sv = alt.startswith('<') or '[' in alt or ']' in alt
            # Lossless prefilter: every existing candidate predicate implies one
            # of these tokens. Detailed AF/transcript parsing is unnecessary for
            # annotations that cannot yield a candidate under the unchanged rules.
            signals = ('HIGH', 'MODERATE', 'splice', 'regulatory_region_ablation',
                       'regulatory_region_amplification', 'TFBS_ablation')
            if fast and raw_annotations and not sv and not any(
                    any(token in item for token in signals) for item in raw_annotations):
                gene_ids = set()
                for item in raw_annotations:
                    fields = item.split('|', gene_index + 1)
                    gene_ids.add(fields[gene_index] if len(fields) > gene_index and fields[gene_index] else 'UNASSIGNED')
                counts['noncandidate_allele_gene_annotations'] += len(gene_ids)
                continue
            groups = defaultdict(list)
            annotations = [dict(zip(keys, item.split('|'))) for item in raw_annotations]
            for a in annotations:
                groups[a.get('Gene') or 'UNASSIGNED'].append(a)
            if not groups:
                groups['UNASSIGNED'] = []
            af = population_max(annotations)
            alt = rec.alts[0]
            sv = alt.startswith('<') or '[' in alt or ']' in alt
            ad = call.get('AD')
            ab = None
            if 'MVA_AB' in rec.header.info:
                raw = rec.info.get('MVA_AB')
                if raw:
                    ab = raw[0] if isinstance(raw, tuple) else raw
            if ab is None and ad and all(x is not None and x >= 0 for x in ad) and sum(ad) > 0:
                ab = ad[1] / sum(ad)
            for gene, entries in groups.items():
                def union(key):
                    return '&'.join(sorted({x.get(key, '') for x in entries} - {''}))
                selected, selection = select_transcript_entry(entries)
                cons = selected.get('Consequence', '') if selected else ''
                all_cons = union('Consequence')
                candidate = (sv or not entries or any(x.get('IMPACT') in ('HIGH', 'MODERATE') for x in entries)
                             or 'splice' in all_cons or any(t in all_cons for t in ('regulatory_region_ablation', 'regulatory_region_amplification', 'TFBS_ablation')))
                if not candidate:
                    counts['noncandidate_allele_gene_annotations'] += 1
                    continue
                relevance = 'SUPPORTED' if selected and (
                    selected.get('MANE_SELECT') or selected.get('MANE_PLUS_CLINICAL')
                ) else 'UNKNOWN'
                selected_priority = _consequence_priority(selected.get('Consequence', '') if selected else '')
                alternate_review = []
                for alt_entry in entries:
                    if alt_entry is selected:
                        continue
                    alt_priority = _consequence_priority(alt_entry.get('Consequence', ''))
                    if alt_priority <= selected_priority:
                        continue
                    alternate_review.append({
                        'transcript': alt_entry.get('Feature', ''),
                        'consequence': alt_entry.get('Consequence', ''),
                        'impact': alt_entry.get('IMPACT', ''),
                        'biotype': alt_entry.get('BIOTYPE', ''),
                        'mane_select': bool(alt_entry.get('MANE_SELECT')),
                        'mane_plus_clinical': bool(alt_entry.get('MANE_PLUS_CLINICAL')),
                        'canonical': str(alt_entry.get('CANONICAL', '')).upper() in ('YES', '1'),
                    })
                row = dict(gene_id=gene, symbol=selected.get('SYMBOL') or union('SYMBOL'), chrom=rec.contig, pos=rec.pos,
                           ref=rec.ref, alt=alt, gt=('|' if call.phased else '/').join('.' if x is None else str(x) for x in gt),
                           build=build, consequence=cons, transcript_relevance=relevance, af=af,
                           population_source='VEP_CACHE', DP=call.get('DP'), GQ=call.get('GQ'), AB=ab,
                           FILTER=';'.join(rec.filter) or '.', PS=call.get('PS'), PID=call.get('PID'), PGT=call.get('PGT'),
                           variant_class='SV' if sv else 'SMALL_VARIANT', end=rec.stop if sv else None,
                           transcripts=selected.get('Feature', ''), selected_transcript=selected.get('Feature', ''),
                           transcript_biotype=selected.get('BIOTYPE', ''), transcript_selection=selection,
                           alternate_transcript_effects=alternate_review,
                           all_transcript_consequences=all_cons,
                           cached_CLIN_SIG=union('CLIN_SIG'),
                           predictors={key: selected.get(key, '') for key in ('SIFT', 'PolyPhen', 'REVEL', 'CADD_PHRED')})
                counts['candidate_allele_gene_rows'] += 1
                yield from_row(row)


def apply_tie_ranges(rows, score_key, rank_key='rank'):
    """Annotate deterministic ranks with scientifically honest exact-score tie ranges."""
    i = 0
    while i < len(rows):
        score = float(rows[i][score_key])
        j = i + 1
        while j < len(rows) and float(rows[j][score_key]) == score:
            j += 1
        lo, hi = i + 1, j
        for k in range(i, j):
            rows[k]['tie_rank_min'] = lo
            rows[k]['tie_rank_max'] = hi
        i = j
    return rows


def freeze_g0(variants, out, config=ScoreConfig(), provenance=None, official=False, phase_evidence=()):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    provenance = provenance or {}
    grouped = defaultdict(list)
    for v in variants:
        grouped[v.gene].append(v)
    if not grouped:
        raise ValueError('No candidate variants; G0 was not frozen')
    evidence, hypotheses, ranks, review_queue, excluded_queue = [], [], [], [], []
    rank_fields = ('gene', 'symbol', 'hypothesis_id', 'inheritance_model',
                   'variant_1', 'variant_2', 'variant_1_score', 'variant_2_score', 'inheritance_score',
                   'case_score', 'case_status', 'phase_status', 'rationale',
                   'technical_rank_tier', 'pair_qc20_status', 'pair_qc20_reasons')
    for gene in sorted(grouped):
        ev, hs = evaluate_gene(grouped[gene], config, phase_evidence)
        evidence.extend(ev)
        hypotheses.extend(hs)

        primary = [h for h in hs
                   if h.get('technical_rank_tier', 'PRIMARY') == 'PRIMARY'
                   and h.get('case_status') not in ('INSUFFICIENT_EVIDENCE', 'CONFLICTING_EVIDENCE')
                   and float(h.get('case_score', 0)) > 0]
        review = [h for h in hs
                  if h.get('technical_rank_tier') == 'TECHNICAL_REVIEW'
                  and h.get('case_status') not in ('INSUFFICIENT_EVIDENCE', 'CONFLICTING_EVIDENCE')
                  and float(h.get('case_score', 0)) > 0]
        excluded = [h for h in hs if h.get('technical_rank_tier') == 'EXCLUDED']

        if primary:
            best = sorted(primary, key=lambda h: (-float(h['case_score']), h['hypothesis_id']))[0]
            ranks.append({k: best.get(k, '') for k in rank_fields})
        if review:
            best_review = sorted(review, key=lambda h: (-float(h['case_score']), h['hypothesis_id']))[0]
            review_queue.append({k: best_review.get(k, '') for k in rank_fields})
        for h in excluded:
            excluded_queue.append({k: h.get(k, '') for k in rank_fields})

    ranks.sort(key=lambda r: (-float(r['case_score']), r['gene']))
    review_queue.sort(key=lambda r: (-float(r['case_score']), r['gene']))
    excluded_queue.sort(key=lambda r: (r['gene'], r['hypothesis_id']))
    standards_review_queue = [
        e for e in evidence
        if e.get('standards_route') != 'SEQUENCE_VARIANT_RESEARCH_TRIAGE'
        or e.get('lof_review_status') != 'NOT_APPLICABLE'
        or e.get('transcript_review_status') != 'NOT_APPLICABLE'
        or e.get('computational_status') != 'NOT_ASSIGNED'
        or e.get('acmg_status') != 'NOT_CLASSIFIED'
    ]
    standards_review_queue.sort(key=lambda e: (e.get('gene',''), e.get('variant_key','')))
    for i, r in enumerate(ranks, 1):
        r['rank'] = i
        r['ranking_status'] = 'G0_COMPLETE_RESEARCH_HEURISTIC' if official else 'PROVISIONAL_FROZEN_INPUT_ONLY'
    apply_tie_ranges(ranks, 'case_score')
    rank_name = '06_genome_only_rank.csv' if official else 'provisional_genome_only_rank.csv'
    write_csv(out / rank_name, ranks)
    write_csv(out / 'technical_review_gene_queue.csv', review_queue)
    write_csv(out / 'technical_excluded_hypotheses.csv', excluded_queue)
    write_csv(out / 'standards_review_variant_queue.csv', standards_review_queue)
    write_csv(out / 'variant_evidence.csv', evidence)
    write_csv(out / 'inheritance_hypotheses.csv', hypotheses)
    write_json(out / 'score_config.json', asdict(config))
    write_json(out / 'g0_provenance.json', dict(provenance,
               created_utc=datetime.now(timezone.utc).isoformat(), phenotype_used=False,
               gene_disease_validity_used=False, official=official,
               phase_evidence=list(phase_evidence), acmg_classification_performed=False,
               variant_count=len(evidence), hypothesis_count=len(hypotheses),
               gene_count=len(ranks), primary_ranked_gene_count=len(ranks),
               technical_review_gene_count=len(review_queue),
               technical_excluded_hypothesis_count=len(excluded_queue),
               standards_review_variant_count=len(standards_review_queue),
               source_code_sha256={p.name: digest(p) for p in (Path(__file__), Path(__file__).with_name('inheritance.py'))}))
    write_json(out / 'audit_report.json', {
        'ranking_status': 'COMPLETE_G0_RESEARCH_HEURISTIC' if official else 'PROVISIONAL_FROZEN_INPUT_ONLY',
        'phenotype_used_in_g0': False, 'gene_specific_scoring': False,
        'acmg_classification_performed': False,
        'unknown_population_alleles': sum(e['population_status'] == 'UNKNOWN' for e in evidence),
        'unresolved_candidate_pairs': sum(h['variant_2'] != '' and h['phase_status'] in ('UNPHASED', 'PHASED_UNCERTAIN') for h in hypotheses),
        'phase_unresolved_pm3_assigned': any(h['phase_status'] in ('UNPHASED', 'PHASED_UNCERTAIN') and h['PM3'] != 'NOT_ASSIGNED' for h in hypotheses),
        'standards_review_variant_count': len(standards_review_queue),
        'clinical_criteria_assigned_in_g0': False,
        'technical_rank_policy': {
            'PRIMARY': ['QC_STRONG', 'QC_PASS'],
            'TECHNICAL_REVIEW': ['QC_REVIEW', 'QC_UNKNOWN'],
            'EXCLUDED': ['QC_FAIL'],
            'review_scores_preserved': True,
            'qc_used_as_pathogenicity_penalty': False,
        },
        'limitations': ['Molecular-impact tiers are research triage, not pathogenicity probabilities',
                        'G0 gene-disease mechanisms remain unknown',
                        'Clinical ACMG/AMP criteria are intentionally not assigned in phenotype-blind G0',
                        'Cached clinical labels and raw predictions are not ACMG classifications',
                        'CNV/SV, noncoding, splice-review and alternate-transcript evidence require dedicated review',
                        'Clinical assay validation and read-level confirmation are outside this VCF-only prioritization run'],
    })
    files = {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()}
    write_json(out / 'G0_FROZEN.json', {'version': config.version, 'rank_file': rank_name,
                                     'official': official, 'files': files})
    # Integrity is enforced by hashes and refusal to reuse output directories;
    # chmod alone is not an integrity boundary on Windows-mounted filesystems.
    return {'out': str(out), 'genes': len(ranks), 'hypotheses': len(hypotheses), 'official': official}


def verify_g0(folder):
    folder = Path(folder)
    manifest = json.loads((folder / 'G0_FROZEN.json').read_text())
    required = {manifest['rank_file'], 'inheritance_hypotheses.csv', 'variant_evidence.csv',
                'score_config.json', 'g0_provenance.json', 'audit_report.json'}
    if not required <= set(manifest.get('files', {})):
        raise ValueError('Incomplete G0 freeze manifest')
    for name, expected in manifest['files'].items():
        if Path(name).name != name or digest(folder / name) != expected:
            raise ValueError('G0 checksum mismatch: ' + name)
    return manifest


def read_csv(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def main(argv=None):
    p = argparse.ArgumentParser(description='Inheritance-aware G0; never reads an active annotation VCF')
    inputs = p.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--frozen-csv', type=Path)
    inputs.add_argument('--frozen-vcf', type=Path)
    inputs.add_argument('--completed-vcf', type=Path)
    p.add_argument('--sha256', help='Required checksum of frozen input')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--config', type=Path)
    p.add_argument('--phase-evidence', type=Path, help='Validated pair-specific phase observations as JSON list')
    p.add_argument('--sample')
    p.add_argument('--build', default='GRCh38')
    p.add_argument('--threads', type=int, default=4, choices=range(1, 17))
    args = p.parse_args(argv)
    c = ScoreConfig(**json.loads(args.config.read_text())) if args.config else ScoreConfig()
    phases = json.loads(args.phase_evidence.read_text()) if args.phase_evidence else []
    if args.completed_vcf:
        with completed_annotation(args.completed_vcf) as (checksum, state):
            counts = defaultdict(int)
            provenance = {'input': str(args.completed_vcf), 'input_sha256': checksum, 'counts': counts,
                          'sample': args.sample or state.get('identity', {}).get('sample'),
                          'scope': 'All qualifying annotated genome-wide candidates; unresolved structural/reference-incompatible calls remain separately unassessed'}
            result = freeze_g0(vcf_variants(args.completed_vcf, provenance['sample'], args.build, counts, threads=args.threads),
                               args.out, c, provenance, official=True, phase_evidence=phases)
    else:
        source = assert_snapshot_path(args.frozen_csv or args.frozen_vcf)
        if not args.sha256:
            p.error('Frozen inputs require --sha256; new snapshots of a live file are forbidden')
        actual = digest(source)
        if actual != args.sha256:
            raise ValueError('Frozen input checksum mismatch')
        variants = csv_variants(source) if args.frozen_csv else vcf_variants(source, args.sample, args.build, threads=args.threads)
        # Exhaust before freezing, then recheck that the supposedly frozen input stayed fixed.
        variants = list(variants)
        if digest(source) != actual:
            raise ValueError('Frozen input changed while reading')
        result = freeze_g0(variants, args.out, c, {'input': str(source), 'input_sha256': actual,
                           'scope': 'Frozen subset only; cannot infer missing alleles or genome-wide ranks'},
                           phase_evidence=phases)
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    main()

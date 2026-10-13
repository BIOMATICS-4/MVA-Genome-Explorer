"""Gene-agnostic, score-faithful Model 6 display projection.

A hypothesis's score, number of alleles, inheritance label and phase must all
come from the SAME ranker row. Other candidate alleles are separate evidence,
not a replacement for the ranked variant(s). No patient/candidate-specific rules.
This does NOT retrain or recalibrate the underlying Exomiser-compatible score.
"""
from __future__ import annotations

import csv
from pathlib import Path

PASS_QC = frozenset({'QC_PASS', 'QC_STRONG'})


def split_keys(value):
    return [x for x in str(value or '').split(';') if x]


def read_tsv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def qc_lookup(rows):
    """Join exact-key QC evidence without silently overriding a conflict.

    The compatibility and regulatory QC sources may include the same allele.
    Identical duplicate labels are allowed; contradictory labels fail closed.
    The quality status is technical evidence, never a pathogenicity assertion.
    """
    result = {}
    for r in rows:
        key, status = r.get('variant_key'), r.get('qc_status')
        if not key or not status:
            continue
        if key in result and result[key] != status:
            raise ValueError(f'Conflicting QC statuses for the same variant: {key}')
        result[key] = status
    return result


def _primary_status(moi, phase, count):
    if moi == 'AR' and count == 1 and phase == 'HOMOZYGOUS_ALT':
        return 'RECESSIVE_HOMOZYGOUS_ALT'
    if moi == 'AR' and count == 2 and phase == 'TRANS':
        return 'RECESSIVE_CONFIRMED_TRANS'
    if moi == 'AD' and count == 1:
        return 'SINGLE_ALLELE_AD_STYLE_SCORE_NOT_DIAGNOSIS'
    return 'HYPOTHESIS_REQUIRES_REVIEW'


def consistent_row(strict, raw_rows, qc):
    """Describe exactly the strict rank; keep distinct AR hypotheses separate.

    Raw pairs failing technical QC remain review-only, not promoted as an
    accepted alternative. Unresolved phase cannot be treated as confirmed trans.
    """
    gene = strict['gene']
    moi = strict.get('moi') or strict.get('best_moi') or ''
    keys = split_keys(strict.get('variant_keys'))
    if not keys:
        raise ValueError(f'Strict row for {gene} has no scoring variant keys')
    primary_phase = strict.get('strict_phase_status') or strict.get('phase_status') or 'NOT_APPLICABLE'
    primary_labels = split_keys(strict.get('strict_qc_variant_statuses'))
    if primary_labels and len(primary_labels) != len(keys):
        raise ValueError(f'Strict row for {gene} has mismatched QC and variant counts')
    primary_qc = primary_labels or [qc.get(k, 'QC_UNKNOWN') for k in keys]
    if any(k in qc and qc[k] != s for k, s in zip(keys, primary_qc)):
        raise ValueError(f'Strict row for {gene} conflicts with supplementary QC records')
    if any(x not in PASS_QC for x in primary_qc):
        raise ValueError(f'Strict row for {gene} contains nonpassing QC status: {primary_qc}')

    result = {
        'gene_rank': strict['gene_rank'],
        'gene': gene,
        'ranked_hypothesis_id': strict.get('id') or f'{gene}_{moi}',
        'ranked_moi': moi,
        'ranked_score': strict['combined_score'],
        'ranked_variant_score': strict.get('variant_score', ''),
        'ranked_phenotype_score': strict.get('phenotype_score', ''),
        'ranked_phase': primary_phase,
        'ranked_variant_count': str(len(keys)),
        'ranked_variant_keys': ';'.join(keys),
        'ranked_effects': strict.get('effects', ''),
        'ranked_qc': ';'.join(primary_qc),
        'ranked_status': _primary_status(moi, primary_phase, len(keys)),
        'additional_hypothesis_moi': '',
        'additional_hypothesis_score': '',
        'additional_hypothesis_phase': '',
        'additional_hypothesis_variant_count': '',
        'additional_hypothesis_variant_keys': '',
        'additional_hypothesis_effects': '',
        'additional_hypothesis_support_sources': '',
        'additional_hypothesis_classification_flag': '',
        'additional_hypothesis_frequency_status': '',
        'additional_hypothesis_qc': '',
        'additional_hypothesis_status': '',
        'review_only_variant_keys': '',
        'review_only_qc': '',
    }
    alts = []
    review = []
    for candidate in raw_rows:
        if candidate.get('gene') != gene or candidate.get('moi') != 'AR':
            continue
        ckeys = split_keys(candidate.get('variant_keys'))
        if len(ckeys) != 2 or len(set(ckeys)) != 2:
            continue
        if set(ckeys) == set(keys):
            continue
        # Only describe pairs that contain the already ranked allele. Do not
        # swap a ranked variant for an unrelated pair for cosmetic reasons.
        if not set(keys).issubset(set(ckeys)):
            continue
        phase = str(candidate.get('phase_status') or '').upper()
        if phase not in ('TRANS', 'UNRESOLVED'):
            continue
        statuses = [qc.get(k, 'QC_UNKNOWN') for k in ckeys]
        if all(s in PASS_QC for s in statuses):
            alts.append((float(candidate['combined_score']), candidate, ckeys, phase, statuses))
        else:
            # A second allele with weak QC is not promoted to a scored pair.
            review.append((candidate, ckeys, statuses))

    if alts:
        _, alt, ckeys, phase, statuses = max(alts, key=lambda x:(x[0], x[1].get('id', '')))
        result.update({
            'additional_hypothesis_moi': 'AR',
            'additional_hypothesis_score': alt['combined_score'],
            'additional_hypothesis_phase': phase,
            'additional_hypothesis_variant_count': str(len(ckeys)),
            'additional_hypothesis_variant_keys': ';'.join(ckeys),
            'additional_hypothesis_effects': alt.get('effects', ''),
            'additional_hypothesis_support_sources': alt.get('pathogenicity_sources', ''),
            'additional_hypothesis_frequency_status': (
                'AT_LEAST_ONE_AF_UNKNOWN' if any(
                    x.strip() in ('', '.') for x in
                    str(alt.get('max_freq_pct', '')).split(';')
                ) else 'AF_REPORTED_FOR_EACH_VARIANT'
            ),
            'additional_hypothesis_classification_flag': (
                'ONE_OR_MORE_ALLELES_NOT_INDEPENDENTLY_CLASSIFIED'
                if any(s != 'CLINVAR_PATHOGENIC' for s in
                       str(alt.get('pathogenicity_sources') or '').split(';'))
                else 'CLINVAR_LABELS_PRESENT_NOT_INDEPENDENTLY_REVIEWED'
            ),
            'additional_hypothesis_qc': ';'.join(statuses),
            'additional_hypothesis_status': (
                'TRANS_PHASE_CONFIRMED_CAUSALITY_UNASSESSED' if phase == 'TRANS'
                else 'CANDIDATE_PAIR_PHASE_UNRESOLVED'
            ),
        })
    if review:
        # Review-only calls are surfaced, never mixed with the strict score.
        seen = set(keys)
        extra = []
        for _, ckeys, statuses in review:
            for k, s in zip(ckeys, statuses):
                if k not in seen:
                    seen.add(k)
                    extra.append((k, s))
        result['review_only_variant_keys'] = ';'.join(k for k, _ in extra)
        result['review_only_qc'] = ';'.join(s for _, s in extra)
    return result


def build_consistent_rows(strict_rows, raw_rows, qc_rows):
    qc = qc_lookup(qc_rows)
    by_gene = {}
    for row in raw_rows:
        by_gene.setdefault(row.get('gene'), []).append(row)
    return [consistent_row(row, by_gene.get(row['gene'], []), qc) for row in strict_rows]


def write_tsv(path, rows):
    if not rows:
        raise ValueError('No strict gene ranking to write')
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

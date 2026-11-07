"""Generic standards-backed evidence adapters.

No challenge genes, coordinates, diseases, or case-specific thresholds belong here.
External assertions are converted to explicit evidence states; they are never treated
as sequence-variant pathogenicity classifications.
"""
from __future__ import annotations
from collections import defaultdict
import csv
from pathlib import Path

HIGH_VALIDITY = {'Definitive', 'Strong', 'Moderate'}
CONFLICT_VALIDITY = {'Disputed', 'Disputed Evidence', 'Refuted', 'Refuted Evidence'}
CONTEXT_ONLY_VALIDITY = {
    'Supportive', 'Limited', 'Animal Model Only', 'No Known Disease Relationship',
    'No known disease relationship',
}

def normalize_disease_id(value: str) -> str:
    value = (value or '').strip()
    if value.lower().startswith('orphanet:'):
        return 'ORPHA:' + value.split(':', 1)[1]
    return value

def moi_models(title: str):
    t = (title or '').strip().lower()
    if 'autosomal recessive' in t:
        return ('AR_HOMOZYGOUS', 'AR_POSSIBLE_COMPOUND_HET')
    if 'autosomal dominant' in t:
        return ('AD_SINGLE_ALLELE', 'DE_NOVO_COMPATIBLE')
    if 'x-linked' in t or 'x linked' in t:
        return ('X_LINKED',)
    return ()

def _disease_ids(row):
    values = {
        row.get('disease_curie', ''),
        row.get('disease_original_curie', ''),
        row.get('submitted_as_disease_id', ''),
    }
    return {normalize_disease_id(v) for v in values if v}

def gencc_validity_rows(path, genes, allowed_diseases=None, version='GENCC_CURRENT_DOWNLOAD'):
    # Aggregate GenCC gene-disease-MOI assertions conservatively.
    # Equivalent identifiers are linked through GenCC disease_curie clusters.
    symbol_to_gene = defaultdict(list)
    for gene, symbol in genes.items():
        if symbol:
            symbol_to_gene[symbol].append(gene)

    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        rows_all = list(csv.DictReader(f))

    aliases_by_canonical = defaultdict(set)
    for row in rows_all:
        canonical = normalize_disease_id(row.get('disease_curie', ''))
        ids = _disease_ids(row)
        if canonical:
            aliases_by_canonical[canonical].update(ids | {canonical})

    grouped = defaultdict(list)
    allowed = set(allowed_diseases) if allowed_diseases is not None else None
    for row in rows_all:
        symbol = row.get('gene_symbol', '')
        if symbol not in symbol_to_gene:
            continue
        models = moi_models(row.get('moi_title', ''))
        if not models:
            continue
        canonical = normalize_disease_id(row.get('disease_curie', ''))
        disease_ids = set(_disease_ids(row))
        if canonical:
            disease_ids |= aliases_by_canonical.get(canonical, set())
        if allowed is not None:
            disease_ids &= allowed
        for gene in symbol_to_gene[symbol]:
            for disease in disease_ids:
                for model in models:
                    grouped[(gene, disease, model)].append(row)

    out = []
    for (gene, disease, model), rows in sorted(grouped.items()):
        classes = {r.get('classification_title', '').strip() for r in rows}
        high_rows = [r for r in rows if r.get('classification_title','').strip() in HIGH_VALIDITY]
        conflict_rows = [r for r in rows if r.get('classification_title','').strip() in CONFLICT_VALIDITY]
        high_submitters = {r.get('submitter_title','').strip() for r in high_rows if r.get('submitter_title','').strip()}
        conflict_submitters = {r.get('submitter_title','').strip() for r in conflict_rows if r.get('submitter_title','').strip()}
        clingen_high = any(r.get('submitter_title','').strip() == 'ClinGen' for r in high_rows)
        clingen_conflict = any(r.get('submitter_title','').strip() == 'ClinGen' for r in conflict_rows)
        independently_conflicting = clingen_conflict or len(conflict_submitters) >= 2
        # Fail closed: GenCC is an aggregation layer, not itself a replacement
        # for ClinGen's gene-disease validity framework. A high-validity ClinGen
        # assertion with matching MOI is sufficient to unlock G1 case-fit.
        # Multiple non-ClinGen high assertions remain useful context but do not
        # silently become ClinGen-equivalent validity.
        support = ('CONFLICTING' if independently_conflicting else
                   'SUPPORTED' if clingen_high else 'UNKNOWN')
        support_basis = ('CLINGEN_CONFLICT_OR_MULTISOURCE_CONFLICT' if independently_conflicting else
                         'CLINGEN_HIGH_VALIDITY_MATCHING_MOI' if clingen_high else
                         'MULTISOURCE_HIGH_CONTEXT_ONLY' if len(high_submitters) >= 2 else
                         'CONTEXT_ONLY')
        out.append({
            'gene': gene,
            'disease_id': disease,
            'inheritance_model': model,
            'inheritance_support': support,
            'source': 'GenCC',
            'version': version,
            'validity_classifications': ';'.join(sorted(c for c in classes if c)),
            'validity_submitters': ';'.join(sorted({r.get('submitter_title','') for r in rows if r.get('submitter_title')})),
            'high_validity_submitters': ';'.join(sorted(high_submitters)),
            'high_validity_submitter_count': len(high_submitters),
            'clingen_high_validity': clingen_high,
            'clingen_conflicting_validity': clingen_conflict,
            'support_basis': support_basis,
            'non_clingen_high_context_only': bool(len(high_submitters) >= 2 and not clingen_high),
            'conflict_submitters': ';'.join(sorted(conflict_submitters)),
            'assertion_count': len(rows),
        })
    return out


def hypothesis_population_compatibility(hypothesis, variant_evidence,
                                        review_af=0.01,
                                        general_high_af_review=0.05):
    """Population-frequency triage for rare-disease prioritization.

    This is deliberately not ACMG/AMP criterion assignment. It prevents an
    obviously common allele from being promoted by phenotype fit while retaining
    the hypothesis for review. Missing/unreliable population frequency remains
    neutral, never evidence of rarity.

    A disease-specific maximum credible AF may be supplied by the caller by
    replacing review_af. The general >5% screen is review-only because ClinGen
    BA1 has explicit exceptions and requires criterion-level review.
    """
    try:
        review_af = float(review_af)
        general_high_af_review = float(general_high_af_review)
    except (TypeError, ValueError):
        raise ValueError('Population review thresholds must be numeric')
    if not (0 <= review_af <= 1 and 0 <= general_high_af_review <= 1):
        raise ValueError('Population review thresholds must be in [0,1]')
    if review_af > general_high_af_review:
        raise ValueError('Rare-disease review threshold cannot exceed general high-AF screen')

    observations = []
    unknown = 0
    for field in ('variant_1', 'variant_2'):
        key = str(hypothesis.get(field, '') or '')
        if not key:
            continue
        e = variant_evidence.get(key, {})
        status = str(e.get('population_status', 'UNKNOWN') or 'UNKNOWN')
        raw = e.get('observed_af')
        af = None
        try:
            if raw not in (None, '', '.', 'nan', 'NaN'):
                af = float(raw)
                if not 0 <= af <= 1:
                    af = None
        except (TypeError, ValueError):
            af = None
        if status == 'OBSERVED' and af is not None:
            observations.append((key, af))
        else:
            unknown += 1

    max_af = max((af for _, af in observations), default=None)
    if max_af is not None and max_af > general_high_af_review:
        state = 'VERY_HIGH_AF_REVIEW'
        eligible = False
        reason = 'Observed population AF exceeds general high-frequency review screen'
    elif max_af is not None and max_af > review_af:
        state = 'RARE_DISEASE_HIGH_AF_REVIEW'
        eligible = False
        reason = 'Observed population AF exceeds configured rare-disease review threshold'
    else:
        state = 'COMPATIBLE_OR_UNKNOWN'
        eligible = True
        reason = ('No observed AF exceeds review threshold; unknown AF remains neutral'
                  if unknown else 'Observed AF does not exceed review threshold')

    return {
        'state': state,
        'rank_eligible': eligible,
        'max_observed_af': max_af,
        'observed_allele_count': len(observations),
        'unknown_allele_count': unknown,
        'review_af': review_af,
        'general_high_af_review': general_high_af_review,
        'reason': reason,
        'formal_population_criterion': 'NOT_ASSIGNED',
        'semantics': 'RESEARCH_REVIEW_GATE_NOT_BA1_BS1_PM2_CLASSIFICATION',
    }

def annotation_disease_ids(path):
    ids = set()
    with Path(path).open(newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('disease_id'):
                ids.add(normalize_disease_id(row['disease_id']))
    return ids

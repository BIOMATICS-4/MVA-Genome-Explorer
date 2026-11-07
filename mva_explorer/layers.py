"""Explicit, additive G1/G2 overlays on a checksum-verified frozen G0.

Evidence is supplied as versioned data. These layers never edit or reclassify G0.
"""
from collections import defaultdict
from pathlib import Path
import math

from .ranking import digest, read_csv, verify_g0, write_csv, write_json, apply_tie_ranges
from .standards import hypothesis_population_compatibility


def overlay(g0, out, evidence, layer='G1', weight=10.0, previous=None, provenance=None):
    if layer not in ('G1', 'G2') or not math.isfinite(weight) or weight < 0:
        raise ValueError('Invalid layer or weight')
    g0, out = Path(g0), Path(out)
    manifest = verify_g0(g0)
    freeze_hash = digest(g0 / 'G0_FROZEN.json')
    rows = read_csv(g0 / 'inheritance_hypotheses.csv')
    genes = {r['gene']: r for r in read_csv(g0 / manifest['rank_file'])}
    variant_evidence = {r['variant_key']: r for r in read_csv(g0 / 'variant_evidence.csv')}
    population_review = []
    genotype_review = []

    def truthy(value):
        return value is True or str(value).strip().lower() in ('1', 'true', 'yes')

    previous_rows = {}
    if layer == 'G2':
        if previous is None:
            raise ValueError('G2 requires a frozen G1 overlay')
        import json
        previous = Path(previous)
        pm = json.loads((previous / 'layer_manifest.json').read_text())
        if pm.get('layer') != 'G1' or pm.get('g0_manifest_sha256') != freeze_hash:
            raise ValueError('G1 was not derived from this G0')
        for name, sha in pm['files'].items():
            if Path(name).name != name or digest(previous / name) != sha:
                raise ValueError('G1 checksum mismatch')
        previous_rows = {r['hypothesis_id']: r for r in read_csv(previous / 'hypothesis_ranks.csv')}
        if set(previous_rows) != {r['hypothesis_id'] for r in rows}:
            raise ValueError('G1 hypothesis identities differ from G0')
    grouped = defaultdict(list)
    used_evidence = []
    for e in evidence:
        e = dict(e)
        if not all(e.get(k) for k in ('gene', 'source', 'version', 'dependency_group')):
            raise ValueError('Layer evidence requires gene, source, version and dependency_group')
        score = float(e['score'])
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError('Layer evidence score must be in [0,1]')
        if layer == 'G1' and not e.get('disease_id'):
            raise ValueError('G1 evidence must identify a disease profile')
        if e.get('inheritance_support', 'UNKNOWN') not in ('SUPPORTED', 'UNKNOWN', 'CONFLICTING'):
            raise ValueError('Invalid inheritance support state')
        grouped[e['gene']].append(e)
        used_evidence.append(e)
    results = []
    for row in rows:
        row = dict(row)
        if row.get('technical_rank_tier', 'PRIMARY') != 'PRIMARY':
            continue
        if row['gene'] not in genes:
            continue
        if layer == 'G2' and not truthy(previous_rows[row['hypothesis_id']].get('G1_primary_rank_eligible', True)):
            # Review-only G1 hypotheses stay in the audit ledger and cannot be
            # reintroduced by a later contextual layer.
            continue
        applicable = [e for e in grouped[row['gene']] if not e.get('inheritance_model')
                      or e['inheritance_model'] == row['inheritance_model']]
        # Choose one disease explanation in G1, never combine unrelated diseases.
        # Correlated systems evidence in G2 is conservatively max-pooled too.
        support_rank = {'CONFLICTING': 0, 'UNKNOWN': 1, 'SUPPORTED': 2}
        best = max(applicable, key=lambda e: (
            support_rank.get(e.get('inheritance_support', 'UNKNOWN'), 1),
            float(e['score']), e.get('disease_id', ''), e['source']
        ), default=None)
        support = best.get('inheritance_support', 'UNKNOWN') if best else 'UNKNOWN'
        similarity = float(best['score']) if best else 0
        before_status = previous_rows[row['hypothesis_id']]['G1_status'] if layer == 'G2' else row['case_status']
        supported = before_status not in ('INSUFFICIENT_EVIDENCE', 'CONFLICTING_EVIDENCE')
        before = float(previous_rows[row['hypothesis_id']]['G1_case_score']) if layer == 'G2' else float(row['case_score'])
        population = hypothesis_population_compatibility(row, variant_evidence) if layer == 'G1' else None
        genotype_eligible = bool(row.get('variant_1')) and supported if layer == 'G1' else True
        primary_eligible = bool(population['rank_eligible']) and genotype_eligible if layer == 'G1' else True
        if layer == 'G1':
            # Standards-first: phenotype case fit is a separate dimension.
            # Never convert it into invented equivalents of molecular-impact points.
            case_fit = similarity if supported and support == 'SUPPORTED' else 0.0
            support_tier = 2 if support == 'SUPPORTED' else 0 if support == 'CONFLICTING' else 1
            row.update({
                'G1_case_score': round(case_fit, 6),
                'G1_case_fit_score': round(case_fit, 6),
                'G1_molecular_score': float(row['case_score']),
                'G1_evidence_score': similarity,
                'G1_contribution': 0.0,
                'G1_inheritance_support': support,
                'G1_disease_id': best.get('disease_id', '') if best else '',
                'G1_source': best['source'] if best else '',
                'G1_status': 'CONFLICTING_EVIDENCE' if support == 'CONFLICTING' else before_status,
                'G1_support_tier': support_tier,
                'G1_population_state': population['state'],
                'G1_population_rank_eligible': population['rank_eligible'],
                'G1_genotype_rank_eligible': genotype_eligible,
                'G1_primary_rank_eligible': primary_eligible,
                'G1_population_max_observed_af': population['max_observed_af'],
                'G1_population_review_reason': population['reason'],
                'G1_population_semantics': population['semantics'],
                'G1_formal_population_criterion': population['formal_population_criterion'],
                'G0_rank': int(genes[row['gene']]['rank']),
                'G0_case_score': float(row['case_score']),
            })
            if not population['rank_eligible']:
                population_review.append(dict(row))
            if not genotype_eligible:
                genotype_review.append(dict(row))
        else:
            # G2 is contextual annotation only. It cannot alter primary G1 order.
            context_score = similarity if supported and support != 'CONFLICTING' else 0.0
            prev = previous_rows[row['hypothesis_id']]
            row.update({
                'G1_case_score': float(prev['G1_case_score']),
                'G1_case_fit_score': float(prev.get('G1_case_fit_score', prev['G1_case_score'])),
                'G1_molecular_score': float(prev.get('G1_molecular_score', row['case_score'])),
                'G1_hypothesis_rank': int(prev['G1_hypothesis_rank']),
                'G2_case_score': float(prev['G1_case_score']),
                'G2_evidence_score': similarity,
                'G2_context_score': round(context_score, 6),
                'G2_contribution': 0.0,
                'G2_inheritance_support': support,
                'G2_disease_id': best.get('disease_id', '') if best else '',
                'G2_source': best['source'] if best else '',
                'G2_status': 'CONFLICTING_EVIDENCE' if support == 'CONFLICTING' else before_status,
                'G0_rank': int(genes[row['gene']]['rank']),
                'G0_case_score': float(row['case_score']),
            })
        results.append(row)
    if layer == 'G1':
        # Keep rejected hypotheses in the auditable hypothesis ledger, but put
        # every population-compatible hypothesis ahead of review-only entries.
        results.sort(key=lambda r: (not truthy(r.get('G1_primary_rank_eligible', True)),
                                   -int(r['G1_support_tier']),
                                   -float(r['G1_case_fit_score']),
                                   -float(r['G1_molecular_score']),
                                   r['hypothesis_id']))
    else:
        results.sort(key=lambda r: (int(r['G1_hypothesis_rank']), r['hypothesis_id']))
    best_genes = {}
    for i, row in enumerate(results, 1):
        row[f'{layer}_hypothesis_rank'] = i
        if layer == 'G1' and not truthy(row.get('G1_primary_rank_eligible', True)):
            continue
        best_genes.setdefault(row['gene'], row)
    previous_genes = {r['gene']: r for r in read_csv(previous / 'gene_rank_comparison.csv')} if previous else {}
    comparisons = []
    for rank, (gene, row) in enumerate(best_genes.items(), 1):
        old = genes[gene]
        comparison = {'gene': gene, 'symbol': old['symbol'], 'G0_rank': int(old['rank']),
                      'G0_case_score': float(old['case_score']), f'{layer}_rank': rank,
                      f'{layer}_case_score': row[f'{layer}_case_score'],
                      f'{layer}_best_hypothesis': row['hypothesis_id'],
                      'rank_change_from_G0': int(old['rank']) - rank,
                      f'{layer}_status': row[f'{layer}_status'],
                      'ranking_status': 'RESEARCH_HEURISTIC' if manifest['official'] else 'PROVISIONAL_FROZEN_SUBSET'}
        if layer == 'G1':
            comparison['G1_case_fit_score'] = row['G1_case_fit_score']
            comparison['G1_molecular_score'] = row['G1_molecular_score']
            comparison['G1_inheritance_support'] = row['G1_inheritance_support']
            comparison['G1_disease_id'] = row['G1_disease_id']
            comparison['primary_score'] = row['G1_case_fit_score']
            comparison['primary_score_layer'] = 'G1_CASE_FIT'
        if previous_genes:
            comparison['G1_rank'] = int(previous_genes[gene]['G1_rank'])
            # Preserve the pre-mechanism gene score even if G2 chooses a
            # different hypothesis from the one that won in G1.
            comparison['G1_case_score'] = float(previous_genes[gene]['G1_case_score'])
            comparison['primary_score'] = comparison['G1_case_score']
            comparison['primary_score_layer'] = 'G1'
            comparison['G2_contribution'] = row['G2_contribution']
            comparison['rank_change_from_G1'] = comparison['G1_rank'] - rank
        comparisons.append(comparison)
    if layer == 'G1':
        i = 0
        while i < len(comparisons):
            a = comparisons[i]
            key = (a.get('G1_inheritance_support','UNKNOWN'),
                   float(a.get('G1_case_fit_score',0)),
                   float(a.get('G1_molecular_score',0)))
            j = i + 1
            while j < len(comparisons):
                b = comparisons[j]
                bkey = (b.get('G1_inheritance_support','UNKNOWN'),
                        float(b.get('G1_case_fit_score',0)),
                        float(b.get('G1_molecular_score',0)))
                if bkey != key:
                    break
                j += 1
            for k in range(i,j):
                comparisons[k]['tie_rank_min'] = i + 1
                comparisons[k]['tie_rank_max'] = j
            i = j
    else:
        # G2 preserves G1 primary ordering; its context score is not a rank score.
        for r in comparisons:
            r['tie_rank_min'] = r[f'{layer}_rank']
            r['tie_rank_max'] = r[f'{layer}_rank']
    # Fail before creating outputs if the baseline changed during calculation.
    verify_g0(g0)
    if digest(g0 / 'G0_FROZEN.json') != freeze_hash:
        raise ValueError('G0 manifest changed during overlay')
    out.mkdir(parents=True, exist_ok=False)
    write_csv(out / 'hypothesis_ranks.csv', results)
    write_csv(out / 'gene_rank_comparison.csv', comparisons)
    if layer == 'G1':
        write_csv(out / 'population_review_hypotheses.csv', population_review)
        write_csv(out / 'genotype_review_hypotheses.csv', genotype_review)
    write_json(out / 'layer_evidence.json', used_evidence)
    if layer == 'G2':
        write_csv(out / 'final_top_candidates.csv', comparisons[:10])
        write_json(out / 'final_audit_report.json', {
            'g0_manifest_sha256': freeze_hash, 'g0_verified': True,
            'official_g0': manifest['official'], 'g0_unchanged': True,
            'case_status_is_not_acmg': True, 'acmg_classification_performed': False,
            'missing_systems_evidence_genes': sorted(set(genes) - set(grouped)),
            'prioritization_only': True, 'g2_weight': weight,
        })
    write_json(out / 'layer_manifest.json', {'layer': layer, 'weight': weight,
               'score_semantics': ('LEXICOGRAPHIC_SUPPORTED_VALIDITY_THEN_CASE_FIT_THEN_MOLECULAR_IMPACT'
                                   if layer == 'G1' else
                                   'CONTEXT_ONLY_PRESERVES_G1_PRIMARY_ORDER'),
               'additive_cross_evidence_scoring': False,
               'population_gate': ('HYPOTHESIS_LEVEL_RESEARCH_REVIEW' if layer == 'G1' else 'INHERITED_FROM_G1'),
               'population_review_hypotheses': (len(population_review) if layer == 'G1' else None),
               'genotype_review_hypotheses': (len(genotype_review) if layer == 'G1' else None),
               'pooling': 'maximum single disease/evidence profile; no correlated summation',
               'g0_manifest_sha256': freeze_hash, 'g0': str(g0.resolve()),
               'provenance': provenance or {},
               'previous_layer_manifest_sha256': digest(previous / 'layer_manifest.json') if previous else None,
               'files': {p.name: digest(p) for p in out.iterdir() if p.is_file()}})
    return comparisons

"""Read-only, outcome-agnostic influence audit of a completed ranking.

Scores are research heuristics. Counterfactuals measure sensitivity, not accuracy.
No candidate identities, coordinates, diagnoses, or optimized thresholds occur here.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import contextmanager
from dataclasses import replace
from datetime import datetime, timezone
import csv
import json
import math
from pathlib import Path

from . import inheritance as inh
from . import qc_integrated_v2 as qc
from .ranking import digest, read_csv, verify_g0, write_csv, write_json

NS = (10, 20, 50)
BAD = {'INSUFFICIENT_EVIDENCE', 'CONFLICTING_EVIDENCE'}


def rank_rows(rows):
    """One consistent, identity-based tie breaker; expose full equal-score ranges."""
    rows = sorted((dict(r) for r in rows), key=lambda r: (-r['score'], r['gene']))
    bounds = {}
    for i, r in enumerate(rows, 1):
        bounds.setdefault(r['score'], [i, i])[1] = i
    for i, r in enumerate(rows, 1):
        r.update(rank=i, tie_rank_min=bounds[r['score']][0], tie_rank_max=bounds[r['score']][1])
    return rows


def compare(before, after, name):
    b = {r['gene']: r for r in before}
    a = {r['gene']: r for r in after}
    rows = []
    for gene in sorted(set(b) | set(a)):
        x, y = b.get(gene, {}), a.get(gene, {})
        rows.append(dict(comparison=name, gene=gene, symbol=y.get('symbol', x.get('symbol', '')),
                         rank_before=x.get('rank'), rank_after=y.get('rank'),
                         score_before=x.get('score'), score_after=y.get('score'),
                         rank_change=x['rank']-y['rank'] if x and y else None,
                         score_change=round(y['score']-x['score'], 6) if x and y else None,
                         tie_min_before=x.get('tie_rank_min'), tie_max_before=x.get('tie_rank_max'),
                         tie_min_after=y.get('tie_rank_min'), tie_max_after=y.get('tie_rank_max'),
                         best_hypothesis_before=x.get('hypothesis_id', ''),
                         best_hypothesis_after=y.get('hypothesis_id', '')))
    metrics = dict(comparison=name, score_changed_genes=sum(r['score_change'] not in (None, 0) for r in rows),
                   rank_changed_genes=sum(r['rank_change'] not in (None, 0) for r in rows),
                   positive_genes_before=sum(r['score'] > 0 for r in before),
                   positive_genes_after=sum(r['score'] > 0 for r in after))
    shared = sorted(set(a) & set(b))
    x = [(b[g]['tie_rank_min']+b[g]['tie_rank_max'])/2 for g in shared]
    y = [(a[g]['tie_rank_min']+a[g]['tie_rank_max'])/2 for g in shared]
    if shared:
        xm, ym = sum(x)/len(x), sum(y)/len(y)
        den = math.sqrt(sum((v-xm)**2 for v in x)*sum((v-ym)**2 for v in y))
        metrics['spearman_average_tied_ranks'] = sum((v-xm)*(w-ym) for v,w in zip(x,y))/den if den else None
    names = lambda genes: [dict(gene=g, symbol=a.get(g, b.get(g))['symbol']) for g in sorted(genes)]
    for n in NS:
        bs, ass = {r['gene'] for r in before[:n]}, {r['gene'] for r in after[:n]}
        cutoff_b, cutoff_a = before[min(n, len(before))-1]['score'], after[min(n, len(after))-1]['score']
        metrics[f'top{n}'] = dict(overlap=len(bs & ass), entered=names(ass-bs), left=names(bs-ass),
                                  cutoff_before=cutoff_b, cutoff_after=cutoff_a,
                                  tie_inclusive_before=sum(r['score'] >= cutoff_b for r in before),
                                  tie_inclusive_after=sum(r['score'] >= cutoff_a for r in after),
                                  after_topN_in_before_cutoff_pool=sum(b[g]['score'] >= cutoff_b for g in ass if g in b),
                                  before_topN_in_after_cutoff_pool=sum(a[g]['score'] >= cutoff_a for g in bs if g in a))
    return rows, metrics


def hypothesis_genes(hypotheses, score_fn=None, allowed=None, production_first=False):
    grouped = defaultdict(list)
    for h in hypotheses:
        grouped[h['gene']].append(h)
    result = []
    for gene, hs in grouped.items():
        candidates = [h for h in hs if allowed is None or allowed(h)]
        scored = [(float(score_fn(h) if score_fn else h['case_score']), h) for h in candidates]
        if not scored:
            scored = [(0.0, hs[0])]
        score, h = scored[0] if production_first else min(scored, key=lambda x: (-x[0], x[1]['hypothesis_id']))
        result.append(dict(gene=gene, symbol=h['symbol'], score=round(score, 6),
                           hypothesis_id=h['hypothesis_id'], inheritance_model=h['inheritance_model'],
                           variant_1=h['variant_1'], variant_2=h['variant_2'], status=h['case_status']))
    return rank_rows(result)


@contextmanager
def replay_engine(stored, integrated=True, technical=True):
    """Run actual production inheritance/QC code with saved QC observations.

    The frozen table retains QC results, not all site measurements. This replays
    their application exactly; it does not claim to remeasure BAM/VCF quality.
    """
    original_ev, original_gene = inh.variant_evidence, inh.evaluate_gene
    original_qc, original_assess = inh.qc_support, qc.assess_allele
    lookup = {e['variant_key']: e for e in stored}
    try:
        if not technical:
            inh.qc_support = lambda v, c: ('PASS', [])
        if integrated:
            def cached(v, *_args, **_kwargs):
                e = lookup[v.key]
                return dict(status=e['qc20_status'], penalty=float(e['qc20_adjustment']),
                            reasons=e['qc20_reasons'].split(';'),
                            near500_count=int(e['qc20_near500_count']),
                            near500_lowab_het_count=int(e['qc20_near500_lowab_het_count']),
                            near100_count=int(e['qc20_near100_count']),
                            near100_indel_count=int(e['qc20_near100_indel_count']),
                            depth_ratio=inh.number(e['qc20_depth_ratio']))
            qc.assess_allele = cached
            qc.install_runtime_hooks(inh, None, None, {})
        yield inh.evaluate_gene
    finally:
        inh.variant_evidence, inh.evaluate_gene = original_ev, original_gene
        inh.qc_support, qc.assess_allele = original_qc, original_assess


def replay(stored, config, phases, mode='baseline'):
    grouped = defaultdict(list)
    for e in stored:
        v = inh.Variant(**json.loads(e['input_evidence']))
        if mode == 'population_masked':
            v = replace(v, af=None, population_status='UNKNOWN', complete_callable_resource=False)
        elif mode == 'transcript_masked':
            v = replace(v, transcript_relevance='UNKNOWN')
        elif mode == 'clinical_masked':
            v = replace(v, reviewed_acmg={}, computational_review={}, cached_clinical='', predictors={})
        grouped[v.gene].append(v)
    hypotheses, evidence = [], []
    with replay_engine(stored, integrated=mode not in ('no_integrated_qc', 'no_technical_qc'),
                       technical=mode != 'no_technical_qc') as evaluate:
        for gene in sorted(grouped):
            ev, hs = evaluate(grouped[gene], config, phases)
            evidence.extend(ev)
            hypotheses.extend(hs)
    return evidence, hypotheses


def variant_genes(evidence, score_fn, eligible=False):
    best = {}
    for e in evidence:
        value = float(score_fn(e)) if not eligible or e['variant_strength'] != 'UNSUPPORTED' else 0.0
        r = dict(gene=e['gene'], symbol=e['symbol'], score=round(value, 6), variant_1=e['variant_key'])
        if e['gene'] not in best or (-r['score'], r['variant_1']) < (-best[e['gene']]['score'], best[e['gene']]['variant_1']):
            best[e['gene']] = r
    return rank_rows(best.values())


def overlay_genes(hypotheses, evidence, weight=10, only=False, base_score=None,
                  mode='standards'):
    """Replay G1 evidence with explicit semantics.

    standards is fail-closed: validity -> phenotype case fit -> molecular impact.
    Cross-evidence point addition is forbidden.
    legacy_additive exists only to reproduce historical frozen audits.
    """
    by_gene = defaultdict(list)
    for e in evidence:
        by_gene[e['gene']].append(e)

    if mode == 'legacy_additive':
        def score(h):
            es = [e for e in by_gene[h['gene']]
                  if not e.get('inheritance_model') or e['inheritance_model'] == h['inheritance_model']]
            e = max(es, key=lambda e: (
                e.get('inheritance_support') != 'CONFLICTING',
                float(e['score']), e.get('disease_id',''), e['source']
            ), default={})
            bonus = weight*float(e.get('score',0)) if (
                (only or h['case_status'] not in BAD)
                and e.get('inheritance_support') != 'CONFLICTING'
            ) else 0
            return (0 if only else base_score(h) if base_score else float(h['case_score'])) + bonus
        return hypothesis_genes(hypotheses, score)

    if mode != 'standards':
        raise ValueError('overlay mode must be standards or legacy_additive')

    grouped_h = defaultdict(list)
    for h in hypotheses:
        grouped_h[h['gene']].append(h)
    rows = []
    support_rank = {'CONFLICTING': 0, 'UNKNOWN': 1, 'SUPPORTED': 2}
    for gene, hs in grouped_h.items():
        candidates = []
        for h in hs:
            es = [e for e in by_gene[gene]
                  if not e.get('inheritance_model') or e['inheritance_model'] == h['inheritance_model']]
            best = max(es, key=lambda e: (
                support_rank.get(e.get('inheritance_support','UNKNOWN'), 1),
                float(e['score']), e.get('disease_id',''), e.get('source','')
            ), default=None)
            support = best.get('inheritance_support','UNKNOWN') if best else 'UNKNOWN'
            viable = h['case_status'] not in BAD
            case_fit = float(best['score']) if best and viable and support == 'SUPPORTED' else 0.0
            molecular = 0.0 if only else float(base_score(h) if base_score else h['case_score'])
            candidates.append((support_rank[support], case_fit, molecular, h, support, best))
        tier, fit, molecular, h, support, best = min(
            candidates,
            key=lambda x: (-x[0], -x[1], -x[2], x[3]['hypothesis_id'])
        )
        rows.append(dict(
            gene=gene, symbol=h['symbol'], score=round(fit,6),
            molecular_score=round(molecular,6), validity_support=support,
            hypothesis_id=h['hypothesis_id'], inheritance_model=h['inheritance_model'],
            variant_1=h['variant_1'], variant_2=h['variant_2'],
            status='CONFLICTING_EVIDENCE' if support == 'CONFLICTING' else h['case_status'],
            disease_id=best.get('disease_id','') if best else ''
        ))
    rows.sort(key=lambda r: (
        -support_rank.get(r['validity_support'],1),
        -r['score'], -r['molecular_score'], r['gene']
    ))
    i = 0
    while i < len(rows):
        key = (rows[i]['validity_support'], rows[i]['score'], rows[i]['molecular_score'])
        j = i + 1
        while j < len(rows) and (
            rows[j]['validity_support'], rows[j]['score'], rows[j]['molecular_score']
        ) == key:
            j += 1
        for k in range(i,j):
            rows[k]['rank'] = k + 1
            rows[k]['tie_rank_min'] = i + 1
            rows[k]['tie_rank_max'] = j
        i = j
    return rows

def phenotype_matrix(phenotype, ontology, annotations):
    """Independent term-level decomposition; validated against frozen G1 evidence."""
    from functools import lru_cache
    parents, aliases, labels = defaultdict(set), {}, {}
    terms, obsolete, current = set(), set(), None
    for raw in Path(ontology).read_text().splitlines():
        line = raw.strip()
        if line.startswith('['): current = None
        elif line.startswith('id: HP:'):
            current = line[4:]; terms.add(current)
        elif current:
            if line.startswith('is_a: '): parents[current].add(line[6:].split()[0])
            elif line.startswith('alt_id: '): aliases[line[8:]] = current
            elif line.startswith('name: '): labels[current] = line[6:]
            elif line == 'is_obsolete: true': obsolete.add(current)
    with Path(phenotype).open() as f:
        patient = sorted({aliases.get(r['hpo_id'], r['hpo_id']) for r in csv.DictReader(f, delimiter='\t')
                          if r.get('role') == 'proband' and r.get('present', 'true').lower() not in ('false', '0', 'no')})
    diseases = defaultdict(set)
    with Path(annotations).open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            t = aliases.get(r['hpo_id'], r['hpo_id'])
            if r['disease_id'] and t in terms and t not in obsolete:
                diseases[r['disease_id']].add(t)
    @lru_cache(None)
    def ancestors(t):
        seen, pending = set(), [t]
        while pending:
            v = pending.pop()
            if v not in seen:
                seen.add(v); pending.extend(parents[v])
        return seen
    counts = Counter()
    for profile in diseases.values():
        counts.update(set().union(*(ancestors(t) for t in profile)))
    def ic(t): return -math.log((counts[t]+1)/(len(diseases)+1))
    matrix, details = {}, {}
    for disease, profile in diseases.items():
        vals, ds = [], []
        profile_anc = set().union(*(ancestors(t) for t in profile))
        total_ic = sum(ic(p) for p in patient if ic(p) > 0)
        for p in patient:
            common = ancestors(p) & profile_anc
            mica = min(common, key=lambda t: (-ic(t), t)) if common else ''
            matched_ic = ic(mica) if mica else 0
            contribution = (min(ic(p), matched_ic) / total_ic) if total_ic > 0 and ic(p) > 0 else 0
            vals.append(contribution)
            ds.append(dict(hpo_id=p, label=labels.get(p, ''), normalized_similarity=contribution,
                           matched_ancestor=mica, ancestor_label=labels.get(mica, ''),
                           patient_ic=ic(p), ancestor_ic=matched_ic,
                           disease_annotation_frequency=counts[p]/len(diseases)))
        matrix[disease], details[disease] = vals, ds
    return patient, matrix, details


def verify_layers(run):
    verify_g0(run/'G0')
    freeze = digest(run/'G0/G0_FROZEN.json')
    folders = [run/'G1'] + sorted(run.glob('G2_mechanism_*'))
    for folder in folders:
        manifest = json.loads((folder/'layer_manifest.json').read_text())
        if manifest['g0_manifest_sha256'] != freeze:
            raise ValueError('Layer belongs to a different G0')
        for name, expected in manifest['files'].items():
            if Path(name).name != name or digest(folder/name) != expected:
                raise ValueError(f'Layer checksum mismatch: {folder/name}')
    return folders


def audit(run, out):
    run, out = Path(run).resolve(), Path(out).resolve()
    if out == run or run in out.parents:
        raise ValueError('Audit output must be outside the frozen run')
    if out.exists():
        raise FileExistsError(out)
    print('Verifying frozen inputs', flush=True)
    layers = verify_layers(run)
    protected = {str(p.relative_to(run)): digest(p) for p in run.rglob('*') if p.is_file()}
    out.mkdir(parents=True)
    for name in ('ranks', 'comparisons', 'source_snapshot'):
        (out/name).mkdir()
    root = Path(__file__).resolve().parents[1]
    source_hashes = {}
    for name in ('inheritance.py', 'ranking.py', 'qc_integrated_v2.py', 'layers.py', 'phenotype.py', 'evidence_audit.py'):
        path = root/'mva_explorer'/name
        source_hashes[name] = digest(path)
        (out/'source_snapshot'/name).write_bytes(path.read_bytes())
    for name in ('run_qc_integrated_v2.py', 'run_g2_mechanism.py', 'run_full_ranking.py', 'run_evidence_audit.py'):
        path = root/'scripts'/name
        source_hashes[name] = digest(path)
        (out/'source_snapshot'/name).write_bytes(path.read_bytes())
    provenance = json.loads((run/'G0/g0_provenance.json').read_text())
    source_matches = {k: source_hashes[k] == v for k, v in provenance['source_code_sha256'].items()}
    if not all(source_matches.values()):
        raise ValueError('Production inheritance/ranking source changed since the run')
    config = inh.ScoreConfig(**json.loads((run/'G0/score_config.json').read_text()))
    stored = read_csv(run/'G0/variant_evidence.csv')
    hs = read_csv(run/'G0/inheritance_hypotheses.csv')
    g0 = read_csv(run/'G0/06_genome_only_rank.csv')
    g1 = read_csv(run/'G1/gene_rank_comparison.csv')
    g1ev = json.loads((run/'G1/layer_evidence.json').read_text())
    g1manifest = json.loads((run/'G1/layer_manifest.json').read_text())
    evmap = {(e['gene'], e['variant_key']): e for e in stored}
    for e in stored: e['parsed_components'] = json.loads(e['components'])
    scenarios, summaries, descriptions = {}, [], {}
    def save(name, ranks, description):
        scenarios[name] = ranks
        descriptions[name] = description
        write_csv(out/'ranks'/f'{name}.csv', ranks)
    def comparison(before, after, name=None):
        name = name or f'{before}__to__{after}'
        rows, summary = compare(scenarios[before], scenarios[after], name)
        write_csv(out/'comparisons'/f'{name}.csv', rows)
        summary.update(before=before, after=after)
        summaries.append(summary)
    def frozen(rows, layer):
        return rank_rows(dict(gene=r['gene'], symbol=r['symbol'],
                              score=float(r['case_score' if layer == 'G0' else f'{layer}_case_score']),
                              production_rank=int(r['rank' if layer == 'G0' else f'{layer}_rank'])) for r in rows)
    save('G0', frozen(g0, 'G0'), 'Frozen G0 scores; harmonized gene-ID tie breaker. Production rank retained.')
    save('G1', frozen(g1, 'G1'), 'Primary pre-mechanism scores; harmonized ties. Production rank retained.')
    for name, rows in [('G0_production', g0), ('G1_production', g1)]:
        layer = name[:2]
        rr = frozen(rows, layer)
        for r in rr: r['rank'] = r['production_rank']
        save(name, sorted(rr, key=lambda r: r['rank']), 'Exact published order, including historical tie breaker.')
    print('Replaying production code against every saved allele and hypothesis', flush=True)
    replay_ev, replay_hs = replay(stored, config, provenance['phase_evidence'])
    actual = {h['hypothesis_id']: h for h in replay_hs}
    fields = ('case_score', 'inheritance_score', 'case_status', 'PM3', 'phase_status', 'compound_het_confirmed',
              'variant_1_score', 'variant_2_score', 'qc_support', 'model_applicability')
    mismatch = []
    def equal(a, b):
        if inh.number(a) is not None and inh.number(b) is not None:
            return abs(float(a)-float(b)) < 1e-8
        return str(a) == str(b)
    for h in hs:
        new = actual.get(h['hypothesis_id'])
        if new is None or any(not equal(new[k], h[k]) for k in fields):
            mismatch.append(h['hypothesis_id'])
    evactual = {(e['gene'], e['variant_key']): e for e in replay_ev}
    ev_mismatch = [key for key, old in evmap.items() if any(not equal(evactual[key][k], old[k])
                    for k in ('variant_score', 'variant_strength', 'qc_support', 'PM2', 'PVS1', 'acmg_status'))]
    if mismatch or ev_mismatch or len(actual) != len(hs):
        raise ValueError(f'Baseline replay failed: {len(mismatch)} hypothesis, {len(ev_mismatch)} allele mismatches')
    replay_ranks = hypothesis_genes(replay_hs, production_first=True)
    if any(a['gene'] != b['gene'] or a['score'] != b['score'] for a,b in zip(replay_ranks, scenarios['G0'])):
        raise ValueError('G0 gene rank replay failed')
    write_json(out/'replay_validation.json', dict(hypotheses_verified=len(hs), alleles_verified=len(stored),
                gene_ranks_verified=len(g0), mismatches=0, production_source_matches=source_matches,
                qc_observation_source='Frozen QC statuses/penalties; no raw-quality remeasurement'))
    del actual, evactual, replay_ev, replay_hs

    # Score build-up is deliberately separated from eligibility and genotype gates.
    chain = [('neutral', (), False), ('consequence', ('consequence',), False),
             ('transcript', ('consequence', 'transcript'), False),
             ('population', ('consequence', 'transcript', 'rarity', 'population_penalty'), False),
             ('base_qc', ('consequence', 'transcript', 'rarity', 'population_penalty', 'qc'), False),
             ('integrated_qc', tuple(stored[0]['parsed_components']), False),
             ('eligibility', tuple(stored[0]['parsed_components']), True)]
    previous = None
    for stage, components, eligible in chain:
        name = 'build_'+stage
        save(name, variant_genes(stored, lambda e, ks=components: max(0, sum(e['parsed_components'][k] for k in ks)), eligible),
             'Diagnostic score build-up on the extracted candidate universe; max allele per gene. '
             + ('Final variant-strength eligibility applied.' if eligible else 'Eligibility gates not yet applied; not a credible candidate ranking.'))
        if previous: comparison(previous, name)
        previous = name
    # In this production run imported assertions are absent; metadata is zero-effect.
    save('clinical_metadata', scenarios['build_eligibility'], 'Cached assertions/predictions and ACMG review flags add no points in this run.')
    comparison('build_eligibility', 'clinical_metadata')
    single = hypothesis_genes(hs, allowed=lambda h: not h['variant_2'])
    save('single_allele_models', single, 'All single-allele/homozygous/de-novo/X/SV models, with production bonuses; pairs excluded.')
    comparison('clinical_metadata', 'single_allele_models')
    comparison('single_allele_models', 'G0', 'compound_pair_handling')
    save('G1_single_allele_models', overlay_genes([h for h in hs if not h['variant_2']],g1ev,g1manifest['weight'],mode='legacy_additive'),
         'G1 recalculated after excluding all pair hypotheses, retaining single-allele models.')
    comparison('G1_single_allele_models','G1','compound_pair_effect_on_G1')
    save('PM3_metadata', scenarios['G0'], 'PM3 review state has no numerical term; all frozen hypotheses NOT_ASSIGNED.')
    comparison('G0', 'PM3_metadata')
    comparison('G0', 'G1', 'phenotype_harmonized_ties')
    comparison('G0_production', 'G1_production', 'phenotype_published_order')
    save('best_after_qc_resort', hypothesis_genes(hs), 'Re-select highest final-score hypothesis after QC hook mutations.')
    comparison('G0', 'best_after_qc_resort', 'qc_hook_order_consistency')

    # Direct ablations keep eligibility/strength unchanged, reselect best hypothesis.
    for component in ('consequence', 'transcript', 'rarity', 'population_penalty', 'qc', 'qc20', 'computational'):
        def score(h, component=component):
            if h['case_status'] in BAD: return 0
            vals = []
            for key in (h['variant_1'], h['variant_2']):
                if key:
                    e = evmap[h['gene'], key]
                    vals.append(max(0, sum(v for k,v in e['parsed_components'].items() if k != component)))
            return max(vals, default=0) + float(h['inheritance_score'])
        name = 'without_points_'+component
        save(name, hypothesis_genes(hs, score), 'Direct score ablation only; preserve production eligibility and strength gates.')
        comparison(name, 'G0', 'direct_effect_'+component)
        save('G1_'+name, overlay_genes(hs,g1ev,g1manifest['weight'],base_score=score,mode='legacy_additive'),
             'G1 after direct component score ablation; retain eligibility and strength gates.')
        comparison('G1_'+name,'G1','direct_effect_'+component+'_on_G1')
    for name, predicate in [('all_inheritance', lambda h: True),
                             ('unresolved_pair', lambda h: bool(h['variant_2']) and h['phase_status'] in ('UNPHASED','PHASED_UNCERTAIN')),
                             ('all_pair', lambda h: bool(h['variant_2']))]:
        key = 'without_bonus_'+name
        save(key, hypothesis_genes(hs, lambda h: float(h['case_score'])-(float(h['inheritance_score']) if predicate(h) else 0)),
             'Remove indicated bonuses and reselect best hypothesis; no change to allele evidence or PM3.')
        comparison(key, 'G0', 'bonus_effect_'+name)
        save('G1_'+key, overlay_genes(hs,g1ev,g1manifest['weight'],base_score=lambda h: float(h['case_score'])-(float(h['inheritance_score']) if predicate(h) else 0),mode='legacy_additive'),
             'Recompute G1 after bonus ablation; no variant evidence or PM3 change.')
        comparison('G1_'+key,'G1','bonus_effect_'+name+'_on_G1')

    # Gate-aware counterfactuals regenerate previously suppressed hypotheses too.
    for mode in ('no_integrated_qc', 'no_technical_qc', 'population_masked', 'transcript_masked', 'clinical_masked'):
        print('Counterfactual: '+mode, flush=True)
        _, alternative = replay(stored, config, provenance['phase_evidence'], mode)
        save(mode, hypothesis_genes(alternative, production_first=True),
             'Full production hypothesis regeneration. ' + {
                 'no_integrated_qc':'Base QC retained; integrated QC penalties and pair restrictions removed.',
                 'no_technical_qc':'Permissive diagnostic upper bound: technical QC support forced PASS; integrated QC removed. Not a proposed operating mode.',
                 'population_masked':'Observed AF masked to UNKNOWN, never zero; removes frequency gates/rarity bonuses/strong rarity support.',
                 'transcript_masked':'Transcript support masked to UNKNOWN, affecting both +3 and STRONG eligibility.',
                 'clinical_masked':'Imported assertions, calibrated review, cached clinical text and predictors removed.'}[mode])
        comparison(mode, 'G0', 'total_effect_'+mode)
        save('G1_'+mode,overlay_genes(alternative,g1ev,g1manifest['weight'],mode='legacy_additive'),
             'Recompute G1 from fully regenerated counterfactual hypotheses.')
        comparison('G1_'+mode,'G1','total_effect_'+mode+'_on_G1')
        del alternative

    # Single components alone: fixed extracted allele universe, no final gates.
    for component in ('consequence', 'transcript', 'rarity', 'qc', 'qc20'):
        name = 'only_'+component
        save(name, variant_genes(stored, lambda e, k=component: e['parsed_components'][k]),
             'Single component only, no final eligibility; extensive ties make exact overlap arbitrary.')
        comparison(name, 'G0', 'sufficiency_'+component+'_G0')
        comparison(name, 'G1', 'sufficiency_'+component+'_G1')

    grouped_h = defaultdict(list)
    for h in hs: grouped_h[h['gene']].append(h)
    def phenotype_ranks(evidence, weight=10, only=False):
        return overlay_genes(hs,evidence,weight,only,mode='legacy_additive')
    recalculated = phenotype_ranks(g1ev, g1manifest['weight'])
    if any(a['gene'] != b['gene'] or abs(a['score']-b['score']) > 1e-6 for a,b in zip(recalculated,scenarios['G1'])):
        raise ValueError('G1 overlay replay failed')
    save('only_phenotype', phenotype_ranks(g1ev, g1manifest['weight'], only=True), 'Phenotype alone, with no variant viability gating. Does not rescue an allele.')
    comparison('only_phenotype', 'G1', 'sufficiency_phenotype')
    resources = g1manifest['provenance']['resource_sha256']
    paths = {}
    for filename, expected in resources.items():
        path = Path(filename)
        if not path.exists() or digest(path) != expected:
            raise ValueError('Phenotype resource differs from frozen run: '+filename)
        paths[path.name] = path
    print('Decomposing phenotype into individual patient terms', flush=True)
    patient, matrix, details = phenotype_matrix(paths['case_phenotype_G1.tsv'], paths['hp.obo'], paths['genes_to_phenotype.txt'])
    for e in g1ev:
        if abs(sum(matrix[e['disease_id']]) - float(e['score'])) > 1e-10:
            raise ValueError('Term-level phenotype decomposition failed')
    term_rows = []
    best_e = {}
    for e in g1ev:
        if e['gene'] not in best_e or (float(e['score']),e['disease_id']) > (float(best_e[e['gene']]['score']),best_e[e['gene']]['disease_id']):
            best_e[e['gene']] = e
    top_genes = {r['gene'] for r in scenarios['G1'][:50]}
    for gene in sorted(top_genes & best_e.keys()):
        e = best_e[gene]
        for detail in details[e['disease_id']]:
            term_rows.append(dict(gene=gene, symbol=grouped_h[gene][0]['symbol'], disease_id=e['disease_id'],
                                  contribution_to_similarity=detail['normalized_similarity'],
                                  contribution_to_score=detail['normalized_similarity']*g1manifest['weight'], **detail))
    write_csv(out/'phenotype_top50_term_details.csv', term_rows)
    for i, term in enumerate(patient):
        alternative = [dict(e, score=(sum(matrix[e['disease_id']])-matrix[e['disease_id']][i])/(len(patient)-1)) for e in g1ev]
        name = 'phenotype_without_'+term.replace(':','_')
        save(name, phenotype_ranks(alternative, g1manifest['weight']), 'Leave one present proband HPO term out; recompute best disease and mean over remaining terms.')
        comparison(name, 'G1', 'term_effect_'+term.replace(':','_'))

    if len(layers) > 1:
        g2folder = layers[-1]
        g2 = read_csv(g2folder/'gene_rank_comparison.csv')
        save('G2', frozen(g2, 'G2'), 'Experimental mechanism scores for inspection only. G1 remains primary.')
        rr = frozen(g2, 'G2')
        for r in rr: r['rank'] = r['production_rank']
        save('G2_production', sorted(rr, key=lambda r:r['rank']), 'Exact published mechanism order, including tie breaker.')
        comparison('G1', 'G2', 'mechanism_harmonized_ties')
        comparison('G1_production', 'G2_production', 'mechanism_published_order')
        g1map = {r['gene']:r for r in g1}
        inspection = []
        for r in g2:
            pre = g1map[r['gene']]
            inspection.append(dict(gene=r['gene'], symbol=r['symbol'], primary_score=float(pre['G1_case_score']),
                                   primary_score_layer='G1', G1_rank=int(pre['G1_rank']),
                                   mechanism_inspection_rank=int(r['G2_rank']),
                                   experimental_G2_score=float(r['G2_case_score']),
                                   mechanism_increment=round(float(r['G2_case_score'])-float(pre['G1_case_score']),6)))
        write_csv(out/'mechanism_inspection_G1_primary.csv', inspection)

    # Diagnostic tables preserve unsupported calls and pairs, not just winners.
    pair_rows = []
    for h in hs:
        if h['variant_2']:
            pair_rows.append({k:h[k] for k in ('gene','symbol','hypothesis_id','variant_1','variant_2','variant_1_score','variant_2_score',
                             'case_score','case_status','phase_status','inheritance_score','PM3','compound_het_confirmed',
                             'second_allele_strength','population_support','pair_qc20_status')})
    write_csv(out/'all_pair_evidence.csv', pair_rows)
    ledger = []
    best_h = {r['hypothesis_id']:r for r in hs}
    for g in g0:
        h = best_h[g['hypothesis_id']]
        for slot in ('variant_1','variant_2'):
            if not h[slot]: continue
            e = evmap[g['gene'],h[slot]]
            ledger.append(dict(gene=g['gene'],symbol=g['symbol'],G0_rank=int(g['rank']),allele_slot=slot,
                               variant_key=h[slot],variant_score=float(e['variant_score']),variant_strength=e['variant_strength'],
                               observed_af=e['observed_af'],population_status=e['population_status'],consequence_terms=e['consequence'],
                               transcript_relevance=e['transcript_relevance'],qc20_status=e['qc20_status'],
                               inheritance_score=float(h['inheritance_score']),G0_score=float(h['case_score']),
                               PM3=h['PM3'],**e['parsed_components']))
    write_csv(out/'G0_winning_hypothesis_components.csv',ledger)
    expected_status = lambda h: ('HIGH_PRIORITY' if float(h['case_score']) >= config.high_priority else
                                  'MODERATE_PRIORITY' if float(h['case_score']) >= config.moderate_priority else 'EXPLORATORY')
    stale = [h for h in hs if h['case_status'] not in BAD and h['case_status'] != expected_status(h)]
    write_csv(out/'stale_hypothesis_statuses.csv',stale)
    write_json(out/'summary.json',summaries)
    write_json(out/'scenario_definitions.json',descriptions)
    counts = dict(records_scanned=provenance['counts']['records_scanned'],allele_gene_rows=len(stored),genes=len(g0),hypotheses=len(hs),
                  variant_states={k:dict(Counter(e[k] for e in stored)) for k in ('variant_strength','population_status','transcript_relevance','qc_support','qc20_status','acmg_status','PM2','PVS1','PP3_BP4')},
                  hypothesis_states={k:dict(Counter(h[k] for h in hs)) for k in ('PM3','phase_status','case_status')},
                  bonus_pairs=sum(float(h['inheritance_score'])>0 for h in pair_rows),
                  unknown_population_bonus_pairs=sum(float(h['inheritance_score'])>0 and 'UNKNOWN:' in h['population_support'] for h in pair_rows),
                  positive_genes=sum(r['score']>0 for r in scenarios['G0']),stale_statuses=len(stale),
                  phenotype_resource_verified=True,phenotype_terms=patient,
                  raw_predictor_nonempty=sum(any(json.loads(e['input_evidence'])['predictors'].values()) for e in stored),
                  cached_clinical_nonempty=sum(bool(e['cached_clinical_context']) for e in stored),
                  reviewed_assertion_nonempty=sum(bool(json.loads(e['input_evidence'])['reviewed_acmg']) for e in stored),
                  computational_review_nonempty=sum(bool(json.loads(e['input_evidence'])['computational_review']) for e in stored))
    counts['qc_reason_counts'] = dict(Counter(reason for e in stored for reason in e['qc20_reasons'].split(';')))
    counts['base_and_integrated_qc_penalties'] = sum(e['parsed_components']['qc']<0 and e['parsed_components']['qc20']<0 for e in stored)
    counts['supported_af_buckets'] = dict(Counter('UNKNOWN' if e['observed_af']=='' else 'RARE' if float(e['observed_af'])<=config.rare_af else 'BETWEEN_RARE_COMMON'
                                                for e in stored if e['variant_strength']!='UNSUPPORTED'))
    gates = []
    for e in stored:
        v = inh.Variant(**json.loads(e['input_evidence']))
        base = inh.variant_evidence(v,config)
        terms = set(v.consequence.replace(',', '&').split('&'))
        reasons = []
        if v.af is not None and v.af>config.common_af: reasons.append('COMMON_AF')
        if not (terms & (inh.LOF|inh.PROTEIN) or any('splice' in t for t in terms) or v.is_sv): reasons.append('NO_SUPPORTED_MOLECULAR_EFFECT')
        if base['qc_support']=='FAIL': reasons.append('BASE_QC_FAIL')
        if e['qc20_status']=='QC_FAIL': reasons.append('INTEGRATED_QC_FAIL')
        if float(e['variant_score'])<=0: reasons.append('NONPOSITIVE_SCORE')
        if v.transcript_relevance=='UNSUPPORTED': reasons.append('UNSUPPORTED_TRANSCRIPT')
        if v.population_status in ('LOW_CALLABILITY','UNRELIABLE_LOCUS'): reasons.append('UNRELIABLE_POPULATION')
        gates.append(dict(gene=e['gene'],symbol=e['symbol'],variant_key=e['variant_key'],
                          variant_strength=e['variant_strength'],base_qc=base['qc_support'],qc20=e['qc20_status'],
                          observed_af=e['observed_af'],population_status=e['population_status'],
                          overlapping_gate_reasons=';'.join(reasons)))
    write_csv(out/'allele_gate_diagnostics.csv',gates)
    counts['overlapping_gate_reason_counts'] = dict(Counter(reason for r in gates for reason in r['overlapping_gate_reasons'].split(';') if reason))
    write_json(out/'counts.json',counts)
    if any(digest(run/name) != expected for name,expected in protected.items()):
        raise ValueError('Frozen input changed during audit')
    write_json(out/'audit_manifest.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),input_run=str(run),
               protected_input_sha256=protected,inputs_unchanged=True,source_sha256=source_hashes,
               numerical_weights_changed=False,clinical_classification_performed=False,
               limitation='Candidate-universe replay. Does not re-read raw genome, remeasure QC, or recover upstream-excluded/noncandidate alleles.',
               files={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file()}))
    print('Audit complete: '+str(out),flush=True)
    return out

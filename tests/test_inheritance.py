from dataclasses import replace
import ast
import re
from pathlib import Path

import pytest

from mva_explorer.inheritance import (
    MODELS, ScoreConfig, Variant, evaluate_gene, phase_pair, variant_evidence,
)


def allele(pos=100, **kw):
    defaults = dict(gene='SYNTHETIC_GENE', symbol='SYNTHETIC', chrom='3', pos=pos,
                    ref='A', alt='T', gt='0/1', consequence='frameshift_variant',
                    transcript_relevance='SUPPORTED', af=0.00001,
                    dp=30, gq=60, ab=0.5, filter='PASS')
    defaults.update(kw)
    return Variant(**defaults)


def pairs(variants, **kwargs):
    return [h for h in evaluate_gene(variants, **kwargs)[1]
            if h['inheritance_model'] == 'AR_POSSIBLE_COMPOUND_HET' and h['variant_2']]


def test_unphased_strong_and_unresolved_second_generic():
    a = allele(af=None)
    b = allele(200, af=None, consequence='missense_variant')
    ev, hs = evaluate_gene([a, b])
    pair = pairs([a, b])[0]
    assert ev[0]['variant_score'] > ev[1]['variant_score']
    assert all(e['population_status'] == 'UNKNOWN' and e['PM2'] == 'NOT_ASSIGNED' for e in ev)
    assert pair['phase_status'] == 'UNPHASED'
    assert pair['PM3'] == 'NOT_ASSIGNED'
    assert pair['compound_het_confirmed'] is False
    assert pair['second_allele_strength'] == 'SECOND_ALLELE_POSSIBLE'
    assert pair['case_status'] == 'HIGH_PRIORITY'
    assert 'PHASE UNRESOLVED' in pair['rationale']
    assert set(h['inheritance_model'] for h in hs) == set(MODELS)


@pytest.mark.parametrize('changes', [dict(gq=2), dict(dp=2), dict(ab=.98)])
def test_generic_numeric_qc_defers_pair_until_technical_review(changes):
    variants = [allele(), allele(200, **changes)]
    ev, hs = evaluate_gene(variants)
    assert next(e for e in ev if e['variant_key'] == variants[1].key)['qc_support'] == 'REVIEW'
    assert not [h for h in hs if h['inheritance_model'] == 'AR_POSSIBLE_COMPOUND_HET' and h['variant_2']]


def test_hard_fail_second_does_not_form_pair():
    variants = [allele(), allele(200, filter='LowQual')]
    ev, hs = evaluate_gene(variants)
    assert next(e for e in ev if e['variant_key'] == variants[1].key)['qc_support'] == 'FAIL'
    assert not [h for h in hs if h['inheritance_model'] == 'AR_POSSIBLE_COMPOUND_HET' and h['variant_2']]


@pytest.mark.parametrize('changes', [dict(consequence='synonymous_variant'),
                                      dict(transcript_relevance='UNSUPPORTED'),
                                      dict(population_status='LOW_CALLABILITY')])
def test_biologically_unsupported_second_no_bonus(changes):
    pair = pairs([allele(), allele(200, **changes)])[0]
    assert pair['second_allele_strength'] == 'SECOND_ALLELE_UNSUPPORTED'
    assert pair['inheritance_score'] == 0
    assert pair['case_status'] == 'INSUFFICIENT_EVIDENCE'


def test_multiple_weak_alleles_cannot_boost_gene():
    strong = allele()
    original = evaluate_gene([strong])[1][0]['case_score']
    variants = [strong] + [allele(p, consequence='synonymous_variant') for p in range(200, 400, 10)]
    assert evaluate_gene(variants)[1][0]['case_score'] == original


def test_unresolved_pair_rarity_does_not_create_bonus():
    strong = pairs([allele(), allele(200)])[0]
    possible = pairs([allele(), allele(200, consequence='missense_variant')])[0]
    assert strong['inheritance_score'] == 0
    assert possible['inheritance_score'] == 0


@pytest.mark.parametrize('gts,expected', [(('0|1', '1|0'), 'IN_TRANS_CONFIRMED'),
                                        (('0|1', '0|1'), 'IN_CIS_CONFIRMED')])
def test_phase_same_block(gts, expected):
    a, b = allele(gt=gts[0], ps='123'), allele(200, gt=gts[1], ps='123')
    pair = pairs([a, b])[0]
    assert pair['phase_status'] == expected
    if expected == 'IN_CIS_CONFIRMED':
        assert pair['inheritance_score'] == 0
        assert not pair['compound_het_confirmed']
        assert pair['PM3'] == 'NOT_ASSIGNED'
    else:
        assert pair['compound_het_confirmed']
        # Confirmed trans alone is not PM3: an independently reviewed P/LP
        # partner is required before PM3 enters review.
        assert pair['PM3'] == 'NOT_ASSIGNED'


@pytest.mark.parametrize('changes', [dict(ps='other'), dict(ps=''), dict(ps='.')])
def test_phase_blocks_not_comparable(changes):
    a = allele(gt='0|1', ps='123')
    b = allele(200, gt='1|0', **changes)
    assert phase_pair(a, b)[0] == 'PHASED_UNCERTAIN'
    assert pairs([a, b])[0]['PM3'] == 'NOT_ASSIGNED'


def test_pgt_pid_and_conflicts():
    a = allele(pgt='0|1', pid='block')
    b = allele(200, pgt='1|0', pid='block')
    assert phase_pair(a, b)[0] == 'IN_TRANS_CONFIRMED'
    a.gt = b.gt = '0|1'; a.ps = b.ps = '123'
    assert phase_pair(a, b)[0] == 'PHASED_UNCERTAIN'
    pair = pairs([a, b])[0]
    assert pair['case_status'] == 'CONFLICTING_EVIDENCE'
    assert pair['PM3'] == 'NOT_ASSIGNED'


def parent(gt, **kw):
    return dict(gt=gt, dp=30, gq=50, filter='PASS', source='synthetic pedigree', **kw)


def test_parental_origin_and_de_novo():
    a = allele(parents={'mother': parent('0/1'), 'father': parent('0/0')})
    b = allele(200, parents={'mother': parent('0/0'), 'father': parent('0/1')})
    assert pairs([a, b])[0]['phase_status'] == 'IN_TRANS_CONFIRMED'
    v = allele(parents={'mother': parent('0/0'), 'father': parent('0/0')})
    h = next(h for h in evaluate_gene([v])[1] if h['inheritance_model'] == 'DE_NOVO_COMPATIBLE')
    assert h['inheritance_score'] == 0
    v.parents['father']['gq'] = 1
    h = next(h for h in evaluate_gene([v])[1] if h['inheritance_model'] == 'DE_NOVO_COMPATIBLE')
    assert h['inheritance_score'] == 0
    assert h['case_status'] == 'INSUFFICIENT_EVIDENCE'


def test_external_phase_requires_provenance_validation_and_exact_pair():
    a, b = allele(), allele(200)
    e = dict(variants=[a.key, b.key], status='IN_TRANS_CONFIRMED', source='validated read phasing', validated=True)
    assert phase_pair(a, b, external=[e])[0] == 'IN_TRANS_CONFIRMED'
    assert phase_pair(a, b, external=[dict(e, validated=False)])[0] == 'UNPHASED'
    assert phase_pair(a, b, external=[dict(e, source='')])[0] == 'UNPHASED'
    assert phase_pair(a, b, external=[dict(e, variants=[a.key, 'different'])])[0] == 'UNPHASED'


def test_duplicate_alleles_not_pairs_and_conflicts_not_hidden():
    a = allele()
    ev, _ = evaluate_gene([a, replace(a)])
    assert len(ev) == 1
    assert not pairs([a, replace(a)])
    ev, _ = evaluate_gene([a, replace(a, gq=3)])
    assert ev[0]['qc_support'] == 'FAIL'
    # Transcript evidence cannot be moved from a mild transcript onto a severe one.
    ev, _ = evaluate_gene([replace(a, transcript_relevance='UNKNOWN'), replace(a, consequence='synonymous_variant')])
    assert ev[0]['transcript_relevance'] == 'UNKNOWN'


def test_overlapping_alleles_need_haplotype_review():
    p = pairs([allele(ref='AT', alt='A'), allele(101)])[0]
    assert p['inheritance_score'] == 0
    assert 'overlapping' in p['rationale']


def test_homozygous_gated_and_variant_score_zygosity_independent():
    a = allele()
    hom = replace(a, gt='1/1', ab=.98)
    assert variant_evidence(a)['variant_score'] == variant_evidence(hom)['variant_score']
    get = lambda v: next(h for h in evaluate_gene([v])[1] if h['inheritance_model'] == 'AR_HOMOZYGOUS')
    assert get(hom)['inheritance_score'] == 0
    assert get(replace(hom, af=.4))['inheritance_score'] == 0
    assert get(replace(hom, ab=.5))['inheritance_score'] == 0
    assert get(replace(hom, consequence='synonymous_variant'))['inheritance_score'] == 0
    # Rarity/absence does not create pathogenicity support by itself.
    assert get(replace(hom, af=None))['inheritance_score'] == get(hom)['inheritance_score']


def test_dominant_single_allele_not_penalized():
    a = allele()
    h = next(h for h in evaluate_gene([a])[1] if h['inheritance_model'] == 'AD_SINGLE_ALLELE')
    assert h['case_score'] == variant_evidence(a)['variant_score']
    assert h['variant_2'] == ''
    assert h['gene_disease_inheritance_support'] == 'UNKNOWN_G0_NOT_QUERIED'


def test_x_sv_and_haploid_other():
    x = allele(chrom='X', gt='1', ab=1)
    hs = evaluate_gene([x])[1]
    assert next(h for h in hs if h['inheritance_model'] == 'X_LINKED')['case_score'] > 0
    assert next(h for h in hs if h['inheritance_model'] == 'AR_HOMOZYGOUS')['case_score'] == 0
    sv = allele(alt='<DEL>', consequence='transcript_ablation', ab=None, end=500)
    hs = evaluate_gene([sv])[1]
    assert next(h for h in hs if h['inheritance_model'] == 'CNV_SV_COMPATIBLE')['case_score'] == 0
    assert variant_evidence(sv)['standards_route'] == 'ACMG_CLINGEN_CNV_SV_SEPARATE'
    assert not pairs([sv, allele(200, gt='1')])[0]['compound_het_confirmed']


def test_missing_af_never_absent_or_zero():
    e = variant_evidence(allele(af=None))
    assert e['observed_af'] is None and e['population_status'] == 'UNKNOWN'
    assert e['PM2'] == 'NOT_ASSIGNED'
    v = allele(af=None, population_status='NOT_OBSERVED_IN_COMPLETE_RESOURCE')
    assert v.population_status == 'UNKNOWN'
    v = allele(af=None, population_status='NOT_OBSERVED_IN_COMPLETE_RESOURCE',
               population_source='complete synthetic resource v1', complete_callable_resource=True)
    assert variant_evidence(v)['PM2'] == 'NOT_ASSIGNED'
    assert v.af is None
    assert variant_evidence(allele(af=0))['population_status'] == 'OBSERVED'


def test_computational_predictions_not_independent_and_cached_labels_not_acmg():
    a = allele(consequence='missense_variant')
    b = replace(a, predictors={'SIFT': 'deleterious', 'PolyPhen': 'probably_damaging', 'REVEL': .99},
                cached_clinical='pathogenic')
    assert variant_evidence(a)['variant_score'] == variant_evidence(b)['variant_score']
    assert variant_evidence(b)['PP3_BP4'] == 'NOT_ASSIGNED'
    assert variant_evidence(b)['acmg_status'] == 'NOT_CLASSIFIED'
    b.computational_review = dict(validated=True, criterion='PP3', tool='synthetic calibrated tool',
                                 rule_id='rule1', version='1', source='synthetic review')
    assert variant_evidence(b)['PP3_BP4'] == 'PP3_REVIEW_REQUIRED'
    assert variant_evidence(b)['components']['computational'] == 0



def test_rarity_is_context_not_positive_score():
    rare = variant_evidence(allele(af=0.00001))
    unknown = variant_evidence(allele(af=None))
    assert rare['components']['rarity'] == 0
    assert rare['variant_score'] == unknown['variant_score']
    assert rare['population_context'] == 'LOW_AF_CONTEXT_ONLY'
    assert unknown['population_context'] == 'UNKNOWN_OR_INTERMEDIATE_CONTEXT'
    assert rare['population_evidence_semantics'] == 'CONTEXT_ONLY_NOT_PATHOGENICITY_IN_G0'


def test_pm3_requires_confirmed_trans_and_reviewed_pathogenic_partner():
    assertion = dict(classification='PATHOGENIC', validated=True, criteria_reviewed=True,
                     source='synthetic expert review', version='1', disease_id='SYN:1')
    a = allele(gt='0|1', ps='123', reviewed_acmg=assertion)
    b = allele(200, gt='1|0', ps='123')
    pair = pairs([a, b])[0]
    assert pair['phase_status'] == 'IN_TRANS_CONFIRMED'
    assert pair['PM3'] == 'NOT_ASSIGNED'

def test_gene_and_coordinate_renaming_invariance():
    baseline = pairs([allele(af=None), allele(200, af=None, consequence='missense_variant')])[0]
    for gene, offset in [('BUB1B', 0), ('UNRELATED_SYNTHETIC', 1000000)]:
        renamed = pairs([allele(100 + offset, gene=gene, af=None),
                         allele(200 + offset, gene=gene, af=None, consequence='missense_variant')])[0]
        for key in ('case_score', 'inheritance_score', 'case_status', 'phase_status', 'PM3',
                    'compound_het_confirmed', 'second_allele_strength', 'rationale'):
            assert renamed[key] == baseline[key]


def test_config_and_input_validation():
    with pytest.raises(ValueError): ScoreConfig(rare_af=.5, common_af=.1)
    with pytest.raises(ValueError): ScoreConfig(strong_pair=float('nan'))
    with pytest.raises(ValueError): allele(af=2)
    with pytest.raises(ValueError): evaluate_gene([allele(), allele(200, gene='OTHER')])


def audit_ranking_source(source):
    tree = ast.parse(source)
    forbidden = {'BUB1B', 'TRIP13', 'CEP57', 'BUB1', 'MAD2L1'}
    for node in ast.walk(tree):
        # Stronger than conditionals: no challenge gene literals anywhere in production scoring.
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            assert not any(re.search(r'\b' + symbol + r'\b', node.value.upper()) for symbol in forbidden), 'challenge gene literal'
            assert not re.search(r'(?<![0-9])[0-9]{7,}(?![0-9])', node.value), 'coordinate-like string literal in scoring source'
        if isinstance(node, (ast.If, ast.IfExp, ast.Compare)):
            expression = node.test if hasattr(node, 'test') else node
            for child in ast.walk(expression):
                if isinstance(child, ast.Constant) and isinstance(child.value, (int, float)):
                    assert abs(child.value) < 1000000, 'challenge-sized coordinate literal in ranking conditional'


def test_no_challenge_specific_production_logic():
    root = Path(__file__).resolve().parents[1]
    for name in ('inheritance.py', 'ranking.py', 'layers.py', 'phenotype.py',
                 'exomiser_compatibility.py', 'compatibility_view.py'):
        audit_ranking_source((root / 'mva_explorer' / name).read_text())
    for name in ('provisional_snapshot_rank.py', 'g1_hpo_rerank.py', 'g2_rerank.py',
                 'run_current_exomiser_compatibility.py'):
        audit_ranking_source((root / 'scripts' / name).read_text())
    for bad in ('if gene == "BUB1B": score += 10', 'if pos == 40456789: score += 10',
                'if key == "GRCh38:15:40456789:A:T": score += 10'):
        with pytest.raises(AssertionError): audit_ranking_source(bad)


def test_cnv_dosage_without_diploid_snv_genotype():
    cnv = allele(alt='<DEL>', gt='', ab=None, dp=None, gq=None, variant_class='CNV',
                 copy_number=1, expected_copy_number=2, structural_qc='PASS')
    ev, hs = evaluate_gene([cnv])
    assert ev[0]['qc_support'] == 'PASS'
    assert next(h for h in hs if h['inheritance_model'] == 'CNV_SV_COMPATIBLE')['case_score'] == 0
    assert ev[0]['standards_route'] == 'ACMG_CLINGEN_CNV_SV_SEPARATE'
    assert all(h['case_status'] == 'INSUFFICIENT_EVIDENCE' and h['inheritance_score'] == 0
               for h in pairs([cnv, allele(200)]))


def test_second_allele_above_generic_rare_threshold_is_context_only():
    pair = pairs([allele(), allele(200, af=.005, consequence='missense_variant')])[0]
    assert pair['inheritance_score'] == 0
    # G0 has no disease-specific maximum credible AF, so 0.005 is not
    # automatically pathogenic or disqualifying; it remains contextual.
    assert pair['second_allele_strength'] == 'SECOND_ALLELE_POSSIBLE'


def test_reviewed_pathogenic_missense_can_be_strong_without_cached_label_shortcut():
    assertion = dict(classification='PATHOGENIC', validated=True, criteria_reviewed=True,
                     disease_id='SYN:1', source='synthetic reviewed assertion', version='1')
    v = allele(consequence='missense_variant', reviewed_acmg=assertion)
    e = variant_evidence(v)
    assert e['variant_strength'] == 'POSSIBLE'
    assert e['acmg_status'] == 'PATHOGENIC'
    assert e['acmg_origin'] == 'IMPORTED_REVIEWED_ASSERTION'
    assert e['components']['consequence'] == ScoreConfig().protein
    assert ScoreConfig().reviewed_pathogenic == 0
    v.reviewed_acmg = dict(assertion, classification='BENIGN')
    assert variant_evidence(v)['variant_strength'] == 'POSSIBLE'
    v.reviewed_acmg = dict(assertion, validated=False)
    assert variant_evidence(v)['acmg_status'] == 'NOT_CLASSIFIED'



def test_noncanonical_splice_requires_prediction_or_rna_support():
    v = allele(
        consequence='splice_region_variant',
        transcript_relevance='SUPPORTED',
        transcript_biotype='protein_coding',
        selected_transcript='TX1',
        transcript_selection='MANE_SELECT',
    )
    e = variant_evidence(v)
    assert e['variant_score'] == 0
    assert e['variant_strength'] == 'UNSUPPORTED'
    assert e['standards_route'] == 'SPLICE_PREDICTION_RNA_REVIEW'
    assert e['splice_review_status'] == 'PREDICTION_OR_RNA_EVIDENCE_REQUIRED'
    assert e['PVS1'] == 'NOT_ASSIGNED'


def test_noncoding_transcript_not_scored_by_coding_heuristic():
    v = allele(
        consequence='splice_region_variant',
        transcript_relevance='SUPPORTED',
        transcript_biotype='lncRNA',
        selected_transcript='TXNC',
        transcript_selection='MANE_SELECT',
    )
    e = variant_evidence(v)
    assert e['variant_score'] == 0
    assert e['variant_strength'] == 'UNSUPPORTED'
    assert e['standards_route'] == 'NONCODING_TRANSCRIPT_SEPARATE'


def test_canonical_splice_donor_remains_discoverable_without_pvs1():
    v = allele(
        consequence='splice_donor_variant',
        transcript_relevance='SUPPORTED',
        transcript_biotype='protein_coding',
        selected_transcript='TX1',
        transcript_selection='MANE_SELECT',
    )
    e = variant_evidence(v)
    assert e['variant_score'] == ScoreConfig().lof
    assert e['variant_strength'] == 'POSSIBLE'
    assert e['standards_route'] == 'SEQUENCE_VARIANT_RESEARCH_TRIAGE'
    assert e['PVS1'] == 'NOT_ASSIGNED'

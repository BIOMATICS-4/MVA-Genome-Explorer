"""Gene-agnostic genome evidence and inheritance hypotheses (research heuristics).

No phenotype, disease whitelist, or automatic ACMG classification belongs here.
See docs/inheritance_engine.md for score definitions and deliberate limitations.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from hashlib import sha256
from itertools import combinations
import json
import math
import re

MODELS = ('AR_HOMOZYGOUS', 'AR_POSSIBLE_COMPOUND_HET', 'AD_SINGLE_ALLELE',
          'DE_NOVO_COMPATIBLE', 'X_LINKED', 'CNV_SV_COMPATIBLE', 'UNRESTRICTED_OTHER')
PHASES = {'IN_TRANS_CONFIRMED', 'IN_CIS_CONFIRMED', 'PHASED_UNCERTAIN',
          'UNPHASED', 'NOT_APPLICABLE'}
POPULATIONS = {'OBSERVED', 'NOT_OBSERVED_IN_COMPLETE_RESOURCE', 'UNKNOWN',
               'LOW_CALLABILITY', 'UNRELIABLE_LOCUS'}
LOF = {'transcript_ablation', 'stop_gained', 'frameshift_variant',
       'splice_acceptor_variant', 'splice_donor_variant'}
PROTEIN = {'missense_variant', 'inframe_insertion', 'inframe_deletion',
           'protein_altering_variant', 'start_lost', 'stop_lost'}
NONCANONICAL_SPLICE = {
    'splice_region_variant', 'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant', 'splice_donor_5th_base_variant',
}
REGULATORY = {'regulatory_region_ablation', 'regulatory_region_amplification', 'TFBS_ablation'}


@dataclass(frozen=True)
class ScoreConfig:
    version: str = 'inheritance-2.1-transcript-specific-research'
    rare_af: float = 0.001
    common_af: float = 0.01
    min_dp: float = 10
    min_gq: float = 20
    het_ab_min: float = 0.2
    het_ab_max: float = 0.8
    hom_ab_min: float = 0.85
    lof: float = 4
    protein: float = 3
    splice_region: float = 0
    other: float = 0
    reviewed_pathogenic: float = 0
    rarity: float = 0
    transcript: float = 0
    qc_failure_penalty: float = 0
    common_penalty: float = 0
    possible_pair: float = 0
    strong_pair: float = 0
    trans: float = 0
    homozygous: float = 0
    homozygous_unknown_rarity: float = 0
    de_novo: float = 0
    high_priority: float = 4
    moderate_priority: float = 2

    def __post_init__(self):
        for k, v in asdict(self).items():
            if k != 'version' and (not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0):
                raise ValueError(f'Invalid scoring configuration: {k}')
        if not 0 <= self.rare_af <= self.common_af <= 1:
            raise ValueError('Frequency thresholds must be ordered within [0,1]')
        if not 0 <= self.het_ab_min < self.het_ab_max <= 1 or not 0 <= self.hom_ab_min <= 1:
            raise ValueError('Invalid allele balance thresholds')
        if self.high_priority < self.moderate_priority:
            raise ValueError('Priority thresholds must be ordered')


def number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def present(value):
    return value is not None and str(value) not in ('', '.', 'None')


def genotype(gt):
    parts = re.split(r'[/|]', gt or '')
    if not parts or any(p not in ('0', '1') for p in parts):
        return 'MISSING_OR_PARTIAL'
    if '1' not in parts:
        return 'NOT_CARRIED'
    if len(parts) == 1:
        return 'HEMIZYGOUS_OR_HAPLOID'
    if len(parts) != 2:
        return 'OTHER_PLOIDY'
    return 'HOMOZYGOUS_ALT' if parts == ['1', '1'] else 'HETEROZYGOUS'


@dataclass
class Variant:
    gene: str
    chrom: str
    pos: int
    ref: str
    alt: str
    gt: str
    build: str = 'GRCh38'
    symbol: str = ''
    consequence: str = ''
    transcript_relevance: str = 'UNKNOWN'
    af: float | None = None
    population_status: str = 'UNKNOWN'
    population_source: str = ''
    complete_callable_resource: bool = False
    dp: float | None = None
    gq: float | None = None
    ab: float | None = None
    filter: str = '.'
    ps: str = ''
    pid: str = ''
    pgt: str = ''
    variant_class: str = 'SMALL_VARIANT'
    end: int | None = None
    copy_number: float | None = None
    expected_copy_number: float | None = None
    structural_qc: str = 'UNKNOWN'
    transcripts: str = ''
    selected_transcript: str = ''
    transcript_biotype: str = ''
    transcript_selection: str = ''
    alternate_transcript_effects: list = field(default_factory=list)
    cached_clinical: str = ''
    predictors: dict = field(default_factory=dict)
    # Imported reviewed calibration; raw predictors alone never contribute points.
    computational_review: dict = field(default_factory=dict)
    reviewed_acmg: dict = field(default_factory=dict)
    parents: dict = field(default_factory=dict)
    annotations: list = field(default_factory=list)
    qc_conflict: bool = False

    def __post_init__(self):
        self.chrom = str(self.chrom).removeprefix('chr')
        self.pos = int(self.pos)
        self.af = number(self.af)
        for key in ('dp', 'gq', 'ab', 'copy_number', 'expected_copy_number'):
            setattr(self, key, number(getattr(self, key)))
        if self.af is not None and not 0 <= self.af <= 1:
            raise ValueError('AF must be in [0,1] or unknown')
        if self.ab is not None and not 0 <= self.ab <= 1:
            raise ValueError('Allele balance must be in [0,1] or unknown')
        if self.af is not None and self.population_status not in ('LOW_CALLABILITY', 'UNRELIABLE_LOCUS'):
            self.population_status = 'OBSERVED'
        if self.population_status not in POPULATIONS:
            self.population_status = 'UNKNOWN'
        if self.population_status == 'OBSERVED' and self.af is None:
            self.population_status = 'UNKNOWN'
        if self.population_status == 'NOT_OBSERVED_IN_COMPLETE_RESOURCE' and not (
                self.complete_callable_resource and self.population_source and self.af is None):
            self.population_status = 'UNKNOWN'
        if self.alt.startswith('<') or '[' in self.alt or ']' in self.alt:
            self.variant_class = 'SV'

    @property
    def key(self):
        suffix = f':END={self.end}' if self.end is not None else ''
        return f'{self.build}:{self.chrom}:{self.pos}:{self.ref}:{self.alt}{suffix}'

    @property
    def is_sv(self):
        return self.variant_class in ('SV', 'CNV')


def qc_support(v, c):
    """Generic research QC triage.

    Only direct call-integrity conflicts are hard failures here. Numeric DP/GQ/AB
    thresholds are not universal clinical standards; without a validated assay
    profile they trigger technical review rather than exclusion. Structural
    dosage calls use their dedicated structural QC path instead of SNV metrics.
    """
    fail, review = [], []
    dosage_call = (v.is_sv and v.copy_number is not None and v.expected_copy_number is not None
                   and v.copy_number >= 0 and v.expected_copy_number > 0
                   and v.copy_number != v.expected_copy_number)
    if v.qc_conflict:
        fail.append('conflicting duplicate observations')
    if v.filter not in ('PASS', '.', ''):
        fail.append('non-PASS FILTER')
    if v.structural_qc == 'FAIL':
        fail.append('structural call QC failed')
    if genotype(v.gt) in ('MISSING_OR_PARTIAL', 'NOT_CARRIED') and not dosage_call:
        fail.append('genotype incomplete or not carried; no supported dosage change')
    if fail:
        return 'FAIL', fail

    if dosage_call and v.structural_qc == 'PASS':
        return 'PASS', []

    if v.dp is None:
        review.append('depth unknown')
    elif v.dp < c.min_dp:
        review.append('low depth')
    if v.gq is None:
        review.append('genotype quality unknown')
    elif v.gq < c.min_gq:
        review.append('low genotype quality')
    if v.ab is None and not v.is_sv:
        review.append('allele balance unknown')
    elif v.ab is not None and not v.is_sv:
        if genotype(v.gt) == 'HETEROZYGOUS' and not c.het_ab_min <= v.ab <= c.het_ab_max:
            review.append('heterozygous allele balance outside generic review range')
        if genotype(v.gt) == 'HOMOZYGOUS_ALT' and v.ab < c.hom_ab_min:
            review.append('homozygous allele balance outside generic review range')

    if review:
        return 'REVIEW', review
    return 'PASS', []


def variant_evidence(v, c=ScoreConfig()):
    terms = set(re.split('[&,]', v.consequence))
    consequence_score = (c.lof if terms & LOF else c.protein if terms & PROTEIN
                         else c.splice_region if terms & NONCANONICAL_SPLICE else c.other)
    assertion = v.reviewed_acmg
    classification = assertion.get('classification', 'NOT_CLASSIFIED')
    curated = (classification in ('PATHOGENIC', 'LIKELY_PATHOGENIC', 'VUS', 'LIKELY_BENIGN', 'BENIGN')
               and assertion.get('validated') is True and assertion.get('criteria_reviewed') is True
               and all(assertion.get(k) for k in ('source', 'version', 'disease_id')))
    pathogenic = curated and classification in ('PATHOGENIC', 'LIKELY_PATHOGENIC')
    benign = curated and classification in ('BENIGN', 'LIKELY_BENIGN')
    # Imported classifications are context, not independent numeric evidence.
    # Their underlying evidence may overlap consequence/population/segregation,
    # so they never add points to this research triage score.
    relevant = v.transcript_relevance == 'SUPPORTED'
    excluded = v.transcript_relevance == 'UNSUPPORTED'
    if excluded:
        consequence_score = c.other
    rare = ((v.population_status == 'OBSERVED' and v.af is not None and v.af <= c.rare_af)
            or v.population_status == 'NOT_OBSERVED_IN_COMPLETE_RESOURCE')
    common = v.population_status == 'OBSERVED' and v.af is not None and v.af > c.common_af
    qc, flags = qc_support(v, c)
    comp = v.computational_review
    calibrated = (comp.get('validated') is True and comp.get('criterion') in ('PP3', 'BP4')
                  and all(comp.get(k) for k in ('rule_id', 'version', 'source', 'tool')))
    # Calibration results enter the ACMG review queue once, not the heuristic score.
    computational = comp['criterion'] + '_REVIEW_REQUIRED' if calibrated else 'NOT_ASSIGNED'
    # Population rarity is context, not standalone pathogenicity evidence.
    # In phenotype-blind G0 there is no disease-specific maximum credible AF,
    # so rarity/absence contributes no positive heuristic points. Clearly common
    # observations may still contradict a rare-Mendelian hypothesis.
    components = {'consequence': consequence_score, 'rarity': 0,
                  'transcript': c.transcript if relevant and not excluded and consequence_score > 0 else 0,
                  'qc': 0, 'population_penalty': 0, 'computational': 0}
    # G0 is molecular discovery, not condition-specific pathogenicity.
    # Imported clinical assertions and computational criteria remain metadata only.
    canonical_molecular_effect = bool(terms & (LOF | PROTEIN))
    noncanonical_splice = bool(terms & NONCANONICAL_SPLICE) and not bool(terms & LOF)
    plausible_effect = canonical_molecular_effect
    # CNV/SV, non-protein-coding transcripts, regulatory effects and
    # non-canonical splice-region calls require dedicated interpretation.
    separate_sv = v.is_sv
    separate_noncoding = bool(v.transcript_biotype) and v.transcript_biotype != 'protein_coding'
    separate_noncoding |= bool(terms & REGULATORY) and not canonical_molecular_effect
    separate_splice_review = noncanonical_splice and not (
        calibrated and comp.get('criterion') == 'PP3'
    )
    credible = (plausible_effect and qc != 'FAIL' and not excluded
                and not separate_sv and not separate_noncoding and not separate_splice_review
                and v.population_status not in ('LOW_CALLABILITY', 'UNRELIABLE_LOCUS'))
    # Avoid ACMG strength terminology in phenotype-blind discovery.
    strength = 'POSSIBLE' if credible else 'UNSUPPORTED'
    return {'gene': v.gene, 'symbol': v.symbol, 'variant_key': v.key,
            'variant_score': round(max(0, sum(components.values())), 3),
            'score_semantics': 'ORDINAL_MOLECULAR_IMPACT_NOT_PATHOGENICITY',
            'molecular_effect_tier': ('PREDICTED_LOF' if terms & LOF else
                                      'PROTEIN_ALTERING' if terms & PROTEIN else
                                      'NONCANONICAL_SPLICE_REVIEW' if terms & NONCANONICAL_SPLICE else
                                      'OTHER'),
            'variant_strength': strength, 'components': components,
            'population_status': v.population_status, 'observed_af': v.af,
            'population_context': ('LOW_AF_CONTEXT_ONLY' if rare else
                                   'HIGH_AF_CONTEXT_ONLY' if common else
                                   'UNKNOWN_OR_INTERMEDIATE_CONTEXT'),
            'population_evidence_semantics': 'CONTEXT_ONLY_NOT_PATHOGENICITY_IN_G0',
            'qc_support': qc, 'qc_flags': flags, 'transcript_relevance': v.transcript_relevance,
            'consequence': v.consequence, 'genotype': genotype(v.gt),
            'cached_clinical_context': v.cached_clinical,
            'computational_dependency_group': 'COMPUTATIONAL_PREDICTION',
            'computational_status': computational, 'acmg_status': classification if curated else 'NOT_CLASSIFIED',
            'acmg_origin': 'IMPORTED_REVIEWED_ASSERTION' if curated else 'NOT_CLASSIFIED',
            'PM2': 'NOT_ASSIGNED',
            'PVS1': 'NOT_ASSIGNED',
            'lof_review_status': 'MECHANISM_TRANSCRIPT_NMD_REVIEW_REQUIRED' if terms & LOF else 'NOT_APPLICABLE',
            'standards_route': ('ACMG_CLINGEN_CNV_SV_SEPARATE' if separate_sv
                                else 'NONCODING_TRANSCRIPT_SEPARATE' if separate_noncoding
                                else 'SPLICE_PREDICTION_RNA_REVIEW' if separate_splice_review
                                else 'SEQUENCE_VARIANT_RESEARCH_TRIAGE'),
            'selected_transcript': v.selected_transcript,
            'transcript_biotype': v.transcript_biotype,
            'transcript_selection': v.transcript_selection,
            'alternate_transcript_effects': v.alternate_transcript_effects,
            'transcript_review_status': ('ALTERNATE_HIGHER_IMPACT_TRANSCRIPT_REVIEW'
                                         if v.alternate_transcript_effects else 'NOT_APPLICABLE'),
            'splice_review_status': ('PREDICTION_OR_RNA_EVIDENCE_REQUIRED'
                                     if separate_splice_review else 'NOT_APPLICABLE'),
            'BA1_BS1': 'NOT_ASSIGNED',
            'PS2_PM6': 'NOT_ASSIGNED',
            'PS3_BS3': 'NOT_ASSIGNED',
            'PP1_BS4': 'NOT_ASSIGNED',
            'PP4': 'NOT_ASSIGNED',
            'input_evidence': asdict(v),
            'PP3_BP4': computational, 'PM3': 'NOT_ASSIGNED'}


def parent_origin(v, c):
    """Only complete, adequate-quality parental calls can resolve parental origin."""
    calls = []
    for role in ('mother', 'father'):
        p = v.parents.get(role, {})
        gt = p.get('gt', '')
        if (not p.get('source') or number(p.get('dp')) is None or number(p.get('gq')) is None
                or number(p.get('dp')) < c.min_dp or number(p.get('gq')) < c.min_gq
                or p.get('filter') != 'PASS' or genotype(gt) in ('MISSING_OR_PARTIAL', 'OTHER_PLOIDY')):
            return None
        calls.append(genotype(gt) != 'NOT_CARRIED')
    return 'mother' if calls == [True, False] else 'father' if calls == [False, True] else 'neither' if calls == [False, False] else 'both'


def phase_pair(a, b, c=ScoreConfig(), external=()):
    if genotype(a.gt) != 'HETEROZYGOUS' or genotype(b.gt) != 'HETEROZYGOUS':
        return 'NOT_APPLICABLE', 'not two diploid heterozygous calls'
    observations = []
    for e in external:
        if (set(e.get('variants', [])) == {a.key, b.key} and e.get('validated') is True
                and e.get('source') and e.get('status') in ('IN_TRANS_CONFIRMED', 'IN_CIS_CONFIRMED')):
            observations.append((e['status'], e['source']))
    if a.chrom == b.chrom and a.build == b.build:
        for tag, g1, g2, p1, p2 in (
                ('GT/PS', a.gt, b.gt, a.ps, b.ps), ('PGT/PID', a.pgt, b.pgt, a.pid, b.pid)):
            if present(p1) and present(p2) and str(p1) == str(p2) and g1 in ('0|1', '1|0') and g2 in ('0|1', '1|0'):
                observations.append(('IN_CIS_CONFIRMED' if g1 == g2 else 'IN_TRANS_CONFIRMED', tag))
        oa, ob = parent_origin(a, c), parent_origin(b, c)
        if {oa, ob} == {'mother', 'father'}:
            observations.append(('IN_TRANS_CONFIRMED', 'parental genotypes'))
    if observations:
        states = {x[0] for x in observations}
        if len(states) > 1:
            return 'PHASED_UNCERTAIN', 'conflicting phase evidence'
        return observations[0][0], '; '.join(sorted({x[1] for x in observations}))
    if any('|' in (x or '') for x in (a.gt, b.gt, a.pgt, b.pgt)):
        return 'PHASED_UNCERTAIN', 'no comparable validated phase block'
    return 'UNPHASED', 'no phase evidence'


def deduplicate(variants):
    """Collapse transcript/record duplicates; conflicting calls fail QC, never add alleles."""
    distinct = {}
    for v in sorted(variants, key=lambda x: (x.key, json.dumps(asdict(x), sort_keys=True))):
        if v.key not in distinct:
            distinct[v.key] = Variant(**asdict(v))
            continue
        prior = distinct[v.key]
        # Equal genomic alleles with different observations are not independent evidence.
        for name in ('gt', 'dp', 'gq', 'ab', 'filter', 'ps', 'pid', 'pgt', 'parents', 'copy_number', 'expected_copy_number', 'structural_qc', 'reviewed_acmg'):
            if getattr(prior, name) != getattr(v, name):
                prior.qc_conflict = True
        # Do not transfer transcript relevance between distinct consequences.
        if prior.consequence != v.consequence or prior.transcript_relevance != v.transcript_relevance:
            prior.transcript_relevance = 'UNKNOWN'
        prior.consequence = '&'.join(sorted(set(re.split('[&,]', prior.consequence + '&' + v.consequence))))
        unreliable = next((x for x in (prior.population_status, v.population_status)
                           if x in ('LOW_CALLABILITY', 'UNRELIABLE_LOCUS')), None)
        if v.af is not None:
            prior.af = max(x for x in (prior.af, v.af) if x is not None)
            prior.population_status = 'OBSERVED'
        if unreliable:
            prior.population_status = unreliable
    return list(distinct.values())


def evaluate_gene(variants, c=ScoreConfig(), phase_evidence=()):
    variants = deduplicate(variants)
    if not variants:
        return [], []
    if len({v.gene for v in variants}) != 1:
        raise ValueError('evaluate_gene requires one gene')
    evidence = {v.key: variant_evidence(v, c) for v in variants}
    hypotheses = []
    gene = variants[0].gene

    def emit(model, members=(), compatible=True, bonus=0, phase='NOT_APPLICABLE',
             phase_source='', second='NOT_APPLICABLE', note='', conflict=False):
        ev = [evidence[v.key] for v in members]
        viable = compatible and bool(ev) and all(e['variant_strength'] != 'UNSUPPORTED' for e in ev)
        if not viable or conflict:
            bonus = 0
        base = max((e['variant_score'] for e in ev), default=0)
        score = round(base + bonus, 3) if viable and not conflict else 0
        status = ('CONFLICTING_EVIDENCE' if conflict else 'INSUFFICIENT_EVIDENCE' if not viable
                  else 'HIGH_PRIORITY' if score >= c.high_priority
                  else 'MODERATE_PRIORITY' if score >= c.moderate_priority else 'EXPLORATORY')
        identity = json.dumps([gene, model, sorted(v.key for v in members)], separators=(',', ':'))
        qc_states = [e.get('qc_support', 'REVIEW') for e in ev]
        technical_tier = ('EXCLUDED' if any(q == 'FAIL' for q in qc_states)
                          else 'TECHNICAL_REVIEW' if any(q in ('REVIEW', 'UNKNOWN') for q in qc_states)
                          else 'PRIMARY')
        row = {'gene': gene, 'symbol': variants[0].symbol,
               'hypothesis_id': sha256(identity.encode()).hexdigest()[:24], 'inheritance_model': model,
               'variant_1': members[0].key if members else '', 'variant_1_score': ev[0]['variant_score'] if ev else '',
               'variant_2': members[1].key if len(members) > 1 else '', 'variant_2_score': ev[1]['variant_score'] if len(ev) > 1 else '',
               'zygosity_pattern': ';'.join(genotype(v.gt) for v in members),
               'phase_status': phase, 'phase_source': phase_source, 'second_allele_strength': second,
               'population_support': ';'.join(e['population_status'] + ':' + e['population_context'] for e in ev),
               'qc_support': ';'.join(e['qc_support'] for e in ev),
               'gene_disease_inheritance_support': 'UNKNOWN_G0_NOT_QUERIED',
               'inheritance_score': bonus, 'case_score': score, 'case_status': status,
               'technical_rank_tier': technical_tier,
               'technical_rank_eligible': technical_tier == 'PRIMARY',
               'model_applicability': 'COMPATIBLE' if compatible else 'NOT_SUPPORTED_BY_OBSERVATIONS',
               'acmg_status_variant_1': ev[0]['acmg_status'] if ev else 'NOT_APPLICABLE',
               'acmg_status_variant_2': ev[1]['acmg_status'] if len(ev) > 1 else 'NOT_APPLICABLE',
               # PM3 is disease-specific and cannot be assigned in phenotype-blind G0.
               'PM3': 'NOT_ASSIGNED',
               'compound_het_confirmed': bool(viable and not conflict and phase == 'IN_TRANS_CONFIRMED'),
               'rationale': note, 'ranking_status': 'PROVISIONAL_RESEARCH_HEURISTIC'}
        hypotheses.append(row)

    for v in variants:
        e = evidence[v.key]
        autosomal = v.chrom in {str(i) for i in range(1, 23)}
        gt = genotype(v.gt)
        if autosomal and gt == 'HOMOZYGOUS_ALT':
            emit('AR_HOMOZYGOUS', (v,), bonus=0,
                 note='Homozygous alternate genotype observed; no pathogenicity bonus in G0 because disease mechanism is unreviewed')
        if autosomal and gt in ('HETEROZYGOUS', 'HOMOZYGOUS_ALT'):
            emit('AD_SINGLE_ALLELE', (v,), note='Single allele compatible; dominant disease mechanism and segregation unreviewed')
        if v.parents:
            de_novo = parent_origin(v, c) == 'neither' and gt in ('HETEROZYGOUS', 'HEMIZYGOUS_OR_HAPLOID')
            emit('DE_NOVO_COMPATIBLE', (v,), compatible=de_novo,
                 bonus=c.de_novo if de_novo and e['qc_support'] == 'PASS' else 0,
                 note='Parental reference calls evaluated; parentage and mosaicism unreviewed; PS2/PM6 not assigned')
        if v.chrom == 'X':
            emit('X_LINKED', (v,), compatible=gt in ('HEMIZYGOUS_OR_HAPLOID', 'HETEROZYGOUS', 'HOMOZYGOUS_ALT'),
                 note='X-linked candidate; sex, PAR, ploidy and disease mechanism require review')
        if v.is_sv:
            emit('CNV_SV_COMPATIBLE', (v,), note='Structural/copy-number candidate; breakpoint, dosage and mechanism review required')
        emit('UNRESTRICTED_OTHER', (v,), note='Unrestricted candidate; no disease inheritance model assumed')
    # Only candidate alleles enter pair enumeration. Unsupported alleles still appear
    # in the variant table and, when paired with a candidate, receive zero bonus.
    for a, b in combinations(variants, 2):
        if a.chrom not in {str(i) for i in range(1, 23)} or a.chrom != b.chrom or a.build != b.build:
            continue
        ea, eb = evidence[a.key], evidence[b.key]
        # Pair construction is downstream of technical QC. Review-only alleles
        # remain visible individually, but do not create combinatorial recessive
        # hypotheses until both calls are technically primary/pass.
        if ea.get('qc_support') != 'PASS' or eb.get('qc_support') != 'PASS':
            continue
        if ea['variant_strength'] == eb['variant_strength'] == 'UNSUPPORTED':
            continue
        first, second = sorted((a, b), key=lambda v: (-evidence[v.key]['variant_score'], v.key))
        weakest = min((ea['variant_strength'], eb['variant_strength']), key={'UNSUPPORTED': 0, 'POSSIBLE': 1, 'STRONG': 2}.get)
        phase, source = phase_pair(a, b, c, phase_evidence)
        hets = genotype(a.gt) == genotype(b.gt) == 'HETEROZYGOUS'
        cis = phase == 'IN_CIS_CONFIRMED'
        # Overlapping sequence representations need normalization/haplotype review.
        overlap = (not a.is_sv and not b.is_sv and a.pos <= b.pos + len(b.ref) - 1
                   and b.pos <= a.pos + len(a.ref) - 1)
        conflict = cis or 'conflicting' in source
        supported = hets and not cis and weakest != 'UNSUPPORTED' and (not overlap or phase == 'IN_TRANS_CONFIRMED')
        # Merely observing two rare/credible alleles in one gene is not positive
        # pathogenicity evidence. Unresolved pairs remain candidates but receive
        # no score bonus. Confirmed-trans evidence contributes only when both
        # alleles have independent support beyond rarity.
        bonus = 0
        note = ('CANDIDATE BIALLELIC HYPOTHESIS — PHASE UNRESOLVED' if phase in ('UNPHASED', 'PHASED_UNCERTAIN')
                else 'Two distinct alleles confirmed in trans; pathogenicity and disease mechanism remain unreviewed' if phase == 'IN_TRANS_CONFIRMED'
                else 'In cis: not two recessive alleles' if cis else 'Not two diploid heterozygous alleles')
        if overlap:
            note += '; overlapping representations require haplotype review'
        emit('AR_POSSIBLE_COMPOUND_HET', (first, second), compatible=supported,
             bonus=bonus, phase=phase, phase_source=source,
             second='SECOND_ALLELE_' + weakest, note=note, conflict=conflict)
    for model in MODELS:
        if not any(h['inheritance_model'] == model for h in hypotheses):
            emit(model, compatible=False, note='No qualifying observed genotype for this model; not biologically excluded')
    return list(evidence.values()), sorted(hypotheses, key=lambda h: (-h['case_score'], h['hypothesis_id']))

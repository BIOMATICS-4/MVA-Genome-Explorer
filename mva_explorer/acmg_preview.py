"""Independent ACMG/AMP-style research evidence preview.

This module is deliberately read-only with respect to MVA ranking, QC, inheritance,
and BioGraph. It translates already-generated evidence plus explicit manual review
inputs into a transparent criterion map. It is not a clinical ACMG classification.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Literal

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

ManualState = Literal['yes', 'no', 'unknown']


class ManualReview(BaseModel):
    variant_key: str | None = None
    confirmed_in_trans: ManualState = 'unknown'
    de_novo: ManualState = 'unknown'
    functional_support: ManualState = 'unknown'
    segregation_support: ManualState = 'unknown'
    phenotype_specific: ManualState = 'unknown'
    lof_mechanism_confirmed: ManualState = 'unknown'


def _split(value):
    return [x.strip() for x in str(value or '').split(';')]


def _clinical_from_source(source: str) -> str:
    s = (source or '').upper()
    if 'CLINVAR_PATHOGENIC' in s and 'LIKELY' not in s:
        return 'Pathogenic'
    if 'CLINVAR_LIKELY_PATHOGENIC' in s:
        return 'Likely pathogenic'
    if 'CLINVAR_VUS' in s or 'UNCERTAIN' in s:
        return 'Uncertain significance'
    if 'CLINVAR_LIKELY_BENIGN' in s:
        return 'Likely benign'
    if 'CLINVAR_BENIGN' in s:
        return 'Benign'
    if 'CONFLICT' in s:
        return 'Conflicting'
    return 'Not available'


def _criterion(code, direction, strength, status, reason, source='pipeline', manual=False):
    return {
        'code': code, 'direction': direction, 'strength': strength,
        'status': status, 'reason': reason, 'source': source, 'manual': manual,
    }


def _path_counts(criteria):
    met = [c for c in criteria if c['status'] == 'MET']
    return {
        'PVS': sum(c['direction']=='pathogenic' and c['strength']=='very_strong' for c in met),
        'PS': sum(c['direction']=='pathogenic' and c['strength']=='strong' for c in met),
        'PM': sum(c['direction']=='pathogenic' and c['strength']=='moderate' for c in met),
        'PP': sum(c['direction']=='pathogenic' and c['strength']=='supporting' for c in met),
        'BA': sum(c['direction']=='benign' and c['strength']=='stand_alone' for c in met),
        'BS': sum(c['direction']=='benign' and c['strength']=='strong' for c in met),
        'BP': sum(c['direction']=='benign' and c['strength']=='supporting' for c in met),
    }


def combine_2015(criteria):
    """Mechanical preview of the 2015 ACMG/AMP combination table.

    Criterion applicability/strength still requires expert review. This function only
    combines criteria already marked MET; it does not decide whether a criterion is valid.
    """
    n = _path_counts(criteria)
    pathogenic_evidence = n['PVS'] + n['PS'] + n['PM'] + n['PP']
    benign_evidence = n['BA'] + n['BS'] + n['BP']
    if pathogenic_evidence and benign_evidence:
        label = 'Uncertain significance'
    elif n['BA'] >= 1 or n['BS'] >= 2:
        label = 'Benign'
    elif (n['BS'] >= 1 and n['BP'] >= 1) or n['BP'] >= 2:
        label = 'Likely benign'
    elif (
        (n['PVS'] >= 1 and (n['PS'] >= 1 or n['PM'] >= 2 or (n['PM'] >= 1 and n['PP'] >= 1) or n['PP'] >= 2))
        or n['PS'] >= 2
        or (n['PS'] >= 1 and (n['PM'] >= 3 or (n['PM'] >= 2 and n['PP'] >= 2) or (n['PM'] >= 1 and n['PP'] >= 4)))
    ):
        label = 'Pathogenic'
    elif (
        (n['PVS'] >= 1 and n['PM'] >= 1)
        or (n['PS'] >= 1 and 1 <= n['PM'] <= 2)
        or (n['PS'] >= 1 and n['PP'] >= 2)
        or n['PM'] >= 3
        or (n['PM'] >= 2 and n['PP'] >= 2)
        or (n['PM'] >= 1 and n['PP'] >= 4)
    ):
        label = 'Likely pathogenic'
    else:
        label = 'Uncertain significance'

    # Research-only directional index to make counterfactual manual edits visible.
    # It is NOT an ACMG clinical point system and never feeds MVA ranking.
    weights = {
        ('pathogenic','very_strong'): 8, ('pathogenic','strong'): 4,
        ('pathogenic','moderate'): 2, ('pathogenic','supporting'): 1,
        ('benign','stand_alone'): -8, ('benign','strong'): -4,
        ('benign','supporting'): -1,
    }
    index = sum(weights.get((c['direction'], c['strength']), 0)
                for c in criteria if c['status'] == 'MET')
    return {'classification': label + ' (preview)', 'counts': n,
            'research_evidence_index': index}


def evaluate_variant(gene_data, variant_key, manual: ManualReview):
    variants = gene_data['variants']
    selected = next((v for v in variants if v['key'] == variant_key), variants[0])
    other = next((v for v in variants if v['key'] != selected['key']), None)
    criteria = []

    effect = (selected.get('effect') or '').lower()
    lof = any(x in effect for x in ('stop_gained','frameshift','splice_acceptor','splice_donor','start_lost'))
    if lof:
        status = ('MET' if manual.lof_mechanism_confirmed == 'yes' else
                  'NOT_MET' if manual.lof_mechanism_confirmed == 'no' else 'REVIEW')
        criteria.append(_criterion(
            'PVS1','pathogenic','very_strong',status,
            'Predicted loss-of-function consequence detected. Full PVS1 requires confirmation that LoF is an established disease mechanism and transcript/NMD caveats are satisfied.',
            'VEP + manual disease-mechanism review', manual.lof_mechanism_confirmed != 'unknown'))

    freq = selected.get('max_freq_pct')
    try: freq_num = float(freq) if freq not in ('', None) else None
    except (TypeError, ValueError): freq_num = None
    if freq_num is not None:
        if freq_num > 5:
            criteria.append(_criterion('BA1','benign','stand_alone','REVIEW',
                f'Observed maximum population frequency is {freq_num:g}%. BA1 requires criterion-specific exception review.',
                selected.get('max_freq_source') or 'population resource'))
        elif freq_num <= 0.01:
            criteria.append(_criterion('PM2','pathogenic','supporting','REVIEW',
                f'Variant is very rare ({freq_num:g}%), but rarity/absence and disease-specific thresholds must be reviewed before assigning PM2.',
                selected.get('max_freq_source') or 'population resource'))
        else:
            criteria.append(_criterion('BS1','benign','strong','REVIEW',
                f'Observed population frequency is {freq_num:g}%; disease-specific maximum credible frequency is required before BS1 can be assigned.',
                selected.get('max_freq_source') or 'population resource'))

    psource = (selected.get('pathogenicity_source') or '').upper()
    if 'CODING_PREDICTOR' in psource or 'PREDICTOR' in psource:
        criteria.append(_criterion('PP3','pathogenic','supporting','REVIEW',
            'Computational predictor evidence is present, but independence/concordance of algorithms must be reviewed.',
            'coding predictor'))

    if manual.confirmed_in_trans == 'yes':
        if other and other.get('clinvar') in ('Pathogenic','Likely pathogenic'):
            criteria.append(_criterion('PM3','pathogenic','moderate','MET',
                'Manual review confirms the selected allele is in trans with another pathogenic/likely pathogenic allele in the recessive pair.',
                'manual phase review', True))
        else:
            criteria.append(_criterion('PM3','pathogenic','moderate','REVIEW',
                'Trans phase was manually confirmed, but the pathogenic status of the other allele is not established in the available evidence.',
                'manual phase review', True))
    elif manual.confirmed_in_trans == 'no':
        criteria.append(_criterion('PM3','pathogenic','moderate','NOT_MET',
            'Manual review indicates the pair is not confirmed in trans.', 'manual phase review', True))
    elif len(variants) > 1:
        criteria.append(_criterion('PM3','pathogenic','moderate','PENDING',
            'A candidate recessive pair exists, but phase is unresolved.', 'MVA inheritance layer'))

    if manual.de_novo == 'yes':
        criteria.append(_criterion('PM6','pathogenic','moderate','MET',
            'De novo occurrence was manually asserted; parentage-confirmed PS2 is not assumed.', 'manual review', True))
    elif manual.de_novo == 'no':
        criteria.append(_criterion('PM6','pathogenic','moderate','NOT_MET',
            'Manual review indicates the variant is not de novo.', 'manual review', True))

    if manual.functional_support == 'yes':
        criteria.append(_criterion('PS3','pathogenic','strong','MET',
            'User confirms a well-established functional study supports a damaging effect.', 'manual review', True))
    elif manual.functional_support == 'no':
        criteria.append(_criterion('PS3','pathogenic','strong','NOT_MET',
            'No qualifying damaging functional evidence asserted.', 'manual review', True))

    if manual.segregation_support == 'yes':
        criteria.append(_criterion('PP1','pathogenic','supporting','MET',
            'User confirms qualifying cosegregation evidence.', 'manual review', True))
    elif manual.segregation_support == 'no':
        criteria.append(_criterion('PP1','pathogenic','supporting','NOT_MET',
            'No qualifying cosegregation evidence asserted.', 'manual review', True))

    if manual.phenotype_specific == 'yes':
        criteria.append(_criterion('PP4','pathogenic','supporting','MET',
            'User confirms that the phenotype is sufficiently specific for PP4 review criteria.', 'manual review', True))
    elif manual.phenotype_specific == 'no':
        criteria.append(_criterion('PP4','pathogenic','supporting','NOT_MET',
            'Phenotype specificity was not considered sufficient for PP4.', 'manual review', True))
    else:
        criteria.append(_criterion('PP4','pathogenic','supporting','REVIEW',
            'MVA phenotype similarity is case-fit evidence only; PP4 is not automatically assigned.', 'HiPhive/MVA phenotype layer'))

    combo = combine_2015(criteria)
    return {
        'selected_variant': selected,
        'criteria': criteria,
        'preview': combo,
        'manual': manual.model_dump(),
        'boundary': 'Independent research preview only; original MVA ranking/QC/BioGraph are unchanged.',
    }


def acmg_router(root: Path):
    router = APIRouter()
    out = root / 'analysis/exomiser_compatibility/current'

    def read_tsv(name):
        p = out / name
        if not p.exists(): return []
        with p.open() as f: return list(csv.DictReader(f, delimiter='\t'))

    def gene_data(gene):
        gene = gene.upper().strip()
        sub = next((r for r in read_tsv('submission_top10_v13.tsv') if r.get('gene','').upper()==gene), None)
        if not sub:
            raise HTTPException(404, f'{gene} is not present in the current audited Top 10.')
        full = read_tsv('full_ranking.tsv')
        pair = next((r for r in full if r.get('id') == sub.get('submission_id')), None)
        rankrow = next((r for r in full if r.get('id') == sub.get('ranking_id')), None)
        evidence_row = pair or rankrow or {}
        keys = [x for x in _split(sub.get('submission_variant_keys') or sub.get('ranking_variant_keys')) if x]
        effects = _split(sub.get('submission_effects') or evidence_row.get('effects'))
        qcs = _split(sub.get('submission_qc_variant_statuses'))
        psources = _split(evidence_row.get('pathogenicity_sources'))
        hgvs = _split(evidence_row.get('hgvs_ids'))
        freqs = _split(evidence_row.get('max_freq_pct'))
        freqsrc = _split(evidence_row.get('max_freq_source'))
        abs_ = _split(evidence_row.get('allele_balance'))
        dps = _split(evidence_row.get('dp'))
        gqs = _split(evidence_row.get('gq'))
        variants=[]
        for i,key in enumerate(keys):
            get=lambda arr: arr[i] if i < len(arr) else ''
            src=get(psources)
            variants.append({
                'key':key,'effect':get(effects),'qc':get(qcs) or 'UNKNOWN',
                'pathogenicity_source':src,'clinvar':_clinical_from_source(src),
                'hgvs':get(hgvs),'max_freq_pct':get(freqs),'max_freq_source':get(freqsrc),
                'allele_balance':get(abs_),'dp':get(dps),'gq':get(gqs),
            })
        return {
            'gene': gene, 'rank': sub.get('gene_rank'),
            'mva_score': sub.get('ranking_combined_score'),
            'ranking_moi': sub.get('ranking_moi'), 'submission_moi': sub.get('submission_moi'),
            'phase': sub.get('submission_phase_status') or evidence_row.get('phase_status') or 'UNKNOWN',
            'representation_status': sub.get('representation_status'),
            'ranking_score_unchanged': sub.get('ranking_score_unchanged'),
            'phenotype_score': evidence_row.get('phenotype_score') or rankrow.get('phenotype_score','') if rankrow else evidence_row.get('phenotype_score',''),
            'variants': variants,
        }

    @router.get('/api/acmg/{gene}')
    def get_map(gene: str):
        data=gene_data(gene)
        manual=ManualReview(variant_key=data['variants'][0]['key'] if data['variants'] else None)
        return {'gene_data':data, **evaluate_variant(data, manual.variant_key, manual)}

    @router.post('/api/acmg/{gene}/preview')
    def recompute(gene: str, body: ManualReview):
        data=gene_data(gene)
        key=body.variant_key or (data['variants'][0]['key'] if data['variants'] else '')
        if not key or key not in {v['key'] for v in data['variants']}:
            raise HTTPException(400, 'Choose a variant belonging to this candidate.')
        return {'gene_data':data, **evaluate_variant(data, key, body)}

    return router

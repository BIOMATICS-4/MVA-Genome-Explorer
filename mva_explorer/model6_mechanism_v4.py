"""Isolated Model 6 V4 mechanism/QC/AF ranking sensitivity.

General method only: no disease gene, phenotype case ID, or known answer names.
This is NOT ACMG variant classification, LoF disease mechanism confirmation,
or calibrated population pathogenicity estimation. V3 and production remain.
"""
from __future__ import annotations
from collections import defaultdict
from mva_explorer.model6_evidence_audited_v2 import (
    finite_score,evidence_priorities,exploratory_gene_priorities)
from mva_explorer.model6_sexaware_v3 import sexaware_evidence_priorities

PREDICTED_LOF=frozenset({
 'frameshift_variant','stop_gained','splice_acceptor_variant','splice_donor_variant'
})
UNPROVEN_PROTEIN=frozenset({
 'stop_lost','start_lost','inframe_insertion','inframe_deletion','protein_altering_variant'
})
MISSENSE=frozenset({'missense_variant'})
QC_REVIEW=frozenset({'QC_REVIEW'})

def evidence_class(effects):
    """Prioritize *predicted* disruption, not ACMG pathogenicity."""
    kinds=[{p.strip() for p in s.replace(',','&').split('&') if p.strip()}
           for s in str(effects or '').split(';') if s.strip()]
    if not kinds:return 5,'UNCLASSIFIED'
    values=[]
    for terms in kinds:
        if terms & PREDICTED_LOF:values.append(0)
        elif terms & UNPROVEN_PROTEIN:values.append(2)
        elif terms & MISSENSE:values.append(3)
        else:values.append(4)
    if all(v==0 for v in values):return 0,'ALL_ALLELES_PREDICTED_LOF'
    if 0 in values:return 1,'MIXED_LOF_AND_OTHER_UNVERIFIED'
    if all(v==2 for v in values):return 2,'INFRAME_OR_OTHER_PROTEIN_ALTERING'
    if all(v==3 for v in values):return 3,'MISSENSE_ONLY_UNCALIBRATED'
    return 4,'UNCLASSIFIED_OR_OTHER_MIXED_EFFECT'

def rerank_hypothesis_records(rows,*,exploratory=False):
    """Take original priority lane as a GATE, then impact class before HPO.
    
    Never infer population rarity from missing AF. AF status is a recorded
    uncertainty; frequency filtering lives in the original compatibility
    preprocessor when a measured frequency is available.
    """
    output=[]
    for record in rows:
        r=dict(record)
        severity,label=evidence_class(r.get('effects'))
        # Regulatory effects remain distinct, and no severity hack may
        # upgrade an unverified genotype/QC or gene–disease relation.
        r['mechanistic_severity_order']=severity
        r['mechanistic_severity_description']=label
        r['AF_measured_for_ranking']=str(r.get('AF_status',''))=='REPORTED_FREQUENCY'
        r['mechanism_warning']='PREDICTION_NOT_CONFIRMED_DISEASE_LOF'
        r['ranking_release_status']='EXPERIMENTAL_UNCALIBRATED'
        output.append(r)

    def key(r):
        base_lane=int(r['lane_number'])
        if exploratory:
            # Honor exploratory mechanism classes (with unresolved phase
            # explicitly labelled, not misrepresented as confirmed).
            effect_lane=0 if 'CODING_OR_ESSENTIAL_SPLICE' in r['effect_class'] else 1 if 'REGULATORY' in r['effect_class'] else 2
            return (
                effect_lane, r['mechanistic_severity_order'],
                r.get('QC_status') in QC_REVIEW,
                not bool(r.get('established_gencc')),
                -finite_score(r.get('gene_phenotype')),
                -finite_score(r.get('variant_score_lower_bound')),
                -finite_score(r.get('disease_HPO_support')),
                r['gene'],r['hypothesis'],r['disease_id'])
        return (
            base_lane, r['mechanistic_severity_order'],
            -finite_score(r.get('gene_phenotype')),
            -finite_score(r.get('variant_score_lower_bound')),
            -finite_score(r.get('disease_HPO_support')),
            r['gene'],r['hypothesis'],r['disease_id'])
    output.sort(key=key)
    seen=set();genes=[]
    for row in output:
        if row['gene'] in seen:continue
        seen.add(row['gene']);r=dict(row)
        r['rank']=len(genes)+1
        if exploratory and r.get('phase') not in {'NOT_APPLICABLE','HOMOZYGOUS_ALT','TRANS'}:
            r['warnings']=str(r.get('warnings') or '')+';UNRESOLVED_GENOTYPE_NOT_CLINICALLY_CONFIRMED'
        genes.append(r)
    return genes,output

def evidence_priorities_v4(ranked,profiles,qc_by_key=None,*,sex_aware=False):
    if sex_aware:_,rows=sexaware_evidence_priorities(ranked,profiles,qc_by_key)
    else:_,rows=evidence_priorities(ranked,profiles,qc_by_key)
    return rerank_hypothesis_records(rows)

def exploratory_priorities_v4(ranked,profiles,qc_by_key=None):
    base=exploratory_gene_priorities(ranked,profiles,qc_by_key)
    return rerank_hypothesis_records(base,exploratory=True)

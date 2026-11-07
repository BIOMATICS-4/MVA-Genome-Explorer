"""Model 6 V4.1 sensitivity: established-QC/MOI gates plus EXISTING
compatibility logit, which already combines phenotype and variant consequence.

Avoid posthoc invented LoF/missense priority cutoffs. Missing AF remains an
optimistic conditional score and must be displayed as unmeasured. Exploratory
order remains V2 gene-HPO-first, with genotypes explicitly labelled.
Never clinically calibrated. No disease-gene identifiers in code.
"""
from mva_explorer.model6_evidence_audited_v2 import (
 evidence_priorities,exploratory_gene_priorities,finite_score
)
from mva_explorer.model6_sexaware_v3 import sexaware_evidence_priorities
from mva_explorer.model6_mechanism_v4 import evidence_class

def rank_existing_hypotheses(rows,*,exploratory=False):
    if exploratory:
        # No severe-effect hard hierarchy; hypothesis discovery should not
        # demote highly phenotypically coherent missense diagnoses.
        ranked=list(rows)
        ranked.sort(key=lambda r:(
          0 if 'CODING_OR_ESSENTIAL_SPLICE' in r.get('effect_class','') else 1 if 'REGULATORY' in r.get('effect_class','') else 2,
          -finite_score(r.get('gene_phenotype')),
          r.get('QC_status')=='QC_REVIEW',
          not (r.get('established_gencc') in (True,'True')),
          -finite_score(r.get('variant_score_lower_bound')),
          -finite_score(r.get('disease_HPO_support')),
          r['gene'],r['hypothesis'],r['disease_id']))
    else:
        ranked=list(rows)
        ranked.sort(key=lambda r:(
          int(r['lane_number']),
          -finite_score(r.get('legacy_combined')),
          -finite_score(r.get('gene_phenotype')),
          -finite_score(r.get('variant_score_lower_bound')),
          r['gene'],r['hypothesis'],r['disease_id']))
    output=[];seen=set()
    for candidate in ranked:
        if candidate['gene'] in seen:continue
        seen.add(candidate['gene'])
        r=dict(candidate)
        r['rank']=len(output)+1
        r['model_variant']='V4.1_CONDITIONAL_LEGACY_COMBINED'
        r['sort_lane']=str(r.get('lane_number'))
        r['effect_evidence_class']=evidence_class(r.get('effects'))[1]
        r['uncertainty']='UNCALIBRATED_PREDICTION'
        if str(r.get('AF_status'))=='UNKNOWN_NOT_MEASURED':
            r['uncertainty']+=';AF_MISSING_SCORE_IS_OPTIMISTIC_CONDITIONAL_UPPER_BOUND'
        if exploratory and str(r.get('phase')) not in {'NOT_APPLICABLE','HOMOZYGOUS_ALT','TRANS'}:
            r['uncertainty']+=';GENOTYPE_PHASE_NOT_CONFIRMED'
        output.append(r)
    return output,ranked

def strict_v41(ranked,profiles,qc_by_key=None,*,sex_aware=False):
    if sex_aware:_,hyp=sexaware_evidence_priorities(ranked,profiles,qc_by_key)
    else:_,hyp=evidence_priorities(ranked,profiles,qc_by_key)
    return rank_existing_hypotheses(hyp)

def exploratory_v41(ranked,profiles,qc_by_key=None):
    rows=exploratory_gene_priorities(ranked,profiles,qc_by_key)
    return rank_existing_hypotheses(rows,exploratory=True)

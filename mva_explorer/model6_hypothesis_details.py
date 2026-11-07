"""Read-only, gene-neutral AD/AR detail from the exact frozen Model 6 input.

A per-allele AD-style score is not evidence that a disease is dominant.
An unphased pair is not confirmed biallelic. Scores are not probabilities.
"""
from math import isclose
from mva_explorer.exomiser_compatibility import (
    _score_variant, _variant_key, combined_score, phase_relation,
)
PASS = {'QC_PASS','QC_STRONG'}
HET = {'0/1','1/0','0|1','1|0'}
HOM = {'1/1','1|1'}


def detail(gene, strict, raw, variant_rows, qc, metadata, profile):
    if gene != strict['gene']:
        raise ValueError('Gene mismatch')
    pheno = float(strict.get('phenotype_score') or strict.get('ranked_phenotype_score'))
    keys = set((strict.get('variant_keys') or strict.get('ranked_variant_keys')).split(';'))
    moi = strict.get('moi') or strict.get('ranked_moi') or strict.get('best_moi')
    chosen_score = float(strict.get('combined_score') or strict.get('ranked_score'))
    seen = {}
    for row in variant_rows:
        if row.get('gene') != gene:
            continue
        key = _variant_key(row)
        if qc.get(key) not in PASS:
            continue
        v = _score_variant(row,metadata,{gene:pheno},set(),profile)
        if v and v['copies'] and (key not in seen or v['variant_score']>seen[key]['variant_score']):
            seen[key] = v
    if not keys or not keys.issubset(seen):
        raise ValueError('Ranked genotype allele missing from QC-passing source evidence')
    ad=[]
    for key,v in seen.items():
        if v['gt'] not in HET:
            continue
        af=v['max_freq_pct']
        if not v['whitelist'] and af is not None and af > float(profile['inheritance_max_percent']['AD']):
            continue
        ad.append(dict(key=key,gt=v['gt'],score=combined_score(pheno,v['variant_score'],profile),
                       allele_score=v['variant_score'],qc=qc[key],
                       af_status=v['population_evidence_status'],selected=moi=='AD' and keys=={key},
                       evidence='AD_STYLE_ONLY_NOT_DISEASE_INHERITANCE'))
    ad.sort(key=lambda x:(-x['score'],x['key']))
    ar=[]
    for r in raw:
        if r.get('gene')!=gene or r.get('moi')!='AR':
            continue
        pair=[k for k in str(r.get('variant_keys','')).split(';') if k]
        if not pair or not set(pair).issubset(seen):
            continue
        if len(pair)==1 and seen[pair[0]]['gt'] in HOM:
            phase='HOMOZYGOUS_ALT'
            strength=seen[pair[0]]['variant_score']
        elif len(pair)==2 and len(set(pair))==2 and all(seen[k]['gt'] in HET for k in pair):
            phase=phase_relation(seen[pair[0]],seen[pair[1]])
            if phase=='CIS': continue
            if r.get('phase_status')=='TRANS' and phase!='TRANS':
                raise ValueError('Unverified trans label')
            if r.get('phase_status') not in (phase,'UNRESOLVED'):
                raise ValueError('Phase provenance conflict')
            strength=sum(seen[k]['variant_score'] for k in pair)/2
        else:
            continue
        score=combined_score(pheno,strength,profile)
        if not isclose(score,float(r['combined_score']),abs_tol=2e-8):
            raise ValueError('AR score does not reproduce from the contributing variants')
        ar.append(dict(keys=pair,genotypes=[seen[k]['gt'] for k in pair],
                       score=score,allele_score=strength,phase=phase,
                       qc=[qc[k] for k in pair],
                       status='REVIEW_ONLY_UNRESOLVED_PHASE' if phase=='UNRESOLVED'
                       else 'GENOTYPE_COMPATIBLE_CAUSALITY_UNASSESSED'))
    ar.sort(key=lambda x:(-x['score'],x['keys']))
    selected=next((x for x in ad if x['selected']),None)
    if selected and not isclose(selected['score'],chosen_score,abs_tol=2e-8):
        raise ValueError('Selected individual AD score does not equal frozen ranked score')
    if moi=='AD' and not selected:
        v=seen[next(iter(keys))] if len(keys)==1 else None
        flag=('AD_STYLE_ON_HOMO_ALT_NEEDS_REVIEW' if v and v['gt'] in HOM
              else 'AD_STYLE_NOT_MATCHED_TO_ELIGIBLE_HET')
    else:
        flag='STYLE_SCORE_NOT_ESTABLISHED_DISEASE_MOI'
    return dict(gene=gene,ranked_moi=moi,ranked_score=chosen_score,
                ranked_keys=sorted(keys),ranked_flag=flag,
                individual_ad=ad,ar=ar,clinical_moi_proven=False)

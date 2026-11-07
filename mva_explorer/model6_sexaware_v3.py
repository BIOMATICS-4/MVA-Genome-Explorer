"""Isolated sex-aware extension for the Model 6 compatibility scorer (research use).

Accept explicitly provided sample sex, never infer it from known causal genes.
For GRCh38 X non-PAR male calls, normalize pseudo-diploid 1/1 representations
as a single alternate X allele for genotype modelling, preserving the
original genotypes and remaining CNV uncertainty. Reject ambiguous male
0/1 on X nonPAR from straightforward hemizygous classification.
Do not apply to PAR1/PAR2, female X or subjects without reported sex.
"""
from collections import defaultdict
from mva_explorer.exomiser_compatibility import rank_variant_rows, genotype_copies
from mva_explorer.model6_evidence_audited_v2 import (
    evidence_priorities,finite_score
)

X_PAR_GRCH38=((10001,2781479),(155701383,156030895))

def sex_ploidy_context(row, sex, reference_build="GRCh38"):
    chrom=str(row.get("chrom","")).upper().removeprefix("CHR")
    pos=int(row["pos"])
    if reference_build!="GRCh38":
        raise ValueError("GRCh38 needed for explicit X PAR intervals")
    if chrom!="X":return "AUTOSOMAL_OR_OTHER"
    if any(a<=pos<=b for a,b in X_PAR_GRCH38):return "X_PSEUDOAUTOSOMAL"
    sex=str(sex or "").upper()
    if sex=="MALE":return "MALE_X_NONPAR"
    if sex=="FEMALE":return "FEMALE_X_NONPAR"
    return "UNKNOWN_SEX_X_NONPAR"

def sexaware_rank_variant_rows(rows,metadata,phenotype,whitelist,profile,
                               *,sample_sex,ar_phase_policy="ALLOW_UNRESOLVED",
                               include_xlinked=True,**kwargs):
    if str(sample_sex).upper() not in {"MALE","FEMALE","UNKNOWN"}:
        raise ValueError("sample sex must be MALE FEMALE or UNKNOWN")
    adjusted=[]
    source_gt={}
    blocked=[]
    for record in rows:
        original=dict(record)
        key=(str(record.get("build") or "GRCh38"),
             str(record["chrom"]).removeprefix("chr"),str(record["pos"]),
             str(record["ref"]),str(record["alt"]),record["gene"])
        source_gt[key]=str(record.get("gt") or "")
        context=sex_ploidy_context(record,sample_sex,key[0])
        if context=="MALE_X_NONPAR":
            raw=str(record.get("gt") or "")
            if raw in ("1/1","1|1","1"):
                original["gt"]="1"
            elif raw in ("0/0","0|0","0"):
                original["gt"]="0"
            else:
                blocked.append({"key":key,"GT":raw,
                    "reason":"male_X_nonPAR_diploid_heterozygote_ambiguous_ploidy_or_mosaicism"})
                continue
        adjusted.append(original)
    ranked=rank_variant_rows(adjusted,metadata,phenotype,whitelist,profile,
                             ar_phase_policy=ar_phase_policy,include_xlinked=include_xlinked,**kwargs)
    result=[]
    for candidate in ranked:
        alleles=candidate["variants"]
        contexts=[sex_ploidy_context(dict(chrom=v["key"].split(":")[1],pos=v["key"].split(":")[2]),sample_sex)
                  for v in alleles]
        if contexts and all(c=="MALE_X_NONPAR" for c in contexts):
            if candidate["moi"] in ("AD","AR"):
                # These are generic inheritance candidates from the compatibility
                # core and MUST NOT masquerade as gene-specific X disease validity.
                continue
            if candidate["moi"] in ("XL","XLR"):
                candidate["phase_status"]="XL_MALE_NONPAR_ALT_GT_COMPATIBLE_CN_UNVERIFIED"
                candidate["sex_ploidy_context"]="MALE_X_NONPAR_GT_NORMALIZED"
        for v in alleles:
            parts=v["key"].split(":")
            original=source_gt.get((parts[0],parts[1],parts[2],parts[3],parts[4],v["gene"]))
            v["source_GT"]=original
            v["genotype_normalized_for_sex"]=bool(original and original!=v["gt"])
        result.append(candidate)
    return result,blocked

def sexaware_evidence_priorities(ranked, profiles, qc_by_key=None):
    _,all_rows=evidence_priorities(ranked,profiles,qc_by_key)
    for r in all_rows:
        if (r["moi"] in ("XL","XLR")
                and r["phase"]=="XL_MALE_NONPAR_ALT_GT_COMPATIBLE_CN_UNVERIFIED"
                and r["established_gencc"]
                and r["QC_status"]!="QC_REVIEW"):
            kind=r["effect_class"]
            if kind=="PREDICTED_CODING_OR_ESSENTIAL_SPLICE":
                r["lane_number"]=0
                r["lane_label"]="XL_MALE_NONPAR_CODING_GT_COMPATIBLE_COPY_NUMBER_NOT_CONFIRMED"
            elif kind=="PREDICTED_REGULATORY_OR_NONCODING":
                r["lane_number"]=2
                r["lane_label"]="XL_MALE_NONPAR_REGULATORY_EFFECT_UNCONFIRMED"
            else:
                r["lane_number"]=3
                r["lane_label"]="XL_MALE_NONPAR_EFFECT_UNCLASSIFIED"
            r["warnings"]+=";SEX_METADATA_MALE_NONPAR_DIPLOID_GT_AS_HEMIZYGOUS;COPY_NUMBER_NOT_ASSAYED"
    all_rows.sort(key=lambda r:(
        r["lane_number"],-r["gene_phenotype"],-r["variant_score_lower_bound"],
        -finite_score(r["disease_HPO_support"]),-r["variant_score_upper_bound"],
        r["gene"],r["hypothesis"],r["disease_id"]))
    result=[];seen=set()
    for row in all_rows:
        if row["gene"] in seen:continue
        seen.add(row["gene"])
        copy=dict(row);copy["rank"]=len(result)+1;result.append(copy)
    return result,all_rows

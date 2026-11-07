"""Model 6 evidence-led experimental ranking, V2 (2026-10-18).

A hypothesis-led PRIORITY ORDER, not a calibrated causal probability.
Immutable original v1.2 scores are carried as evidence, not refitted.
Variants evaluated identically regardless of known causal gene or case ID.
"""
from __future__ import annotations
from collections import defaultdict
import math

CODING_OR_SPLICE = frozenset((
 "frameshift_variant","stop_gained","stop_lost","start_lost",
 "splice_acceptor_variant","splice_donor_variant",
 "missense_variant","protein_altering_variant","inframe_deletion",
 "inframe_insertion"
))
REGULATORY = frozenset((
 "5_prime_UTR_variant","3_prime_UTR_variant","upstream_gene_variant",
 "downstream_gene_variant","non_coding_transcript_exon_variant",
 "non_coding_transcript_variant","regulatory_region_variant",
 "intron_variant"
))
CONFIRMED_PHASE = {"NOT_APPLICABLE","HOMOZYGOUS_ALT","TRANS"}
MIN_DP=10
MIN_GQ=20

def effect_class(candidate):
    effects = set()
    for v in candidate["variants"]:
        effect=v.get("effect","")
        effects.update(part for part in effect.replace(",", "&").split("&") if part)
    if effects & CODING_OR_SPLICE:
        return "PREDICTED_CODING_OR_ESSENTIAL_SPLICE"
    if effects & REGULATORY:
        return "PREDICTED_REGULATORY_OR_NONCODING"
    return "UNCLASSIFIED_EFFECT"

def quality_status(candidate, external_qc=None):
    """Explicit measured failures only, missing read evidence stays unknown."""
    states=[]
    for v in candidate["variants"]:
        key=v["key"]
        known=(external_qc or {}).get(key)
        if known is not None:
            states.append(known)
            continue
        dp=v.get("dp");gq=v.get("gq")
        if dp is None or gq is None:
            states.append("QC_UNKNOWN")
        elif dp < MIN_DP or gq < MIN_GQ:
            states.append("QC_REVIEW_MEASURED_LOW")
        else:
            states.append("QC_BASIC_THRESHOLDS_ONLY")
    if any(x in ("QC_FAIL","QC_REVIEW","QC_REVIEW_MEASURED_LOW") for x in states):
        return "QC_REVIEW"
    if all(x in ("QC_PASS","QC_STRONG") for x in states):
        return "QC_INDEPENDENTLY_SUPPORTED"
    return "QC_PARTIAL_OR_UNKNOWN"

def genotype_review(candidate):
    moi=candidate["moi"]
    phase=candidate["phase_status"]
    if moi=="AR":
        return phase not in ("HOMOZYGOUS_ALT","TRANS")
    if moi in ("XL","XLR"):
        # XL requires sex, ploidy, X inactivation/segregation; not supplied here.
        return True
    if moi=="AD":
        return phase!="NOT_APPLICABLE"
    return True

def ranking_lane(candidate, profile, qc_status):
    # Source and diagnosis intentionally ABSENT from decision logic.
    if genotype_review(candidate):
        return 4,"MOI_PHASE_OR_SEX_PLOIDY_UNVERIFIED"
    if qc_status=="QC_REVIEW":
        return 5,"KNOWN_QC_CONCERN"
    if not profile.get("established_gene_disease_moi"):
        return 3,"GENE_DISEASE_VALIDITY_REVIEW"
    kind=effect_class(candidate)
    if kind=="PREDICTED_CODING_OR_ESSENTIAL_SPLICE":
        return 0,"ESTABLISHED_MOI_PREDICTED_CODING_OR_SPLICE"
    if kind=="PREDICTED_REGULATORY_OR_NONCODING":
        return 2,"REGULATORY_MECHANISM_REVIEW_NOT_CAUSAL_PROOF"
    return 3,"UNCLASSIFIED_VARIANT_MECHANISM"

def finite_score(value):
    try:
        value=float(value)
        return value if math.isfinite(value) else 0.0
    except (TypeError,ValueError):
        return 0.0

def evidence_priorities(ranked, profiles, qc_by_key=None):
    """Return ALL eligible hypotheses in a gene-distinct headline and review lanes.

    Sorting: biological evidence class > gene-level phenotype > conservative
    variant lower-bound if available > disease phenotype (last-resort tiebreak).
    Disease HPO is deliberately NOT a second vote (same correlated source).
    """
    by_hyp={(r["gene"],r["moi"]):r for r in ranked}
    gene_mappings=defaultdict(set)
    for subject in ranked:
        for allele in subject["variants"]:gene_mappings[allele["key"]].add(subject["gene"])
    hypotheses=[]
    for profile in profiles:
        key=(profile["gene"],profile["moi"])
        candidate=by_hyp.get(key)
        if candidate is None:
            raise ValueError("Disease–genotype profile missing exact gene/MOI hypothesis")
        allele_keys=";".join(v["key"] for v in candidate["variants"])
        if profile["genotype_variant_keys"] != allele_keys:
            raise ValueError("Mismatched genotype alleles in disease profile")
        qc=quality_status(candidate,qc_by_key)
        lane,label=ranking_lane(candidate,profile,qc)
        lower=min((finite_score(v.get("variant_score_lower_bound")) for v in candidate["variants"]),default=0)
        variant_value=finite_score(candidate["variant_score"])
        missing_AF=any(v.get("max_freq_pct") is None for v in candidate["variants"])
        text_name=str(profile.get("disease_name") or "").lower()
        moi=candidate["moi"]
        conflicting_label=(moi=="AD" and "recessive" in text_name) or (moi=="AR" and "dominant" in text_name)
        multiple_gene_mappings=any(len(gene_mappings[v["key"]])>1 for v in candidate["variants"])
        row={
          "gene":candidate["gene"],"hypothesis":candidate["id"],
          "moi":candidate["moi"],"disease":profile["disease_name"],
          "disease_id":profile["disease_id"],
          "gene_phenotype":finite_score(candidate["phenotype_score"]),
          "disease_HPO_support":profile["hpo_disease_profile_score"],
          "legacy_combined":finite_score(candidate["combined_score"]),
          "variant_score_upper_bound":variant_value,
          "variant_score_lower_bound":lower,
          "effect_class":effect_class(candidate),
          "effects":";".join(v.get("effect","") for v in candidate["variants"]),
          "variant_keys":allele_keys,"genotype":";".join(str(v.get("gt","")) for v in candidate["variants"]),
          "phase":candidate["phase_status"],"QC_status":qc,
          "AF_status":"UNKNOWN_NOT_MEASURED" if missing_AF else "REPORTED_FREQUENCY",
          "variant_provenance":";".join(v.get("pathogenicity_source","") for v in candidate["variants"]),
          "lane_number":lane,"lane_label":label,
          "established_gencc":profile["established_gene_disease_moi"],
          "gencc_validity":profile["gencc_validity_classes"],
          "source_variant_count":len(candidate["variants"]),
          "warnings":";".join([
            *(['POPULATION_AF_UNMEASURED_UPPER_BOUND_NOT_RARITY'] if missing_AF else []),
            *(['QC_NOT_VERIFIED_FROM_PATIENT_READS'] if qc=="QC_PARTIAL_OR_UNKNOWN" else []),
            *(['REGULATORY_PREDICTION_NOT_VALIDATED_FUNCTIONAL_EFFECT'] if lane==2 else []),
            *(['SEX_PLOIDY_OR_PHASE_NOT_RESOLVED'] if lane==4 else []),
            *(['MULTIPLE_GENE_MAPPINGS_SAME_ALLELE'] if multiple_gene_mappings else []),
            *(['DISEASE_LABEL_MOI_CONFLICT_CHECK_GENCC'] if conflicting_label else []),
          ])
        }
        hypotheses.append(row)
    # Determine a genotype-matched disease hypothesis for each gene.
    hypotheses.sort(key=lambda r:(
        r["lane_number"],-r["gene_phenotype"],
        -r["variant_score_lower_bound"],
        -(finite_score(r["disease_HPO_support"])),
        -r["variant_score_upper_bound"],
        r["gene"],r["hypothesis"],r["disease_id"]))
    seen=set();genes=[]
    for r in hypotheses:
        if r["gene"] in seen:continue
        seen.add(r["gene"])
        new=dict(r);new["rank"]=len(genes)+1;genes.append(new)
    return genes,hypotheses

def exploratory_gene_priorities(ranked, profiles, qc_by_key=None):
    """Separate unconfirmed *hypothesis-generation* order.

    Phase-unresolved AR and X-linked hypotheses can be SHOWN here, not
    upgraded to confirmed disease genotypes. Within a gene, prefer the
    documented disease/MOI hypothesis over a convenient high-scoring AD-style
    single-allele hypothesis if that gene's disease is recessive.
    Across genes: predicted coding/splice vs regulatory-mechanism REVIEW,
    then gene-wide phenotype match; no target-gene knowledge used.
    """
    _,all_hyp=evidence_priorities(ranked,profiles,qc_by_key)
    by_gene=defaultdict(list)
    for row in all_hyp:by_gene[row["gene"]].append(row)

    def mechanism_tier(row):
        if row["effect_class"]=="PREDICTED_CODING_OR_ESSENTIAL_SPLICE":
            return 0
        if row["effect_class"]=="PREDICTED_REGULATORY_OR_NONCODING":
            return 1
        return 2
    selected=[]
    for gene,options in by_gene.items():
        # Evidence relation belongs to the EXACT gene/MOI hypothesis.
        chosen=min(options,key=lambda row:(
            not row["established_gencc"],
            mechanism_tier(row),
            row["QC_status"]=="QC_REVIEW",
            row["phase"] not in CONFIRMED_PHASE,
            -(finite_score(row["disease_HPO_support"])),
            -row["variant_score_lower_bound"],
            row["hypothesis"],row["disease_id"]))
        selected.append(dict(chosen))
    selected.sort(key=lambda row:(
        mechanism_tier(row),
        -row["gene_phenotype"],
        row["QC_status"]=="QC_REVIEW",
        not row["established_gencc"],
        -row["variant_score_lower_bound"],
        -finite_score(row["disease_HPO_support"]),
        row["gene"],row["hypothesis"]))
    for i,row in enumerate(selected,1):
        row["rank"]=i
        row["ranking_mode"]="EXPLORATORY_NO_CONFIRMED_GENOTYPE_CLAIM"
        if row["phase"] not in CONFIRMED_PHASE:
            row["warnings"]+=";GENOTYPE_NOT_CONFIRMED_FOR_DISEASE"
    return selected

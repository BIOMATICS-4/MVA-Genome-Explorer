"""Isolated Exomiser 13.2 compatibility scoring.

This module intentionally does not alter the MVA standards-first G0/G1/G2 path.
It reproduces generic Exomiser-style filtering and ranking for compatibility/QC.
"""
from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROFILE = ROOT / "config/exomiser_compatibility.json"


def number(value):
    try:
        if value in (None, "", "."):
            return None
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def load_profile(path=DEFAULT_PROFILE):
    profile = json.loads(Path(path).read_text())
    if profile.get("profile") != "exomiser-13.2-compatibility":
        raise ValueError("Not an Exomiser compatibility profile")
    return profile


def frequency_score(max_frequency_percent, profile):
    value = number(max_frequency_percent)
    if value is None or value <= 0:
        return 1.0
    if value > float(profile["frequency_filter_max_percent"]):
        return 0.0
    cfg = profile["frequency_score"]
    return max(0.0, float(cfg["intercept"]) -
               float(cfg["exponential_coefficient"]) * math.exp(value))


def combined_score(phenotype_score, variant_score, profile):
    cfg = profile["combined_score_logit"]
    z = (float(cfg["intercept"]) +
         float(cfg["phenotype_coefficient"]) * float(phenotype_score) +
         float(cfg["variant_coefficient"]) * float(variant_score))
    if z >= 50:
        return 1.0
    if z <= -50:
        return 0.0
    return 1.0 / (1.0 + math.exp(-z))


def effect_fallback_score(effect, profile):
    scores = profile["pathogenicity_fallback_scores"]
    terms = [x for x in re.split(r"[&,]", effect or "") if x]
    values = [float(scores[x]) for x in terms if x in scores]
    return max(values) if values else None


def is_regulatory_effect(effect, profile):
    terms = {x for x in re.split(r"[&,]", effect or "") if x}
    return bool(terms & set(profile.get("regulatory_effects", [])))


def regulatory_gene_linkage(effect, distance, profile):
    """Return a gene-linkage weight/class for a regulatory locus.

    ReMM estimates regulatory-locus pathogenicity. This separate layer estimates
    how defensible the selected VEP transcript/gene assignment is, using only
    consequence class and VEP DISTANCE; it never uses gene identity.
    """
    cfg = profile.get("regulatory_gene_linkage", {})
    terms = {x for x in re.split(r"[&,]", effect or "") if x}

    overlap = cfg.get("transcript_overlap_scores", {})
    overlap_hits = [(term, number(overlap.get(term))) for term in terms if term in overlap]
    overlap_hits = [(term, score) for term, score in overlap_hits if score is not None]
    if overlap_hits:
        term, score = max(overlap_hits, key=lambda x: x[1])
        return float(score), f"TRANSCRIPT_OVERLAP:{term}"

    value = number(distance)
    if value is not None:
        value = abs(value)

    def distance_rule(rules, prefix):
        if value is None:
            score = float(cfg.get("missing_distance_score", 0.0))
            return score, f"{prefix}_DISTANCE_UNKNOWN"
        for rule in rules:
            maximum = rule.get("max")
            if maximum is None or value <= float(maximum):
                return float(rule["score"]), str(rule["class"])
        return 0.0, f"{prefix}_UNLINKED"

    if "upstream_gene_variant" in terms:
        return distance_rule(cfg.get("upstream_distance_bp", []), "UPSTREAM")
    if "downstream_gene_variant" in terms:
        return distance_rule(cfg.get("downstream_distance_bp", []), "DOWNSTREAM")
    if "regulatory_region_variant" in terms:
        score = float(cfg.get("unlinked_regulatory_region_score", 0.0))
        return score, "UNLINKED_REGULATORY_REGION"

    score = float(cfg.get("missing_distance_score", 0.0))
    return score, "REGULATORY_LINKAGE_UNRESOLVED"


def genotype_copies(gt):
    parts = re.split(r"[/|]", gt or "")
    if not parts or any(x not in ("0", "1") for x in parts):
        return 0
    return sum(x == "1" for x in parts)


def read_phenotype_scores(path):
    """Read only the gene and Exomiser phenotype columns, avoiding huge evidence fields."""
    scores = {}
    with Path(path).open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        gene_key = "GENE_SYMBOL" if "GENE_SYMBOL" in header else "gene"
        score_key = ("EXOMISER_GENE_PHENO_SCORE"
                     if "EXOMISER_GENE_PHENO_SCORE" in header else "hiphive_score")
        gi, si = header.index(gene_key), header.index(score_key)
        max_index = max(gi, si)
        for line in handle:
            parts = line.rstrip("\n").split("\t", max_index + 1)
            if len(parts) <= max_index:
                continue
            gene = parts[gi]
            score = number(parts[si])
            if gene and score is not None:
                scores[gene] = max(scores.get(gene, 0.0), score)
    return scores


def read_variant_metadata(path):
    out = {}
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["hgvs_id"]] = {
                "max_freq_pct": number(row.get("max_freq_pct")),
                "max_freq_source": row.get("max_freq_source", ""),
                "revel": number(row.get("revel")),
                "mvp": number(row.get("mvp")),
                "pred_max": number(row.get("pred_max")),
                # ReMM is already on a 0-1 pathogenicity-like scale and is
                # intended for regulatory/non-coding Mendelian variants.
                "remm": number(row.get("remm") or row.get("remm_score")),
                "regulatory_pred": number(
                    row.get("regulatory_pred") or row.get("noncoding_pred")
                ),
            }
    return out


def read_pathogenic_whitelist(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = set()
    with opener(path, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                rows.add((fields[0].removeprefix("chr"), fields[1], fields[2], fields[3]))
    return rows


def _variant_key(row):
    if row.get("variant_key"):
        return row["variant_key"]
    build = row.get("build") or "GRCh37"
    chrom = str(row["chrom"]).removeprefix("chr")
    return f"{build}:{chrom}:{row['pos']}:{row['ref']}:{row['alt']}"


def _score_variant(row, metadata, phenotype, whitelist, profile):
    effect = row.get("effect") or row.get("consequence", "")
    base = effect_fallback_score(effect, profile)
    hgvs = row.get("hgvs_id", "")
    meta = metadata.get(hgvs, {})
    freq = number(row.get("max_freq_pct"))
    if freq is None:
        freq = number(meta.get("max_freq_pct"))
    chrom = str(row["chrom"]).removeprefix("chr")
    white = (
        str(row.get("pathogenic_whitelist", "")).lower() in ("1", "true", "yes")
        or str(row.get("clinvar_pathogenic_whitelist", "")).lower() in ("1", "true", "yes")
        or (chrom, str(row["pos"]), row["ref"], row["alt"]) in whitelist
    )
    if (freq is not None and
            freq > float(profile["frequency_filter_max_percent"]) and not white):
        return None

    pred = number(row.get("pred_max"))
    if pred is None:
        pred = number(meta.get("pred_max"))

    regulatory_pred = number(row.get("remm") or row.get("regulatory_pred"))
    if regulatory_pred is None:
        regulatory_pred = number(meta.get("remm"))
    if regulatory_pred is None:
        regulatory_pred = number(meta.get("regulatory_pred"))

    linkage_score = None
    linkage_class = "NOT_APPLICABLE"
    regulatory_effective_score = None
    distance_to_transcript = number(
        row.get("distance_to_transcript") or row.get("vep_distance") or row.get("distance")
    )
    if is_regulatory_effect(effect, profile):
        linkage_score, linkage_class = regulatory_gene_linkage(
            effect, distance_to_transcript, profile
        )

    # Coding/splice consequences use the calibrated fallback/predictor route.
    # Regulatory/non-coding consequences require explicit non-coding evidence
    # (e.g. ReMM) or a strict Pathogenic ClinVar whitelist override. This keeps
    # generic upstream/UTR annotations from receiving an invented rescue score.
    if base is None:
        if white and profile.get("clinvar_pathogenic_whitelist_override", False):
            pathogenicity = 1.0
            pathogenicity_source = "CLINVAR_PATHOGENIC"
        elif is_regulatory_effect(effect, profile) and regulatory_pred is not None:
            regulatory_effective_score = regulatory_pred * float(linkage_score or 0.0)
            if regulatory_effective_score <= 0:
                return None
            pathogenicity = regulatory_effective_score
            pathogenicity_source = "REMM_X_GENE_LINKAGE"
        else:
            return None
    else:
        pathogenicity = max(base, pred if pred is not None else -1.0)
        pathogenicity_source = (
            "CODING_PREDICTOR" if pred is not None and pred > base
            else "CONSEQUENCE_FALLBACK"
        )
        if white and profile.get("clinvar_pathogenic_whitelist_override", False):
            pathogenicity = max(pathogenicity, 1.0)
            if pathogenicity >= 1.0:
                pathogenicity_source = "CLINVAR_PATHOGENIC"
    # Missing AF is not a measured zero. The legacy numeric compatibility
    # score uses the neutral multiplicative factor 1, which is an optimistic
    # upper bound when the observed population frequency is not available.
    # Export that uncertainty explicitly; never call it proof of rarity.
    if freq is None:
        af_status = "UNKNOWN_NOT_MEASURED"
        score_lower = 0.0  # unknown AF could exceed the population gate
        score_upper = pathogenicity
    else:
        af_status = "REPORTED_ZERO_COVERAGE_UNVERIFIED" if freq == 0 else "REPORTED_NONZERO"
        score_lower = score_upper = pathogenicity * frequency_score(freq, profile)
    return {
        "population_evidence_status": af_status,
        "variant_score_lower_bound": score_lower,
        "variant_score_upper_bound": score_upper,
        "score_bound_semantics": (
            "CONDITIONAL_UPPER_BOUND_AF_UNOBSERVED" if freq is None
            else "CONDITIONAL_SCORE_WITH_REPORTED_AF"
        ),
        "gene": row["gene"],
        "key": _variant_key(row),
        "hgvs_id": hgvs,
        "copies": genotype_copies(row.get("gt", "")),
        "gt": row.get("gt", ""),
        "phase_set": row.get("phase_set", ""),
        "phase_gt": row.get("phase_gt", ""),
        "dp": number(row.get("dp")),
        "gq": number(row.get("gq")),
        "allele_balance": number(row.get("allele_balance")),
        "effect": row.get("effect") or row.get("consequence", ""),
        "pathogenicity_score": pathogenicity,
        "pathogenicity_source": pathogenicity_source,
        "regulatory_predictor": regulatory_pred,
        "regulatory_linkage_score": linkage_score,
        "regulatory_linkage_class": linkage_class,
        "regulatory_effective_score": regulatory_effective_score,
        "distance_to_transcript": distance_to_transcript,
        "max_freq_pct": freq,
        "max_freq_source": row.get("max_freq_source") or meta.get("max_freq_source", ""),
        "whitelist": white,
        "phenotype_score": float(phenotype.get(row["gene"], 0.0)),
        "variant_score": pathogenicity * frequency_score(freq, profile),
    }


def phase_relation(left, right):
    """Return TRANS, CIS or UNRESOLVED for two heterozygous variants.

    PGT is only interpreted when both calls share the same non-empty PID phase
    block. Same phased alternate haplotype means cis; opposite haplotypes means
    trans. Missing/different phase blocks are unresolved, not assumed trans.
    """
    left_set = str(left.get("phase_set") or "")
    right_set = str(right.get("phase_set") or "")
    left_gt = str(left.get("phase_gt") or "")
    right_gt = str(right.get("phase_gt") or "")
    valid = {"0|1": 1, "1|0": 0}
    if not left_set or left_set != right_set:
        return "UNRESOLVED"
    if left_gt not in valid or right_gt not in valid:
        return "UNRESOLVED"
    return "CIS" if valid[left_gt] == valid[right_gt] else "TRANS"


def rank_variant_rows(rows, metadata, phenotype, whitelist, profile=None,
                      inheritance_mode="ANY", *, qc_status_by_key=None,
                      ar_phase_policy="ALLOW_UNRESOLVED", include_xlinked=False,
                      ar_pair_aggregation="MEAN"):
    """Return Exomiser-compatibility gene/MOI candidates from variant-gene rows.

    inheritance_mode is a case-level constraint, not a gene rule. Use ANY when
    inheritance is unknown; AD or AR when the case definition explicitly
    provides that inheritance model.
    """
    profile = profile or load_profile()
    if ar_phase_policy not in {"ALLOW_UNRESOLVED", "CONFIRMED_ONLY"}:
        raise ValueError("Unsupported AR phase policy")
    if ar_pair_aggregation not in {"MEAN", "BOTTLENECK_MIN"}:
        raise ValueError("Unsupported AR allele evidence aggregation")
    inheritance_mode = str(inheritance_mode or "ANY").upper()
    if inheritance_mode not in {"ANY", "AD", "AR"}:
        raise ValueError("inheritance_mode must be ANY, AD or AR")
    by_gene = defaultdict(list)
    for row in rows:
        if qc_status_by_key is not None and qc_status_by_key.get(_variant_key(row)) not in {"QC_PASS", "QC_STRONG"}:
            # Apply strict technical QC BEFORE genotype pairing and candidate
            # selection, not after choosing one best-scoring pair per gene.
            continue
        scored = _score_variant(row, metadata, phenotype, whitelist, profile)
        if scored is not None and scored["copies"] > 0:
            by_gene[scored["gene"]].append(scored)

    limits = profile["inheritance_max_percent"]
    candidates = []
    for gene, variants in by_gene.items():
        pheno = max(v["phenotype_score"] for v in variants)
        ad = [v for v in variants if v["copies"] >= 1 and
              (v["whitelist"] or v["max_freq_pct"] is None or
               v["max_freq_pct"] <= float(limits["AD"]))]
        if ad and inheritance_mode in {"ANY", "AD"}:
            best = max(ad, key=lambda v: (v["variant_score"], v["key"]))
            candidates.append(_candidate(gene, "AD", pheno, best["variant_score"],
                                         [best], profile))

        hom = [v for v in variants if v["copies"] >= 2 and
               (v["whitelist"] or v["max_freq_pct"] is None or
                v["max_freq_pct"] <= float(limits["AR_HOM_ALT"]))]
        ar_options = []
        if hom:
            best = max(hom, key=lambda v: (v["variant_score"], v["key"]))
            ar_options.append((best["variant_score"], [best], "HOMOZYGOUS_ALT"))

        het = {}
        for variant in variants:
            if variant["copies"] != 1:
                continue
            if (not variant["whitelist"] and variant["max_freq_pct"] is not None and
                    variant["max_freq_pct"] > float(limits["AR_COMP_HET"])):
                continue
            current = het.get(variant["key"])
            if current is None or variant["variant_score"] > current["variant_score"]:
                het[variant["key"]] = variant
        ordered = sorted(het.values(), key=lambda v: (-v["variant_score"], v["key"]))
        if len(ordered) >= 2:
            for i, left in enumerate(ordered[:-1]):
                for right in ordered[i + 1:]:
                    relation = phase_relation(left, right)
                    if relation == "CIS" or (
                        ar_phase_policy == "CONFIRMED_ONLY" and relation != "TRANS"
                    ):
                        continue
                    # The original Exomiser-compatible MEAN is preserved. A
                    # bottleneck sensitivity alternative avoids a strong allele
                    # compensating for a weak second allele. Neither metric is
                    # a calibrated probability of biallelic pathogenicity.
                    if ar_pair_aggregation == "BOTTLENECK_MIN":
                        value = min(left["variant_score"], right["variant_score"])
                    else:
                        value = (left["variant_score"] + right["variant_score"]) / 2.0
                    ar_options.append((value, [left, right], relation))
        if ar_options and inheritance_mode in {"ANY", "AR"}:
            value, contributing, phase_status = max(
                ar_options,
                key=lambda item: (item[0], 1 if item[2] == "TRANS" else 0)
            )
            candidates.append(_candidate(
                gene, "AR", pheno, value, contributing, profile,
                phase_status=phase_status
            ))

        # An X-linked gene is NOT automatically AD or AR. Explicitly include
        # an X-linked single-allele hypothesis when requested by a genotype-
        # aware analysis. Without sex/ploidy and copy-number context, neither
        # heterozygote nor hemizygote causality is asserted. No gene-name logic.
        if include_xlinked and inheritance_mode == "ANY":
            chrom_x=[v for v in variants if v["key"].split(":")[1].upper().removeprefix("CHR")=="X"]
            if chrom_x:
                best_x=max(chrom_x,key=lambda v:(v["variant_score"],v["key"]))
                xgt=str(best_x.get("gt") or "")
                if xgt=="1":
                    xstatus="XL_HAPLOID_ALT_CONSISTENT_PLOIDY_UNVERIFIED"
                elif xgt in ("0/1","1/0","0|1","1|0"):
                    xstatus="XL_HETEROZYGOUS_DISEASE_MECHANISM_UNASSESSED"
                elif xgt in ("1/1","1|1"):
                    xstatus="XL_HOMOZYGOUS_ALT_MECHANISM_UNASSESSED"
                else:
                    xstatus="X_LINKED_PLOIDY_AND_DIAGNOSIS_UNRESOLVED"
                candidates.append(_candidate(
                    gene,"XL",pheno,best_x["variant_score"],[best_x],profile,
                    phase_status=xstatus
                ))
                # X-linked recessive may be supported by one affected X copy
                # in a haploid context, or by two alternate copies. Ordinary
                # diploid heterozygotes are retained as review-only by the XL
                # hypothesis, not promoted as XLR causality.
                if xgt in ("1", "1/1", "1|1"):
                    candidates.append(_candidate(
                        gene,"XLR",pheno,best_x["variant_score"],[best_x],profile,
                        phase_status="XLR_HAPLOID_OR_BIALLELIC_ALT_PLOIDY_UNVERIFIED"
                    ))

    candidates.sort(key=lambda x: (-x["combined_score"], -x["phenotype_score"],
                                   -x["variant_score"], x["id"]))
    for rank, candidate in enumerate(candidates, 1):
        candidate["rank"] = rank
    return candidates


def _candidate(gene, moi, phenotype_score, variant_score, variants, profile,
               phase_status="NOT_APPLICABLE"):
    return {
        "id": f"{gene}_{moi}",
        "gene": gene,
        "moi": moi,
        "phenotype_score": phenotype_score,
        "variant_score": variant_score,
        "combined_score": combined_score(phenotype_score, variant_score, profile),
        "phase_status": phase_status,
        "frequency_evidence_status": ";".join(v.get("population_evidence_status", "") for v in variants),
        "variant_score_lower_bound": sum(v.get("variant_score_lower_bound", 0.0) for v in variants) / len(variants),
        "variant_score_upper_bound": sum(v.get("variant_score_upper_bound", v["variant_score"]) for v in variants) / len(variants),
        "allelic_bottleneck_support": min(v["variant_score"] for v in variants),
        "allelic_mean_support": sum(v["variant_score"] for v in variants) / len(variants),
        "variants": variants,
    }


def rank_compatibility_files(variant_gene_path, metadata_path, phenotype_path,
                             whitelist_path, profile_path=DEFAULT_PROFILE,
                             inheritance_mode="ANY"):
    profile = load_profile(profile_path)
    metadata = read_variant_metadata(metadata_path)
    phenotype = read_phenotype_scores(phenotype_path)
    whitelist = read_pathogenic_whitelist(whitelist_path)
    with Path(variant_gene_path).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rank_variant_rows(
        rows, metadata, phenotype, whitelist, profile,
        inheritance_mode=inheritance_mode
    )


def write_ranking(path, candidates):
    fields = ["rank", "id", "gene", "moi", "combined_score", "phenotype_score",
              "variant_score", "phase_status", "variant_keys", "hgvs_ids", "effects",
              "pathogenicity_sources", "regulatory_predictors",
              "regulatory_linkage_scores", "regulatory_linkage_classes",
              "regulatory_distances", "regulatory_effective_scores",
              "phase_sets", "phase_gts", "allele_balance", "dp", "gq",
              "max_freq_pct", "max_freq_source", "whitelist",
              "frequency_evidence_status", "variant_score_lower_bound",
              "variant_score_upper_bound"]
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for item in candidates:
            variants = item["variants"]
            writer.writerow({
                "rank": item["rank"], "id": item["id"], "gene": item["gene"],
                "moi": item["moi"], "combined_score": f"{item['combined_score']:.12f}",
                "phenotype_score": f"{item['phenotype_score']:.12f}",
                "variant_score": f"{item['variant_score']:.12f}",
                "phase_status": item.get("phase_status", "NOT_APPLICABLE"),
                "variant_keys": ";".join(v["key"] for v in variants),
                "hgvs_ids": ";".join(v["hgvs_id"] for v in variants),
                "effects": ";".join(v["effect"] for v in variants),
                "pathogenicity_sources": ";".join(
                    v.get("pathogenicity_source", "") for v in variants
                ),
                "regulatory_predictors": ";".join(
                    "" if v.get("regulatory_predictor") is None
                    else f"{v['regulatory_predictor']:.12g}" for v in variants
                ),
                "regulatory_linkage_scores": ";".join(
                    "" if v.get("regulatory_linkage_score") is None
                    else f"{v['regulatory_linkage_score']:.12g}" for v in variants
                ),
                "regulatory_linkage_classes": ";".join(
                    str(v.get("regulatory_linkage_class") or "") for v in variants
                ),
                "regulatory_distances": ";".join(
                    "" if v.get("distance_to_transcript") is None
                    else str(int(v["distance_to_transcript"])) for v in variants
                ),
                "regulatory_effective_scores": ";".join(
                    "" if v.get("regulatory_effective_score") is None
                    else f"{v['regulatory_effective_score']:.12g}" for v in variants
                ),
                "phase_sets": ";".join(str(v.get("phase_set") or "") for v in variants),
                "phase_gts": ";".join(str(v.get("phase_gt") or "") for v in variants),
                "allele_balance": ";".join("" if v.get("allele_balance") is None else str(round(v["allele_balance"], 6)) for v in variants),
                "dp": ";".join("" if v.get("dp") is None else str(int(v["dp"])) for v in variants),
                "gq": ";".join("" if v.get("gq") is None else str(int(v["gq"])) for v in variants),
                "max_freq_pct": ";".join("" if v["max_freq_pct"] is None
                                         else f"{v['max_freq_pct']:.12g}" for v in variants),
                "max_freq_source": ";".join(v["max_freq_source"] for v in variants),
                "whitelist": ";".join(str(int(v["whitelist"])) for v in variants),
                "frequency_evidence_status": item["frequency_evidence_status"],
                "variant_score_lower_bound": f"{item['variant_score_lower_bound']:.12f}",
                "variant_score_upper_bound": f"{item['variant_score_upper_bound']:.12f}",
            })


def _population_max_from_vep(entries):
    """Return maximum ancestry-specific gnomAD AF as percent plus source field.

    Aggregate gnomAD exome/genome AF is used only when no ancestry-specific
    value is available, so a populated aggregate field is never mistaken for
    missing frequency evidence.
    """
    specific = []
    aggregate = []
    for entry in entries:
        for key, raw in entry.items():
            if key in ("gnomADe_AF", "gnomADg_AF"):
                target = aggregate
            elif key.startswith(("gnomADe_", "gnomADg_")) and key.endswith("_AF"):
                target = specific
            else:
                continue
            for token in re.split(r"[&/,]", raw or ""):
                value = number(token)
                if value is None or not 0 <= value <= 1:
                    continue
                target.append((value * 100.0, key))
    pool = specific or aggregate
    return max(pool, default=(None, ""), key=lambda x: (x[0], x[1]))


def _strict_pathogenic_clinvar(entries):
    labels = []
    for entry in entries:
        raw = entry.get("CLIN_SIG", "")
        labels.extend(x.strip().lower() for x in re.split(r"[&|/]", raw) if x.strip())
    if not labels:
        return False
    bad = ("benign", "uncertain", "conflict")
    return any("pathogenic" in x for x in labels) and not any(
        any(token in x for token in bad) for x in labels
    )


def _genomic_hgvs(chrom, pos, ref, alt):
    """Create a genomic HGVS identifier from a normalized VCF allele."""
    chrom = str(chrom).removeprefix("chr")
    pos = int(pos)
    ref, alt = str(ref), str(alt)
    if len(ref) == len(alt) == 1:
        return f"chr{chrom}:g.{pos}{ref}>{alt}"
    # Trim common suffix then prefix while preserving coordinates.
    r, a = ref, alt
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    shift = 0
    while r and a and r[0] == a[0]:
        r, a = r[1:], a[1:]
        shift += 1
    start = pos + shift
    if not r and a:
        left = start - 1
        right = start
        return f"chr{chrom}:g.{left}_{right}ins{a}"
    if r and not a:
        end = start + len(r) - 1
        site = str(start) if start == end else f"{start}_{end}"
        return f"chr{chrom}:g.{site}del"
    end = start + len(r) - 1
    site = str(start) if start == end else f"{start}_{end}"
    if len(r) == len(a) == 1:
        return f"chr{chrom}:g.{start}{r}>{a}"
    return f"chr{chrom}:g.{site}delins{a}"


def _select_vep_compat_entry(entries, profile):
    """Transcript-consistent VEP proxy for Jannovar-style effect selection."""
    if not entries:
        return None
    tiers = []
    for entry in entries:
        mane = bool(entry.get("MANE_SELECT") or entry.get("MANE_PLUS_CLINICAL") or entry.get("MANE"))
        canonical = str(entry.get("CANONICAL", "")).upper() in ("YES", "1")
        tier = 3 if mane else 2 if canonical else 1
        score = effect_fallback_score(entry.get("Consequence", ""), profile)
        impact = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}.get(entry.get("IMPACT", ""), 0)
        tiers.append((tier, score if score is not None else -1.0, impact,
                      entry.get("Feature", ""), entry))
    best_tier = max(x[0] for x in tiers)
    eligible = [x for x in tiers if x[0] == best_tier]
    chosen = max(eligible, key=lambda x: (x[1], x[2], x[3]))
    entry = chosen[-1]
    if chosen[1] >= 0 or is_regulatory_effect(entry.get("Consequence", ""), profile):
        return entry
    return None


def extract_vcf_compatibility_rows(
        vcf_path, sample=None, profile=None, threads=2, build="GRCh38"):
    """Stream a completed VEP VCF into compatibility variant-gene rows.

    Population frequency is the maximum ancestry-specific gnomAD exome/genome AF,
    expressed as percent, with the source population retained.
    """
    import pysam
    profile = profile or load_profile()
    rows = []
    whitelist = set()
    with pysam.VariantFile(str(vcf_path), threads=threads) as vf:
        samples = list(vf.header.samples)
        if sample is None:
            if len(samples) != 1:
                raise ValueError("Explicit sample required for multi-sample VCF")
            sample = samples[0]
        if sample not in samples:
            raise ValueError("Selected sample not found in VCF")
        if "CSQ" not in vf.header.info:
            raise ValueError("VEP CSQ annotation required")
        keys = vf.header.info["CSQ"].description.split("Format: ", 1)[-1].split("|")
        effect_tokens = tuple(profile["pathogenicity_fallback_scores"])
        regulatory_tokens = tuple(profile.get("regulatory_effects", []))
        for record in vf:
            filters = set(record.filter.keys())
            if filters and filters != {"PASS"}:
                continue
            call = record.samples[sample]
            gt = call.get("GT") or ()
            if not gt or not any(a is not None and a > 0 for a in gt):
                continue
            if len(record.alts or ()) != 1:
                continue
            alt = record.alts[0]
            raw_annotations = record.info.get("CSQ", ()) or ()
            # Preserve both calibrated coding/splice effects and regulatory
            # effects that may later receive explicit ReMM/ClinVar evidence.
            if raw_annotations and not any(
                    any(token in item for token in effect_tokens + regulatory_tokens)
                    for item in raw_annotations):
                continue
            annotations = [dict(zip(keys, item.split("|"))) for item in raw_annotations]
            groups = defaultdict(list)
            for entry in annotations:
                symbol = entry.get("SYMBOL") or entry.get("Gene")
                if symbol:
                    groups[symbol].append(entry)
            max_freq_pct, max_freq_source = _population_max_from_vep(annotations)
            carried_gt = ("|" if call.phased else "/").join(
                "." if a is None else str(a) for a in gt
            )
            ad = call.get("AD") or ()
            ad_total = sum((x or 0) for x in ad) if ad else 0
            allele_balance = ((ad[1] or 0) / ad_total
                              if len(ad) > 1 and ad_total else None)
            phase_set = call.get("PID") if "PID" in call else None
            phase_gt = call.get("PGT") if "PGT" in call else None
            for gene, entries in groups.items():
                selected = _select_vep_compat_entry(entries, profile)
                if selected is None:
                    continue
                effect = selected.get("Consequence", "")
                pathogenic = _strict_pathogenic_clinvar(entries)
                if (effect_fallback_score(effect, profile) is None
                        and not is_regulatory_effect(effect, profile)
                        and not pathogenic):
                    continue
                chrom = str(record.contig).removeprefix("chr")
                if pathogenic:
                    whitelist.add((chrom, str(record.pos), record.ref, alt))
                rows.append({
                    "build": build,
                    "chrom": chrom,
                    "pos": str(record.pos),
                    "ref": record.ref,
                    "alt": alt,
                    "hgvs_id": _genomic_hgvs(chrom, record.pos, record.ref, alt),
                    "gene": gene,
                    "effect": effect,
                    "gt": carried_gt,
                    "phase_set": "" if phase_set is None else str(phase_set),
                    "phase_gt": "" if phase_gt is None else str(phase_gt),
                    "dp": call.get("DP"),
                    "gq": call.get("GQ"),
                    "allele_balance": "" if allele_balance is None else allele_balance,
                    "site_mq": (record.info.get("MQ") if "MQ" in record.header.info else None),
                    "site_qd": (record.info.get("QD") if "QD" in record.header.info else None),
                    "site_fs": (record.info.get("FS") if "FS" in record.header.info else None),
                    "site_mq_rank_sum": (record.info.get("MQRankSum") if "MQRankSum" in record.header.info else None),
                    "site_read_pos_rank_sum": (record.info.get("ReadPosRankSum") if "ReadPosRankSum" in record.header.info else None),
                    "max_freq_pct": "" if max_freq_pct is None else max_freq_pct,
                    "max_freq_source": max_freq_source,
                    "clinvar_pathogenic_whitelist": int(pathogenic),
                    "existing_variation": selected.get("Existing_variation", ""),
                    "selected_transcript": selected.get("Feature", ""),
                    "distance_to_transcript": selected.get("DISTANCE", ""),
                    "transcript_selection": (
                        "MANE" if selected.get("MANE_SELECT") or selected.get("MANE_PLUS_CLINICAL") or selected.get("MANE")
                        else "CANONICAL" if str(selected.get("CANONICAL", "")).upper() in ("YES", "1")
                        else "OTHER"
                    ),
                })
    # One variant-gene observation only.
    dedup = {}
    for row in rows:
        key = (row["chrom"], row["pos"], row["ref"], row["alt"], row["gene"])
        old = dedup.get(key)
        if old is None:
            dedup[key] = row
            continue
        new_score = effect_fallback_score(row["effect"], profile) or -1
        old_score = effect_fallback_score(old["effect"], profile) or -1
        if new_score > old_score:
            dedup[key] = row
    return list(dedup.values()), whitelist


def write_variant_rows(path, rows):
    fields = ["build", "chrom", "pos", "ref", "alt", "hgvs_id", "gene",
              "effect", "gt", "phase_set", "phase_gt", "dp", "gq", "allele_balance",
              "site_mq", "site_qd", "site_fs", "site_mq_rank_sum",
              "site_read_pos_rank_sum", "max_freq_pct", "max_freq_source",
              "clinvar_pathogenic_whitelist", "existing_variation", "selected_transcript",
              "distance_to_transcript", "transcript_selection", "remm", "remm_source"]
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

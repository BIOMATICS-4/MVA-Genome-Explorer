"""Fast gene-agnostic technical/genetic QC for MVA Genome Explorer.

No disease gene, candidate gene, variant coordinate, expected answer, leaderboard
result, or competitor result is encoded in this module.

The QC is computed from the patient's own annotated VCF and the existing blind G0
variant-evidence universe. G1 phenotype contribution is re-applied only after the
QC-adjusted G0 scores have been calculated.
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path
import csv, json, math, statistics
from typing import Any, Dict, Optional, Tuple

DEFAULT_CONFIG = {
    "gq_fail": 20.0,
    "gq_review": 30.0,
    "dp_fail": 8.0,
    "dp_review": 12.0,

    "het_ab_fail_low": 0.15,
    "het_ab_fail_high": 0.85,
    "het_ab_review_low": 0.25,
    "het_ab_review_high": 0.75,
    "het_ab_strong_low": 0.30,
    "het_ab_strong_high": 0.70,

    "hom_alt_ab_fail": 0.60,
    "hom_alt_ab_review": 0.80,
    "hom_alt_ab_strong": 0.90,

    "mq_fail": 40.0,
    "qd_fail": 2.0,
    "fs_fail": 60.0,
    "sor_review": 3.0,
    "mq_rank_sum_fail": -12.5,
    "read_pos_rank_sum_fail": -8.0,

    "local_window_bp": 500,
    "tight_window_bp": 100,
    "local_review_count": 25,
    "local_fail_count": 80,
    "tight_review_count": 8,
    "tight_indel_review_count": 4,
    "low_ab_cluster_fraction_review": 0.35,
    "low_ab_cluster_fraction_fail": 0.50,

    "depth_ratio_review_high": 4.0,
    "depth_ratio_review_low": 0.25,
    "depth_ratio_fail_high": 8.0,
    "depth_ratio_fail_low": 0.125,

    "sequence_flank_bp": 25,
    "homopolymer_review": 8,
    "entropy_review": 1.20,

    "penalty": {
        "QC_STRONG": 0.0,
        "QC_PASS": 0.0,
        # Technical evidence was not supplied (for example a public VCF with
        # GT only). Keep this distinct from observed weak/conflicting evidence.
        "QC_UNASSESSED": -2.0,
        "QC_UNKNOWN": -4.0,
        "QC_REVIEW": -12.0,
        "QC_FAIL": -30.0,
    },
}

def load_config(path: Optional[str] = None) -> Dict[str, Any]:
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))
    if path:
        user = json.loads(Path(path).read_text())
        for k, v in user.items():
            if k == "penalty":
                cfg["penalty"].update(v)
            else:
                cfg[k] = v
    return cfg

def parse_variant_id(value: str) -> Optional[Tuple[str, int, str, str]]:
    if not value or str(value).lower() == "nan":
        return None
    parts = str(value).split(":", 4)
    if len(parts) != 5:
        return None
    _, chrom, pos, ref, alt = parts
    try:
        pos_i = int(pos)
    except ValueError:
        return None
    return chrom.removeprefix("chr"), pos_i, ref, alt

def _num(v):
    if v in (None, "", ".", "nan", "NaN"):
        return None
    if isinstance(v, (tuple, list)):
        if not v:
            return None
        v = v[0]
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except Exception:
        return None

def _gt_label(gt):
    if not gt:
        return "UNKNOWN"
    vals = [x for x in gt if x is not None]
    if not vals:
        return "UNKNOWN"
    nonref = sum(1 for x in vals if x > 0)
    if nonref == 0:
        return "HOM_REF"
    if nonref == len(vals):
        return "HOMOZYGOUS_ALT" if len(set(vals)) == 1 else "MULTIALLELIC_ALT"
    return "HETEROZYGOUS"

def _ab_from_sample(sample):
    try:
        ad = sample.get("AD")
    except Exception:
        ad = None
    if not ad:
        return None
    vals = [x for x in ad if x is not None]
    total = sum(vals)
    if total <= 0:
        return None
    return sum(vals[1:]) / total

def _is_pass_filter(keys):
    keys = [str(x) for x in keys]
    return not keys or keys == ["PASS"] or set(keys) <= {"PASS", "."}

def _safe_info(record, key):
    try:
        return record.info.get(key)
    except Exception:
        return None

def _entropy(seq: str):
    seq = "".join(c for c in seq.upper() if c in "ACGT")
    if not seq:
        return None
    out = 0.0
    n = len(seq)
    for b in "ACGT":
        p = seq.count(b) / n
        if p:
            out -= p * math.log2(p)
    return out

def _max_homopolymer(seq: str) -> int:
    best = cur = 0
    prev = ""
    for c in seq.upper():
        if c == prev and c in "ACGT":
            cur += 1
        elif c in "ACGT":
            cur = 1
        else:
            cur = 0
        prev = c
        best = max(best, cur)
    return best

def _variant_from_input_evidence(row):
    """Fallback for ranking exports that carry coordinates only in input_evidence JSON."""
    for name in ("input_evidence", "evidence", "raw_evidence"):
        raw = row.get(name)
        if not raw:
            continue
        try:
            obj = json.loads(raw)
        except Exception:
            continue
        if not isinstance(obj, dict):
            continue
        chrom = obj.get("chrom")
        pos = obj.get("pos")
        ref = obj.get("ref")
        alt = obj.get("alt")
        if chrom not in (None, "") and pos not in (None, "") and ref and alt:
            try:
                return (str(chrom).removeprefix("chr"), int(pos), str(ref), str(alt))
            except Exception:
                pass
    return None

def detect_variant_column(fieldnames, sample_rows):
    """Detect the coordinate column by schema and actual GRCh38-like content."""
    preferred = (
        "variant", "variant_id", "variant_key", "variant_identity",
        "normalized_variant", "canonical_variant", "variant_1"
    )
    names = [x for x in (fieldnames or []) if x]
    for name in preferred:
        if name in names:
            if any(parse_variant_id(r.get(name, "")) for r in sample_rows):
                return name

    best_name, best_hits = None, 0
    for name in names:
        hits = sum(1 for r in sample_rows if parse_variant_id(r.get(name, "")))
        if hits > best_hits:
            best_name, best_hits = name, hits
    return best_name if best_hits else None

def read_variant_evidence(path: Path):
    rows = []
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        for row in reader:
            rows.append(row)

    probe = rows[: min(len(rows), 2000)]
    variant_column = detect_variant_column(fieldnames, probe)

    keys = set()
    for row in rows:
        key = parse_variant_id(row.get(variant_column, "")) if variant_column else None
        if key is None:
            key = _variant_from_input_evidence(row)
        if key:
            keys.add(key)

    if not keys:
        raise RuntimeError(
            "QC20 could not discover any normalized variant coordinates in "
            f"{path}. Columns were: {fieldnames}. "
            "Aborting instead of assigning QC_UNKNOWN to every hypothesis."
        )

    return rows, keys, variant_column

def _build_positions(keys):
    by_chrom = defaultdict(list)
    key_at = defaultdict(list)
    for key in keys:
        chrom, pos, ref, alt = key
        by_chrom[chrom].append(pos)
        key_at[(chrom, pos)].append(key)
    for chrom in by_chrom:
        by_chrom[chrom] = sorted(set(by_chrom[chrom]))
    return by_chrom, key_at

def _init_stats(keys):
    return {
        k: {
            "exact_found": False,
            "duplicate_exact": 0,
            "site_filter": None,
            "QUAL": None, "QD": None, "MQ": None, "FS": None, "SOR": None,
            "MQRankSum": None, "ReadPosRankSum": None,
            "GT": None, "DP": None, "GQ": None, "AD": None, "AB": None,
            "near500_count": 0,
            "near500_lowab_het_count": 0,
            "near100_count": 0,
            "near100_indel_count": 0,
            "local_dps": [],
            "sequence_entropy": None,
            "max_homopolymer": None,
        } for k in keys
    }

def scan_vcf_for_qc(vcf_path: Path, keys, sample_name=None, threads=4,
                    config=None, reference_fasta=None):
    """One streaming pass. Detailed parsing only occurs inside candidate windows."""
    import pysam

    cfg = config or DEFAULT_CONFIG
    by_chrom, key_at = _build_positions(keys)
    stats = _init_stats(keys)

    vf = pysam.VariantFile(str(vcf_path), threads=max(1, int(threads)))
    samples = list(vf.header.samples)
    if not samples:
        raise RuntimeError("Annotated VCF contains no sample.")
    if sample_name:
        if sample_name not in samples:
            raise RuntimeError(f"Sample {sample_name!r} not found. Available: {samples}")
        sample = sample_name
    elif len(samples) == 1:
        sample = samples[0]
    else:
        raise RuntimeError("Multiple VCF samples found; pass --sample.")

    w = int(cfg["local_window_bp"])
    tight = int(cfg["tight_window_bp"])

    for rec in vf:
        chrom = str(rec.chrom).removeprefix("chr")
        pos_list = by_chrom.get(chrom)
        if not pos_list:
            continue

        left = bisect_left(pos_list, rec.pos - w)
        right = bisect_right(pos_list, rec.pos + w)
        if left == right:
            continue

        sample_data = rec.samples[sample]
        gt = sample_data.get("GT")
        carried = bool(gt and any(x is not None and x > 0 for x in gt))
        if not carried:
            continue

        gt_label = _gt_label(gt)
        ab = _ab_from_sample(sample_data)
        dp = _num(sample_data.get("DP"))
        is_indel = bool(rec.alts and any(len(rec.ref) != len(a) for a in rec.alts if a))
        abnormal_het_ab = (
            gt_label == "HETEROZYGOUS"
            and ab is not None
            and (ab < cfg["het_ab_review_low"] or ab > cfg["het_ab_review_high"])
        )

        for cpos in pos_list[left:right]:
            dist = abs(rec.pos - cpos)
            for key in key_at[(chrom, cpos)]:
                st = stats[key]
                st["near500_count"] += 1
                if abnormal_het_ab:
                    st["near500_lowab_het_count"] += 1
                if dp is not None and len(st["local_dps"]) < 128:
                    st["local_dps"].append(dp)
                if dist <= tight:
                    st["near100_count"] += 1
                    if is_indel:
                        st["near100_indel_count"] += 1

        if (chrom, rec.pos) not in key_at:
            continue

        filter_keys = list(rec.filter.keys())
        for alt in rec.alts or ():
            exact_key = (chrom, rec.pos, rec.ref, alt)
            if exact_key not in stats:
                continue
            st = stats[exact_key]
            st["exact_found"] = True
            st["duplicate_exact"] += 1
            st["site_filter"] = ";".join(filter_keys) if filter_keys else "PASS"
            st["QUAL"] = _num(rec.qual)
            for k in ("QD", "MQ", "FS", "SOR", "MQRankSum", "ReadPosRankSum"):
                st[k] = _num(_safe_info(rec, k))
            st["GT"] = gt_label
            st["DP"] = dp
            st["GQ"] = _num(sample_data.get("GQ"))
            try:
                ad = sample_data.get("AD")
                st["AD"] = ",".join("" if x is None else str(x) for x in ad) if ad else None
            except Exception:
                st["AD"] = None
            st["AB"] = ab

    vf.close()

    if reference_fasta:
        try:
            fa = pysam.FastaFile(str(reference_fasta))
            contigs = set(fa.references)
            flank = int(cfg["sequence_flank_bp"])
            for key, st in stats.items():
                chrom, pos, ref, alt = key
                cname = chrom if chrom in contigs else ("chr" + chrom if "chr" + chrom in contigs else None)
                if not cname:
                    continue
                start = max(0, pos - 1 - flank)
                end = pos - 1 + max(len(ref), 1) + flank
                seq = fa.fetch(cname, start, end)
                st["sequence_entropy"] = _entropy(seq)
                st["max_homopolymer"] = _max_homopolymer(seq)
            fa.close()
        except Exception:
            # Optional evidence: do not fail the run if FASTA/context is unavailable.
            pass

    return sample, stats

def classify_variant_qc(key, st, cfg=None):
    cfg = cfg or DEFAULT_CONFIG
    fail, review, unassessed, positive = [], [], [], []

    if not st["exact_found"]:
        return "QC_UNKNOWN", ["EXACT_VCF_RECORD_NOT_FOUND"]

    if st["duplicate_exact"] > 1:
        fail.append("DUPLICATE_NORMALIZED_EVENT")

    filt = (st.get("site_filter") or "PASS").split(";")
    if not _is_pass_filter(filt):
        fail.append("NONPASS_FILTER")

    if st.get("QD") is not None and st["QD"] < cfg["qd_fail"]:
        fail.append("QD_BELOW_THRESHOLD")
    if st.get("MQ") is not None and st["MQ"] < cfg["mq_fail"]:
        fail.append("MQ_BELOW_THRESHOLD")
    if st.get("FS") is not None and st["FS"] > cfg["fs_fail"]:
        fail.append("STRAND_BIAS_FS")
    if st.get("SOR") is not None and st["SOR"] > cfg["sor_review"]:
        review.append("STRAND_BIAS_SOR")
    if st.get("MQRankSum") is not None and st["MQRankSum"] < cfg["mq_rank_sum_fail"]:
        fail.append("MQ_RANK_SUM")
    if st.get("ReadPosRankSum") is not None and st["ReadPosRankSum"] < cfg["read_pos_rank_sum_fail"]:
        fail.append("READ_POSITION_BIAS")

    gq, dp, ab, gt = st.get("GQ"), st.get("DP"), st.get("AB"), st.get("GT")

    if gq is not None:
        if gq < cfg["gq_fail"]:
            fail.append("LOW_GQ")
        elif gq < cfg["gq_review"]:
            review.append("MARGINAL_GQ")
        else:
            positive.append("GQ_SUPPORTED")
    else:
        unassessed.append("GQ_UNAVAILABLE")

    if dp is not None:
        if dp < cfg["dp_fail"]:
            fail.append("LOW_DP")
        elif dp < cfg["dp_review"]:
            review.append("MARGINAL_DP")
        else:
            positive.append("DP_SUPPORTED")
    else:
        unassessed.append("DP_UNAVAILABLE")

    if gt == "HETEROZYGOUS":
        if ab is None:
            unassessed.append("HET_ALLELE_BALANCE_UNAVAILABLE")
        elif ab < cfg["het_ab_fail_low"] or ab > cfg["het_ab_fail_high"]:
            fail.append("HET_ALLELE_BALANCE_CONFLICT")
        elif ab < cfg["het_ab_review_low"] or ab > cfg["het_ab_review_high"]:
            review.append("HET_ALLELE_BALANCE_REVIEW")
        elif cfg["het_ab_strong_low"] <= ab <= cfg["het_ab_strong_high"]:
            positive.append("HET_ALLELE_BALANCE_STRONG")

    elif gt == "HOMOZYGOUS_ALT":
        if ab is None:
            unassessed.append("HOM_ALT_ALLELE_BALANCE_UNAVAILABLE")
        elif ab < cfg["hom_alt_ab_fail"]:
            fail.append("HOM_ALT_ALLELE_BALANCE_CONFLICT")
        elif ab < cfg["hom_alt_ab_review"]:
            review.append("HOM_ALT_ALLELE_BALANCE_REVIEW")
        elif ab >= cfg["hom_alt_ab_strong"]:
            positive.append("HOM_ALT_ALLELE_BALANCE_STRONG")

    n500 = int(st["near500_count"] or 0)
    n100 = int(st["near100_count"] or 0)
    low_het = int(st["near500_lowab_het_count"] or 0)
    frac = (low_het / n500) if n500 else 0.0

    if n500 >= cfg["local_fail_count"] and frac >= cfg["low_ab_cluster_fraction_fail"]:
        fail.append("DENSE_LOW_AB_LOCAL_CLUSTER")
    else:
        if n500 >= cfg["local_review_count"]:
            review.append("LOCAL_VARIANT_CLUSTER")
        if n500 and frac >= cfg["low_ab_cluster_fraction_review"]:
            review.append("LOW_AB_HET_CLUSTER")

    if n100 >= cfg["tight_review_count"]:
        review.append("TIGHT_VARIANT_CLUSTER")
    if int(st["near100_indel_count"] or 0) >= cfg["tight_indel_review_count"]:
        review.append("COMPLEX_INDEL_CONTEXT")

    local_dps = [x for x in st["local_dps"] if x is not None]
    if dp is not None and len(local_dps) >= 8:
        med = statistics.median(local_dps)
        if med > 0:
            ratio = dp / med
            st["depth_local_median"] = med
            st["depth_ratio"] = ratio
            if ratio > cfg["depth_ratio_fail_high"] or ratio < cfg["depth_ratio_fail_low"]:
                fail.append("EXTREME_LOCAL_DEPTH_OUTLIER")
            elif ratio > cfg["depth_ratio_review_high"] or ratio < cfg["depth_ratio_review_low"]:
                review.append("LOCAL_DEPTH_OUTLIER")

    ent = st.get("sequence_entropy")
    hp = st.get("max_homopolymer")
    if ent is not None and ent < cfg["entropy_review"]:
        review.append("LOW_SEQUENCE_COMPLEXITY")
    if hp is not None and hp >= cfg["homopolymer_review"]:
        review.append("HOMOPOLYMER_CONTEXT")

    if fail:
        return "QC_FAIL", sorted(set(fail + review + unassessed))
    if review:
        return "QC_REVIEW", sorted(set(review + unassessed))

    # Missing metrics are absence of evidence, not observed technical weakness.
    # Preserve this as a separate state so sparse public/legacy VCFs can remain
    # provisionally rankable without being promoted to PASS/STRONG.
    if unassessed:
        return "QC_UNASSESSED", sorted(set(unassessed + positive))

    strong_ab = (
        "HET_ALLELE_BALANCE_STRONG" in positive
        or "HOM_ALT_ALLELE_BALANCE_STRONG" in positive
    )
    if strong_ab and "GQ_SUPPORTED" in positive and "DP_SUPPORTED" in positive:
        return "QC_STRONG", sorted(set(positive))
    if positive:
        return "QC_PASS", sorted(set(positive))
    return "QC_UNKNOWN", ["INSUFFICIENT_TECHNICAL_METRICS"]

def _write_csv(path: Path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields, seen = [], set()
        for row in rows:
            for k in row:
                if k not in seen:
                    seen.add(k)
                    fields.append(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        if rows:
            w.writerows(rows)

def _row_variant_key(row, variant_column=None):
    key = parse_variant_id(row.get(variant_column, "")) if variant_column else None
    return key if key is not None else _variant_from_input_evidence(row)

def _row_variant_text(row, key, variant_column=None):
    if variant_column and row.get(variant_column):
        return row.get(variant_column)
    if key:
        chrom, pos, ref, alt = key
        return f"GRCh38:{chrom}:{pos}:{ref}:{alt}"
    return ""

def _first_present(row, names, default=""):
    for name in names:
        if name in row and row.get(name) not in (None, ""):
            return row.get(name)
    return default

def build_variant_qc_rows(evidence_rows, stats, cfg, variant_column=None):
    out, status_by_key, reasons_by_key = [], {}, {}
    for row in evidence_rows:
        key = _row_variant_key(row, variant_column)
        if not key:
            continue
        st = stats.get(key, {})
        status, reasons = classify_variant_qc(key, st, cfg)
        status_by_key[key] = status
        reasons_by_key[key] = reasons
        out.append({
            "gene": row.get("gene", ""),
            "symbol": row.get("symbol", ""),
            "variant": _row_variant_text(row, key, variant_column),
            "original_score": _first_present(row, ("variant_score", "score", "evidence_score"), ""),
            "overall_qc_status": status,
            "qc_penalty": cfg["penalty"][status],
            "qc_reasons": ";".join(reasons),
            "GT": st.get("GT"),
            "DP": st.get("DP"),
            "GQ": st.get("GQ"),
            "AD": st.get("AD"),
            "allele_balance": st.get("AB"),
            "FILTER": st.get("site_filter"),
            "QD": st.get("QD"),
            "MQ": st.get("MQ"),
            "FS": st.get("FS"),
            "SOR": st.get("SOR"),
            "MQRankSum": st.get("MQRankSum"),
            "ReadPosRankSum": st.get("ReadPosRankSum"),
            "near500_count": st.get("near500_count"),
            "near500_lowab_het_count": st.get("near500_lowab_het_count"),
            "near100_count": st.get("near100_count"),
            "near100_indel_count": st.get("near100_indel_count"),
            "depth_local_median": st.get("depth_local_median"),
            "depth_ratio": st.get("depth_ratio"),
            "sequence_entropy": st.get("sequence_entropy"),
            "max_homopolymer": st.get("max_homopolymer"),
        })
    return out, status_by_key, reasons_by_key

def rerank_hypotheses(g0_path: Path, g1_path: Optional[Path],
                      status_by_key, reasons_by_key, cfg):
    with g0_path.open(newline="", encoding="utf-8") as f:
        g0 = list(csv.DictReader(f))

    g1_by_id = {}
    if g1_path and g1_path.exists():
        with g1_path.open(newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                g1_by_id[row.get("hypothesis_id", "")] = row

    ranked = []
    for row in g0:
        v1 = parse_variant_id(row.get("variant_1", ""))
        v2 = parse_variant_id(row.get("variant_2", ""))

        keys = [k for k in (v1, v2) if k]
        statuses = [status_by_key.get(k, "QC_UNKNOWN") for k in keys]
        qc_adjust = sum(cfg["penalty"][s] for s in statuses)

        # Pair-level rule: a proposed two-allele mechanism cannot keep its
        # inheritance bonus if either required allele is technically unresolved.
        inh = _num(row.get("inheritance_score")) or 0.0
        if v2 and any(s in ("QC_REVIEW", "QC_FAIL", "QC_UNKNOWN") for s in statuses):
            qc_adjust -= inh

        raw_g0 = _num(row.get("case_score")) or 0.0
        qc_g0 = max(0.0, raw_g0 + qc_adjust)

        g1 = g1_by_id.get(row.get("hypothesis_id", ""), {})
        g1_contribution = _num(g1.get("G1_contribution")) or 0.0
        qc_g1 = max(0.0, qc_g0 + g1_contribution)

        if "QC_FAIL" in statuses:
            pair_status = "QC_FAIL"
        elif "QC_REVIEW" in statuses:
            pair_status = "QC_REVIEW"
        elif "QC_UNKNOWN" in statuses:
            pair_status = "QC_UNKNOWN"
        elif "QC_UNASSESSED" in statuses:
            pair_status = "QC_UNASSESSED"
        elif statuses and all(s == "QC_STRONG" for s in statuses):
            pair_status = "QC_STRONG"
        else:
            pair_status = "QC_PASS"

        reasons = []
        for key in keys:
            reasons.extend(reasons_by_key.get(key, []))

        item = dict(row)
        item.update({
            "allele1_qc_status": status_by_key.get(v1, "QC_UNKNOWN") if v1 else "NOT_APPLICABLE",
            "allele2_qc_status": status_by_key.get(v2, "QC_UNKNOWN") if v2 else "NOT_APPLICABLE",
            "pair_qc_status": pair_status,
            "pair_qc_reasons": ";".join(sorted(set(reasons))),
            "raw_G0_score": raw_g0,
            "qc_adjustment": qc_adjust,
            "QC_G0_score": qc_g0,
            "G1_contribution_preserved": g1_contribution,
            "QC_G1_score": qc_g1,
        })
        ranked.append(item)

    ranked.sort(
        key=lambda r: (
            -float(r["QC_G1_score"]),
            -float(r["QC_G0_score"]),
            r.get("symbol", ""),
            r.get("hypothesis_id", ""),
        )
    )
    for i, row in enumerate(ranked, 1):
        row["QC_hypothesis_rank"] = i
    return ranked

def best_gene_rows(hypothesis_rows):
    best = {}
    for row in hypothesis_rows:
        gene = row.get("gene") or row.get("symbol")
        if gene not in best:
            best[gene] = row

    rows = []
    for row in best.values():
        rows.append({
            "gene": row.get("gene", ""),
            "symbol": row.get("symbol", ""),
            "best_hypothesis_id": row.get("hypothesis_id", ""),
            "inheritance_model": row.get("inheritance_model", ""),
            "variant_1": row.get("variant_1", ""),
            "variant_2": row.get("variant_2", ""),
            "pair_qc_status": row.get("pair_qc_status", ""),
            "pair_qc_reasons": row.get("pair_qc_reasons", ""),
            "raw_G0_score": row.get("raw_G0_score", 0),
            "qc_adjustment": row.get("qc_adjustment", 0),
            "QC_G0_score": row.get("QC_G0_score", 0),
            "G1_contribution": row.get("G1_contribution_preserved", 0),
            "QC_G1_score": row.get("QC_G1_score", 0),
        })

    rows.sort(
        key=lambda r: (
            -float(r["QC_G1_score"]),
            -float(r["QC_G0_score"]),
            r["symbol"],
        )
    )
    for i, row in enumerate(rows, 1):
        row["QC_gene_rank"] = i
    return rows

def write_results(out_dir: Path, variant_rows, hypothesis_rows, gene_rows, summary):
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(out_dir / "variant_qc.csv", variant_rows)
    _write_csv(out_dir / "hypothesis_qc.csv", hypothesis_rows)
    _write_csv(out_dir / "top10_hypotheses.csv", hypothesis_rows[:10])
    _write_csv(out_dir / "top10_genes.csv", gene_rows[:10])

    # Explicit blind-QC freeze marker before interpretation of the preserved G1 addition.
    freeze = {
        "stage": "QC_G0_FROZEN",
        "gene_specific_qc": False,
        "phenotype_used_to_construct_qc": False,
        "variant_qc_count": len(variant_rows),
        "hypothesis_count": len(hypothesis_rows),
    }
    (out_dir / "QC_G0_FROZEN.json").write_text(json.dumps(freeze, indent=2, sort_keys=True))
    (out_dir / "qc_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

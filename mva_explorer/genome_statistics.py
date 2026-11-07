"""Descriptive genome statistics derived from the completed local VEP annotation.

This module is intentionally non-ranking. It returns aggregate counts only and
never exposes or serves the raw patient VCF.
"""
from __future__ import annotations
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import threading

from fastapi import APIRouter, HTTPException

VERSION = "genome-statistics-0.1.2"
_POOL = ThreadPoolExecutor(max_workers=1)
_JOB = None
_LOCK = threading.Lock()

PRIMARY_CONSEQUENCE_PRIORITY = [
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_region_variant",
    "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant",
    "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "regulatory_region_variant", "intergenic_variant",
]

CODING = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant", "splice_region_variant",
    "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
}


def _write_json(path: Path, value):
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2))
    tmp.replace(path)


def _now():
    return datetime.now(timezone.utc).isoformat()


def _consequence_tokens(annotations):
    out=set()
    for a in annotations:
        for token in str(a.get("Consequence") or "").split("&"):
            token=token.strip()
            if token: out.add(token)
    return out


def primary_consequence(annotations):
    tokens=_consequence_tokens(annotations)
    for term in PRIMARY_CONSEQUENCE_PRIORITY:
        if term in tokens:
            return term
    return sorted(tokens)[0] if tokens else "unannotated"


def primary_location(annotations):
    tokens=_consequence_tokens(annotations)
    if tokens & CODING:
        return "Coding exon / splice"
    if "5_prime_UTR_variant" in tokens:
        return "5′ UTR"
    if "3_prime_UTR_variant" in tokens:
        return "3′ UTR"
    if "non_coding_transcript_exon_variant" in tokens:
        return "Non-coding exon"
    if "intron_variant" in tokens:
        return "Intronic"
    if "upstream_gene_variant" in tokens:
        return "Upstream"
    if "downstream_gene_variant" in tokens:
        return "Downstream"
    if "intergenic_variant" in tokens:
        return "Intergenic"
    if any("regulatory_region" in t or "TF" in t for t in tokens):
        return "Regulatory / other"
    return "Other"


def regulatory_overlays(annotations):
    """Return non-exclusive regulatory feature labels.

    A variant can therefore remain 'Intronic' in its primary location while also
    being counted as 'Enhancer' here.
    """
    labels=set()
    for a in annotations:
        feature_type=str(a.get("Feature_type") or "")
        biotype=str(a.get("BIOTYPE") or "").strip().lower()
        cons=str(a.get("Consequence") or "").lower()
        if feature_type != "RegulatoryFeature" and "regulatory_region" not in cons and "tf_binding" not in cons:
            continue
        b=biotype.replace("_"," ")
        if "enhancer" in b:
            labels.add("Enhancer")
        elif "promoter" in b:
            labels.add("Promoter")
        elif "ctcf" in b:
            labels.add("CTCF binding")
        elif "open chromatin" in b:
            labels.add("Open chromatin")
        elif "tf" in b or "transcription factor" in b:
            labels.add("TF binding")
        else:
            labels.add("Other regulatory")
    return labels


def _population_max(annotations):
    values=[]
    for a in annotations:
        for key,val in a.items():
            if key.startswith(("gnomADe_","gnomADg_")) and key.endswith("_AF"):
                for token in re.split(r"[&/]", str(val or "")):
                    try:
                        n=float(token)
                        if 0 <= n <= 1: values.append(n)
                    except (TypeError,ValueError):
                        pass
    return max(values) if values else None


def rarity_bucket(annotations):
    af=_population_max(annotations)
    if af is None:
        return "Unknown / not found"
    if af <= 0.0001:
        return "Very rare variants ≤0.01%"
    if af <= 0.001:
        return "Rare variants >0.01% to 0.1%"
    if af <= 0.01:
        return "Low-frequency variants >0.1% to 1%"
    return "Common >1%"


def source_identity(out: Path):
    status=out/"annotation_status.json"
    marker=out/"vep.done.json"
    if not status.exists():
        return None
    state=json.loads(status.read_text())
    sha=None
    if marker.exists():
        try: sha=json.loads(marker.read_text()).get("outputs",{}).get("annotated.vcf.gz")
        except Exception: pass
    return {
        "annotation_status":state.get("status"),
        "input_sha256":state.get("identity",{}).get("input_sha256"),
        "annotation_sha256":sha,
        "annotation_version":state.get("version"),
        "sample":state.get("identity",{}).get("sample"),
    }


def build_summary(out: Path):
    import pysam
    annotated=out/"annotated.vcf.gz"
    if not annotated.exists():
        raise ValueError("Completed annotated VCF not found.")
    identity=source_identity(out)
    if not identity or identity.get("annotation_status") != "ANNOTATION_COMPLETE":
        raise ValueError("Genome annotation is not complete.")

    status_path=out/"genome_statistics_status.json"
    summary_path=out/"genome_statistics.json"
    status={"status":"RUNNING","stage":"Scanning VEP annotations","percent":0.0,
            "records_scanned":0,"updated":_now(),"source":identity,"version":VERSION}
    _write_json(status_path,status)

    locations=Counter(); consequences=Counter(); regulatory=Counter(); rarity=Counter()
    records=0
    with pysam.VariantFile(str(annotated)) as vf:
        if "CSQ" not in vf.header.info:
            raise ValueError("Annotated VCF does not contain VEP CSQ annotations.")
        keys=vf.header.info["CSQ"].description.split("Format: ",1)[-1].split("|")
        total_hint=None
        ann_status=out/"annotation_status.json"
        try:
            st=json.loads(ann_status.read_text())
            total_hint=st.get("annotation_summary",{}).get("carried_annotated_alleles") or st.get("annotated_records_exported")
        except Exception:
            pass
        for records,rec in enumerate(vf,1):
            anns=[dict(zip(keys,item.split("|"))) for item in rec.info.get("CSQ",())]
            locations[primary_location(anns)] += 1
            consequences[primary_consequence(anns)] += 1
            rarity[rarity_bucket(anns)] += 1
            for label in regulatory_overlays(anns): regulatory[label] += 1
            if records % 100000 == 0:
                pct=min(99.5, (records/total_hint*100.0) if total_hint else 0.0)
                status.update(percent=round(pct,1),records_scanned=records,updated=_now())
                _write_json(status_path,status)

    ref_audit={}
    try: ref_audit=json.loads((out/"reference_audit.json").read_text())
    except Exception: pass
    summary={
        "version":VERSION,"generated_utc":_now(),"source":identity,
        "scope":{
            "annotated_reference_compatible_records":records,
            "input_records_seen":ref_audit.get("records_seen"),
            "reference_compatible_records":ref_audit.get("reference_compatible_records"),
            "unresolved_or_reference_incompatible_records":sum(v for k,v in ref_audit.items() if k not in {"records_seen","carried_records","reference_compatible_records","no_called_alt_in_selected_sample"} and isinstance(v,int))
        },
        "primary_location":dict(locations),
        "primary_consequence":dict(consequences),
        "regulatory_overlays":dict(regulatory),
        "rarity":dict(rarity),
    }
    _write_json(summary_path,summary)
    _write_json(status_path,{"status":"READY","stage":"Complete","percent":100.0,
                             "records_scanned":records,"updated":_now(),"source":identity,"version":VERSION})
    return summary


def genome_statistics_router(root: Path):
    router=APIRouter()
    out=root/"analysis/genome_only"

    @router.get("/api/genome-statistics/status")
    def status():
        identity=source_identity(out)
        sp=out/"genome_statistics_status.json"
        summary=out/"genome_statistics.json"
        if summary.exists():
            try:
                data=json.loads(summary.read_text())
                if data.get("source",{}).get("annotation_sha256") == (identity or {}).get("annotation_sha256"):
                    return {"status":"READY","stage":"Complete","percent":100.0,"source":identity,"version":VERSION}
            except Exception:
                pass
        if sp.exists():
            try:return json.loads(sp.read_text())
            except Exception: pass
        return {"status":"NOT_STARTED","stage":"Waiting","percent":0.0,"source":identity,"version":VERSION}

    @router.post("/api/genome-statistics/start", status_code=202)
    def start():
        global _JOB
        identity=source_identity(out)
        if not identity or identity.get("annotation_status") != "ANNOTATION_COMPLETE":
            raise HTTPException(409,"Complete genome annotation first.")
        current=status()
        if current.get("status") == "READY": return current
        with _LOCK:
            if _JOB is not None and not _JOB.done(): return {**current,"status":"RUNNING"}
            def work():
                try:return build_summary(out)
                except Exception as exc:
                    _write_json(out/"genome_statistics_status.json",{"status":"FAILED","stage":"Failed","error":str(exc),"updated":_now(),"source":identity,"version":VERSION})
                    raise
            _JOB=_POOL.submit(work)
        return {"status":"RUNNING","stage":"Scanning VEP annotations","percent":0.0,"source":identity,"version":VERSION}

    @router.get("/api/genome-statistics/summary")
    def summary():
        if status().get("status") != "READY":
            raise HTTPException(409,"Genome Statistics is not ready yet.")
        p=out/"genome_statistics.json"
        if not p.exists(): raise HTTPException(404)
        return json.loads(p.read_text())

    return router

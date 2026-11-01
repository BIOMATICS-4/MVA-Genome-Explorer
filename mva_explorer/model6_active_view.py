"""Read-only promoted Model 6 v1.2 research ranking and variant-backed details.

This does not overwrite frozen compatibility v1.2 ranking artifacts. It serves a
version-labelled ACTIVE research revision from atomic, input-hash-protected
snapshot data. Fail closed when source inputs are modified or missing.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from fastapi import APIRouter, HTTPException

RELEASE_LABEL = "Model 6 v1.2.0"
METHOD_NOTE = (
    "Research hypothesis prioritisation only. Unresolved recessive phase is "
    "not a confirmed diagnosis. Missing population AF and incomplete QC "
    "must not be interpreted as evidence of rarity or validity."
)

def _read_tsv(path):
    with path.open(newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))

def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def _keys(row):
    return [x for x in str(row.get("variant_keys") or "").split(";") if x]

def _load_active(root):
    directory = root / "analysis/exomiser_compatibility/current/model6_active_v12"
    path = directory / "manifest.json"
    if not path.is_file():
        raise HTTPException(404, "Promoted Model 6 result is not installed.")
    data = json.loads(path.read_text())
    if data.get("display_version") != "1.2.0" or data.get("active") is not True:
        raise HTTPException(409, "Model 6 active release metadata invalid.")
    project_root = root.parent
    for relative, digest in data.get("verify_input_sha256", {}).items():
        src = project_root / relative
        if (not src.is_file()) or _sha(src) != digest:
            raise HTTPException(409, f"Model 6 source changed: {relative}; recalculate first.")
    for filename, digest in data.get("verify_result_sha256", {}).items():
        file = directory / filename
        if (not file.is_file()) or _sha(file) != digest:
            raise HTTPException(409, f"Model 6 active result file corrupted: {filename}")
    return directory, data

def _view_rows(directory, mode):
    filename = "exploratory_all_genes.tsv" if mode == "exploratory" else "restrictive_all_genes.tsv"
    rows = _read_tsv(directory / filename)
    if len(rows)<10 or any(int(r["rank"]) != i+1 for i,r in enumerate(rows)):
        raise HTTPException(409, "Model 6 stored gene ranks have gaps or invalid ordering.")
    if len({r["gene"] for r in rows}) != len(rows):
        raise HTTPException(409, "Model 6 stored ranking duplicates a gene.")
    return rows

def _qc_map(root):
    base=root/"analysis/exomiser_compatibility/current"
    result={}
    p=base/"compatibility_variant_qc.tsv"
    if p.is_file():
        for row in _read_tsv(p):
            key=row.get("variant_key")
            if key and row.get("qc_status"):
                result[key]=row["qc_status"]
    p=base/"regulatory_variant_qc.tsv"
    if p.is_file():
        for row in _read_tsv(p):
            key=":".join(str(row.get(k,"")) for k in ("build","chrom","pos","ref","alt"))
            if row.get("qc_status") and key not in result:
                result[key]=row["qc_status"]
    return result

def _variant_map(root):
    base=root/"analysis/exomiser_compatibility/current"
    result={}
    p=base/"variant_gene.tsv"
    for row in _read_tsv(p):
        key=":".join(str(row.get(k,"")).removeprefix("chr") if k=="chrom" else str(row.get(k,""))
                     for k in ("build","chrom","pos","ref","alt"))
        result.setdefault((row.get("gene"),key),row)
    return result

def model6_active_router(root):
    router=APIRouter()

    @router.get("/api/compatibility/model6-current-status")
    def current_status():
        _,meta=_load_active(root)
        return dict(display_name=RELEASE_LABEL,
            display_version=meta["display_version"],
            research_mode=meta["default_ranking_mode"],
            scoring_revision=meta["scoring_revision"],
            restrictions=METHOD_NOTE,
            installation=meta["installed_at"])

    @router.get("/api/compatibility/model6-current-top10")
    def current_top10(mode: str="exploratory"):
        if mode not in ("exploratory","restrictive"):
            raise HTTPException(400,"Mode must be exploratory or restrictive")
        directory,meta=_load_active(root)
        ranked=_view_rows(directory,mode)
        qc=_qc_map(root)
        output=[]
        for row in ranked[:10]:
            keys=_keys(row)
            if not keys:
                raise HTTPException(409,f"No variants recorded for {row['gene']}")
            phase=row.get("phase") or "UNRESOLVED"
            clinical_status=("AR_TWO_VARIANTS_PHASE_UNRESOLVED" if
                row.get("moi")=="AR" and len(keys)>1 and phase=="UNRESOLVED"
                else "COMPUTATIONAL_HYPOTHESIS_NOT_DIAGNOSIS")
            output.append({
                "gene_rank":int(row["rank"]),
                "gene":row["gene"],
                "ranked_hypothesis_id":row["hypothesis"],
                "ranked_moi":row["moi"],
                "ranked_score":row["legacy_combined"],
                "ranked_variant_score":row.get("variant_score_upper_bound",""),
                "ranked_phenotype_score":row["gene_phenotype"],
                "ranked_phase":phase,
                "ranked_status":clinical_status,
                "ranked_variant_keys":";".join(keys),
                "ranked_variant_count":len(keys),
                "ranked_effects":row.get("effects",""),
                "ranked_qc":";".join(qc.get(k,"QC_UNASSESSED") for k in keys),
                "disease_id":row.get("disease_id",""),
                "disease_name":row.get("disease",""),
                "known_supported_gene_moi":row.get("moi",""),
                "moi_evidence_caveat":METHOD_NOTE+" "+row.get("warnings",""),
                "phase_review_required":phase=="UNRESOLVED",
                "ranking_mode":mode.upper(),
                "display_model":RELEASE_LABEL,
                "internal_revision":meta["scoring_revision"]
            })
        return output

    @router.get("/api/compatibility/model6-current-inheritance-details/{gene}")
    def current_inheritance_details(gene: str, mode: str="exploratory"):
        if mode not in ("exploratory","restrictive"):
            raise HTTPException(400,"Mode must be exploratory or restrictive")
        directory,_=_load_active(root)
        candidates=_view_rows(directory,mode)[:10]
        row=next((r for r in candidates if r["gene"]==gene),None)
        if row is None:
            raise HTTPException(404,"Gene is not in the current Model 6 Top 10.")
        raw=_variant_map(root)
        qc=_qc_map(root)
        keys=_keys(row)
        effects=(row.get("effects") or "").split(";")
        gt=(row.get("genotype") or "").split(";")
        variants=[]
        for i,key in enumerate(keys):
            variant=raw.get((gene,key))
            if variant is None:
                raise HTTPException(409,"Selected variant does not exist in source candidate rows.")
            actual=variant.get("gt") or ""
            claimed=gt[i] if i<len(gt) else ""
            if claimed and actual and claimed!=actual:
                raise HTTPException(409,"Selected genotype disagrees with source input.")
            variants.append({
                "key":key,"gt":actual or claimed,
                "effect":effects[i] if i<len(effects) else variant.get("effect",""),
                "qc":qc.get(key,"QC_UNASSESSED"),
                "source":"ORIGINAL_CACHED_VARIANT_GENE_ROW"
            })
        return {
            "gene":gene,
            "ranked_moi":row["moi"],
            "ranked_score":float(row["legacy_combined"]),
            "ranked_keys":keys,
            "clinical_moi_proven":False,
            "ranking_mode":mode.upper(),
            "ranked_phase":row.get("phase",""),
            "ranked_disease_id":row.get("disease_id",""),
            "ranked_disease":row.get("disease",""),
            "hypotheses":[{
                "moi":row["moi"],
                "phase":row.get("phase",""),
                "score":float(row["legacy_combined"]),
                "keys":keys,
                "variants":variants,
                "disease_id":row.get("disease_id",""),
                "status":row.get("lane_label",""),
                "warnings":row.get("warnings","")
            }],
            "individual_ad":[],
            "ar":[],
            "variant_info":variants,
            "notes":METHOD_NOTE
        }
    return router

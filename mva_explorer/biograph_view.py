"""Bridge to the independently frozen local BioGraph Engine.

BioGraph is an adjunct interpretation service. It never feeds scores back into
MVA ranking, QC, inheritance, or submission outputs.
"""
from __future__ import annotations

import json
import os
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

from fastapi import APIRouter, HTTPException, Query

BIOGRAPH_URL = os.environ.get("MVA_BIOGRAPH_URL", "http://127.0.0.1:8000").rstrip("/")


def _get(path: str, params: dict | None = None, timeout: float = 30.0):
    url = BIOGRAPH_URL + path
    if params:
        url += "?" + urlencode(params)
    req = Request(url, headers={"Accept": "application/json"})
    try:
        with urlopen(req, timeout=timeout) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", "replace")
        try:
            detail = json.loads(detail).get("detail", detail)
        except Exception:
            pass
        raise HTTPException(exc.code, str(detail)) from exc
    except (URLError, TimeoutError, ConnectionError) as exc:
        raise HTTPException(
            503,
            "BioGraph Engine is not available. Start the integrated workspace "
            "with: bash start_integrated.sh",
        ) from exc


def _enc(value: str) -> str:
    return quote(value, safe=":-_.~")


def biograph_router():
    router = APIRouter(prefix="/api/knowledge", tags=["biological-context"])

    @router.get("/status")
    def status():
        try:
            root = _get("/", timeout=2.0)
            return {
                "status": "CONNECTED",
                "mode": "POST_RANKING_ADJUNCT",
                "role": "POST_RANKING_INTERPRETATION_ONLY",
                "ranking_effect": "NONE",
                "ranking_influence": "NONE",
                "biograph_url": BIOGRAPH_URL,
                "engine": root,
            }
        except HTTPException:
            return {
                "status": "OFFLINE",
                "mode": "POST_RANKING_ADJUNCT",
                "role": "POST_RANKING_INTERPRETATION_ONLY",
                "ranking_effect": "NONE",
                "ranking_influence": "NONE",
                "biograph_url": BIOGRAPH_URL,
                "detail": (
                    "BioGraph Engine is not available. "
                    "Start the integrated workspace with: bash start_integrated.sh"
                ),
            }

    @router.get("/search")
    def search(q: str = Query(...), hpo_only: bool = False):
        return _get("/search", {"q": q, "hpo_only": str(hpo_only).lower()})

    @router.get("/genes/search")
    def gene_search(q: str = Query(...)):
        return _get("/genes/search", {"q": q})

    @router.get("/disease/{disease_id}/expand")
    def disease_expand(disease_id: str):
        return _get(f"/disease/{_enc(disease_id)}/expand")

    @router.get("/disease/{disease_id}/phenotypes")
    def disease_phenotypes(disease_id: str):
        return _get(f"/disease/{_enc(disease_id)}/phenotypes")

    @router.get("/phenotype/{phenotype_id}/expand")
    def phenotype_expand(phenotype_id: str):
        return _get(f"/phenotype/{_enc(phenotype_id)}/expand")

    @router.get("/go/{go_id}/expand")
    def go_expand(go_id: str):
        return _get(f"/go/{_enc(go_id)}/expand")

    @router.get("/gene/{gene_id}/expand")
    def gene_expand(gene_id: str):
        return _get(f"/gene/{_enc(gene_id.upper())}/expand")

    @router.get("/gene/{gene_id}/go_terms")
    def gene_go_terms(gene_id: str):
        return _get(f"/gene/{_enc(gene_id.upper())}/go_terms")

    @router.get("/gene/{gene_id}/diseases")
    def gene_diseases(gene_id: str):
        return _get(f"/gene/{_enc(gene_id.upper())}/diseases")

    @router.get("/pathway/{pathway_id}/expand")
    def pathway_expand(pathway_id: str, disease_id: str | None = None):
        params = {"disease_id": disease_id} if disease_id else None
        return _get(f"/pathway/{_enc(pathway_id)}/expand", params)

    @router.get("/gene/{gene_id}/metadata")
    def gene_metadata(gene_id: str):
        return _get(f"/gene/{_enc(gene_id.upper())}/metadata")

    @router.get("/pathway/{pathway_id}/metadata")
    def pathway_metadata(pathway_id: str):
        return _get(f"/pathway/{_enc(pathway_id)}/metadata")

    @router.get("/disease/{disease_id}/metadata")
    def disease_metadata(disease_id: str):
        return _get(f"/disease/{_enc(disease_id)}/metadata")

    @router.get("/phenotype/{phenotype_id}/metadata")
    def phenotype_metadata(phenotype_id: str):
        return _get(f"/phenotype/{_enc(phenotype_id)}/metadata")

    @router.get("/go/{go_id}/metadata")
    def go_metadata(go_id: str):
        return _get(f"/go/{_enc(go_id)}/metadata")

    @router.get("/gene/{gene_id}/context")
    def gene_context(gene_id: str):
        gene = gene_id.upper()
        diseases = _get(f"/gene/{_enc(gene)}/diseases")
        go_terms = _get(f"/gene/{_enc(gene)}/go_terms")
        pathways = _get(f"/gene/{_enc(gene)}/expand")
        return {
            "gene": gene,
            "role": "POST_RANKING_INTERPRETATION_ONLY",
            "ranking_effect": "NONE",
            "diseases": diseases.get("nodes", []),
            "go_terms": go_terms.get("nodes", []),
            "pathways": pathways.get("nodes", []),
            "counts": {
                "diseases": diseases.get("total_count", 0),
                "go_terms": go_terms.get("total_count", 0),
                "pathways": pathways.get("total_count", 0),
            },
        }

    return router

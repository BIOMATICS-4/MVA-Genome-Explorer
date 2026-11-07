"""Read-only validation summary API for the unified MVA workspace."""
from pathlib import Path
import json

from fastapi import APIRouter, HTTPException


def validation_router(project_root: Path):
    router = APIRouter(prefix="/api/validation", tags=["validation"])
    validation_root = project_root / "benchmarks" / "exomiser_compatibility"

    @router.get("/summary")
    def summary():
        path = validation_root / "cross_benchmark_validation_20261007.json"
        if not path.exists():
            raise HTTPException(404, "Cross-benchmark validation summary is not available.")
        return json.loads(path.read_text())

    @router.get("/pank2")
    def pank2():
        path = (
            validation_root
            / "confirmed_case_pank2_20261007"
            / "mva_v12_spikein"
            / "validation_summary.json"
        )
        if not path.exists():
            raise HTTPException(404, "PANK2 validation summary is not available.")
        return json.loads(path.read_text())

    return router

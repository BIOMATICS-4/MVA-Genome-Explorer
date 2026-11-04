"""Data-free regression checks for Model 6 v1.2 saved-result provenance.

No gated genomic records or patient-level details are needed.
"""
import csv
import hashlib
import json

from fastapi import FastAPI
from fastapi.testclient import TestClient

from mva_explorer.model6_active_view import model6_active_router


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(tmp_path):
    root = tmp_path / "data"
    directory = root / "analysis/exomiser_compatibility/current/model6_active_v12"
    directory.mkdir(parents=True)
    fields = [
        "gene", "rank", "hypothesis", "moi", "phase", "variant_keys",
        "effects", "genotype", "legacy_combined", "gene_phenotype",
        "variant_score_upper_bound", "QC_status", "warnings",
    ]
    p = directory / "exploratory_all_genes.tsv"
    with p.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for i in range(1, 11):
            writer.writerow({
                "gene": f"TESTGENE{i}", "rank": i,
                "hypothesis": f"TESTGENE{i}_AR", "moi": "AR",
                "phase": "UNRESOLVED", "variant_keys": f"TEST:{i}:A:T;TEST:{i}:C:G",
                "effects": "missense_variant;missense_variant",
                "genotype": "0/1;0/1", "legacy_combined": str(1-i/20),
                "gene_phenotype": str(1-i/15),
                "variant_score_upper_bound": "0.7",
                "QC_status": "QC_INDEPENDENTLY_SUPPORTED",
                "warnings": "PHASE_UNRESOLVED",
            })
    source = root / "source_record.txt"
    source.write_text("CHANGED AFTER FREEZE")
    manifest = {
        "display_version": "1.2.0", "active": True,
        "scoring_revision": "frozen-test", "installed_at": "2026-10-18",
        "default_ranking_mode": "EXPLORATORY",
        "verify_input_sha256": {"data/source_record.txt": hashlib.sha256(b"ORIGINAL").hexdigest()},
        "verify_result_sha256": {p.name: _sha(p)},
    }
    (directory / "manifest.json").write_text(json.dumps(manifest))
    app = FastAPI()
    app.include_router(model6_active_router(root))
    return TestClient(app), p


def test_saved_results_visible_but_changed_live_source_fails_closed(tmp_path):
    client, _ = _fixture(tmp_path)
    current = client.get("/api/compatibility/model6-current-top10")
    assert current.status_code == 409
    assert "source changed" in current.text
    saved = client.get("/api/compatibility/model6-saved-top10")
    assert saved.status_code == 200
    rows = saved.json()
    assert len(rows) == 10
    assert rows[0]["archived_snapshot"] is True
    assert rows[0]["gene"] == "TESTGENE1"
    assert rows[0]["ranked_phase"] == "UNRESOLVED"
    assert "HISTORICAL SNAPSHOT" in rows[0]["moi_evidence_caveat"]
    details = client.get("/api/compatibility/model6-saved-inheritance-details/TESTGENE1")
    assert details.status_code == 200
    assert details.json()["archived_snapshot"] is True
    assert all(v["source"].startswith("FROZEN_MODEL6") for v in details.json()["variant_info"])


def test_saved_results_reject_modified_frozen_output(tmp_path):
    client, p = _fixture(tmp_path)
    p.write_text(p.read_text() + "TAMPERED\n")
    assert client.get("/api/compatibility/model6-saved-top10").status_code == 409
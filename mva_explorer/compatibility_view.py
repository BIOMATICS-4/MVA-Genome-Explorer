"""Read-only API and page for the isolated Exomiser compatibility analysis."""
import csv
import html
import json
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, HTMLResponse, Response
import io
from .model6_consistency import read_tsv, build_consistent_rows, qc_lookup
from .model6_moi_context import load_gencc, audit_gene_moi
from .model6_hypothesis_details import detail as inheritance_detail
from .exomiser_compatibility import load_profile, read_variant_metadata


def compatibility_router(root):
    router = APIRouter()
    out = root / "analysis/exomiser_compatibility/current"

    def _read_tsv(name):
        path = out / name
        if not path.exists():
            return []
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def _gene_top10():
        rows = _read_tsv("full_ranking.tsv")
        seen, genes = set(), []
        for row in rows:
            if row["gene"] in seen:
                continue
            seen.add(row["gene"])
            genes.append({
                "gene_rank": len(genes) + 1,
                "gene": row["gene"],
                "best_moi": row["moi"],
                "hypothesis_rank": row["rank"],
                "combined_score": row["combined_score"],
                "phenotype_score": row["phenotype_score"],
                "variant_score": row["variant_score"],
                "phase_status": row.get("phase_status", ""),
                "variant_keys": row["variant_keys"],
                "effects": row["effects"],
                "max_freq_pct": row["max_freq_pct"],
                "max_freq_source": row["max_freq_source"],
            })
            if len(genes) == 10:
                break
        return genes

    def _table(rows):
        if not rows:
            return "No completed ranking yet."
        headers = list(rows[0])
        value = "<table><tr>" + "".join(
            f"<th>{html.escape(str(k))}</th>" for k in headers
        ) + "</tr>"
        for row in rows:
            cells = []
            for k in headers:
                raw = str(row.get(k, ""))
                shown = html.escape(raw)
                if k == "gene" and raw:
                    shown = f'<a target="_top" href="/?tab=evidence&amp;gene={html.escape(raw, quote=True)}">{shown}</a>'
                cells.append(f"<td>{shown}</td>")
            value += "<tr>" + "".join(cells) + "</tr>"
        return value + "</table>"

    @router.get("/api/exomiser-compatibility/status")
    @router.get("/api/compatibility/status")
    def status():
        p = out / "summary.json"
        return json.loads(p.read_text()) if p.exists() else {"status":"NOT_RUN"}

    @router.get("/api/exomiser-compatibility/top10")
    @router.get("/api/compatibility/top10")
    def top10():
        rows = _read_tsv("top10.tsv")
        if not rows:
            raise HTTPException(404, "Compatibility ranking has not been run.")
        return rows

    @router.get("/api/exomiser-compatibility/top10-genes")
    @router.get("/api/compatibility/top10-genes")
    def top10_genes():
        rows = _gene_top10()
        if not rows:
            raise HTTPException(404, "Compatibility ranking has not been run.")
        return rows

    @router.get("/api/exomiser-compatibility/top10-qc-strict")
    @router.get("/api/compatibility/top10-qc-strict")
    def top10_qc_strict():
        rows = _read_tsv("qc_strict_top10_genes.tsv")
        if not rows:
            raise HTTPException(404, "QC-strict compatibility view has not been generated.")
        return rows

    @router.get("/api/exomiser-compatibility/model6-consistent-top10")
    @router.get("/api/compatibility/model6-consistent-top10")
    def model6_consistent_top10():
        """Live score-faithful view: one score refers to one exact genotype hypothesis.

        Always recompute from the current input artifacts so the Results page
        cannot silently display a projection from an earlier patient run.
        """
        required = ("qc_strict_top10_genes.tsv", "full_ranking.tsv",
                    "compatibility_variant_qc.tsv")
        if any(not (out / name).exists() for name in required):
            raise HTTPException(404, "Model 6 strict ranking or QC is not available.")
        qc_rows = read_tsv(out / "compatibility_variant_qc.tsv")
        regulatory = out / "regulatory_variant_qc.tsv"
        if regulatory.exists():
            for row in read_tsv(regulatory):
                row["variant_key"] = ":".join(
                    row[column] for column in ("build", "chrom", "pos", "ref", "alt")
                )
                qc_rows.append(row)
        try:
            results = build_consistent_rows(
                read_tsv(out / "qc_strict_top10_genes.tsv"),
                read_tsv(out / "full_ranking.tsv"), qc_rows
            )
            gencc = root / 'standards/gencc-submissions.csv'
            gencc_data = load_gencc(gencc) if gencc.exists() else {}
            for row in results:
                row.update(audit_gene_moi(row['gene'], row['ranked_moi'], gencc_data))
            return results
        except ValueError as exc:
            raise HTTPException(409, f"Model 6 evidence mismatch: {exc}") from exc

    @router.get("/api/compatibility/model6-inheritance-details/{gene}")
    def model6_inheritance_details(gene: str):
        """Read-only, on-demand scores. No ranking mutation or MOI diagnosis."""
        strict = next((r for r in _read_tsv("qc_strict_top10_genes.tsv")
                       if r.get("gene") == gene), None)
        if strict is None:
            raise HTTPException(404, "Gene absent from current strict Top 10")
        mandatory = ("variant_gene.tsv", "predictor_metadata.tsv", "full_ranking.tsv",
                     "compatibility_variant_qc.tsv")
        if any(not (out / name).exists() for name in mandatory):
            raise HTTPException(404, "Model 6 source inputs unavailable")
        raw = [r for r in _read_tsv("full_ranking.tsv") if r.get("gene") == gene]
        with (out / "variant_gene.tsv").open(newline="") as handle:
            variants = [r for r in csv.DictReader(handle, delimiter="\t")
                        if r.get("gene") == gene]
        # The two technical QC sources overlap. Never let the last file read
        # silently override a disagreeing QC classification for one allele.
        qc_records = _read_tsv("compatibility_variant_qc.tsv")
        extra = out / "regulatory_variant_qc.tsv"
        if extra.exists():
            # Regulatory QC is locus-level: its column is best_gene, not gene.
            # Restrict to loci actually present for this requested gene.
            requested_keys = {":".join(v[k] for k in ("build", "chrom", "pos", "ref", "alt"))
                              for v in variants}
            for r in read_tsv(extra):
                key = ":".join(r[k] for k in ("build", "chrom", "pos", "ref", "alt"))
                if key in requested_keys and r.get("qc_status"):
                    qc_records.append({"variant_key": key, "qc_status": r["qc_status"]})
        try:
            qcs = qc_lookup(qc_records)
            return inheritance_detail(gene, strict, raw, variants, qcs,
                                      read_variant_metadata(out / "predictor_metadata.tsv"),
                                      load_profile())
        except (ValueError, KeyError) as e:
            raise HTTPException(409, f"Score provenance mismatch: {e}") from e

    @router.get("/api/exomiser-compatibility/submission-top10")
    @router.get("/api/compatibility/submission-top10")
    def submission_top10():
        rows = _read_tsv("submission_top10_v13.tsv")
        if not rows:
            raise HTTPException(404, "Audited submission representation has not been generated.")
        return rows

    @router.get("/api/exomiser-compatibility/populations")
    @router.get("/api/compatibility/populations")
    def populations():
        p = out / "population_summary.json"
        if not p.exists():
            raise HTTPException(404, "Population audit has not been generated.")
        return json.loads(p.read_text())

    @router.get("/api/exomiser-compatibility/qc-audit")
    @router.get("/api/compatibility/qc-audit")
    def qc_audit():
        p = out / "qc_audit.json"
        if not p.exists():
            raise HTTPException(404, "Compatibility QC audit has not been generated.")
        return json.loads(p.read_text())

    @router.get("/api/exomiser-compatibility/download/{name}")
    def download(name: str):
        allowed = {
            "summary.json","top10.tsv","full_ranking.tsv","variant_gene.tsv",
            "population_summary.json","change_report.json","ar_candidate_qc.tsv",
            "qc_strict_ranking.tsv","qc_strict_top10_genes.tsv","submission_top10_v13.tsv","qc_strict_summary.json",
            "compatibility_variant_qc.tsv","regulatory_variant_qc.tsv","qc_audit.json",
            "linkage_v1_summary.json","linkage_v12_summary.json",
            "model6_consistent_top10.tsv",
        }
        if name not in allowed:
            raise HTTPException(404)
        if name == "model6_consistent_top10.tsv":
            # Generate from current inputs, never serve stale results from a
            # previous run of the same project with a different patient.
            rows = model6_consistent_top10()
            stream = io.StringIO()
            if rows:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            return Response(
                content=stream.getvalue(), media_type="text/tab-separated-values",
                headers={"Content-Disposition": 'attachment; filename="model6_consistent_top10.tsv"'},
            )
        p = out / name
        if not p.exists():
            raise HTTPException(404)
        return FileResponse(p, filename=name)

    @router.get("/compatibility", response_class=HTMLResponse)
    def page():
        s = status()
        raw_rows = _gene_top10() if s.get("status") == "COMPLETE" else []
        qc_rows = _read_tsv("qc_strict_top10_genes.tsv")
        try:
            submission_rows = model6_consistent_top10()
        except HTTPException as exc:
            if exc.status_code != 404:
                raise  # A QC/evidence conflict must not silently disappear.
            submission_rows = []
        pop_path = out / "population_summary.json"
        pop = json.loads(pop_path.read_text()) if pop_path.exists() else {}
        links = " ".join(
            f'<a href="/api/exomiser-compatibility/download/{name}">{name}</a>'
            for name in (
                "qc_strict_top10_genes.tsv","submission_top10_v13.tsv","qc_strict_ranking.tsv","qc_strict_summary.json",
                "compatibility_variant_qc.tsv","regulatory_variant_qc.tsv","top10.tsv","full_ranking.tsv",
                "ar_candidate_qc.tsv","qc_audit.json","linkage_v1_summary.json","linkage_v12_summary.json",
                "population_summary.json","summary.json","change_report.json",
                "model6_consistent_top10.tsv"
            )
            if (out / name).exists() or (name == "model6_consistent_top10.tsv" and bool(submission_rows))
        )
        return f"""<!doctype html><html><meta charset="utf-8"><title>Exomiser compatibility</title>
<style>body{{font:15px system-ui;background:#f3f6fa;color:#000;max-width:1200px;margin:40px auto;padding:24px}}section{{background:white;padding:22px;border-radius:12px;margin:18px 0;overflow:auto}}.warn{{background:#fff0cf;padding:14px;border-radius:8px}}table{{border-collapse:collapse;width:100%;font-size:13px}}td,th{{padding:7px;border-bottom:1px solid #ddd;text-align:left}}pre{{white-space:pre-wrap}}a{{color:#000;text-decoration:underline;margin-right:12px}}</style>
<h1>MVA Exomiser-compatible ranking</h1>
<p class="warn">Research compatibility profile only. Not ACMG classification or a clinical diagnosis. The default MVA standards-first ranking remains separate.</p>
<section><h2>QC-strict Top 10 genes</h2><p>Current v1.2.0 model: Exomiser-compatible biological prioritization plus preserved MVA fast_qc technical QC. Regulatory/noncoding loci use raw ReMM evidence multiplied by a separate gene-linkage confidence based on selected MANE/canonical transcript overlap or VEP distance; raw ReMM is retained for audit. Two-allele AR hypotheses require confirmed trans phase. QC_REVIEW/FAIL/UNKNOWN calls are excluded from this strict headline view.</p>{_table(qc_rows)}</section>
<section><h2>Evidence-consistent v1.2 result</h2><p>The ranked hypothesis score and displayed contributing variants always come from the SAME hypothesis. Any alternative recessive pair has its OWN score and phase flag. An AD-style singleton label is a computational ranking scenario, not an established disease inheritance mechanism or diagnosis. Unresolved phase remains unresolved. The archived v13 export is preserved in Downloads solely for auditing.</p>{_table(submission_rows)}</section>
<section><h2>Raw compatibility Top 10 genes</h2><p>Preserved for transparent comparison, including phase-unresolved hypotheses. For regulatory calls, raw ReMM, linkage confidence, VEP distance and the linkage-weighted effective score are retained separately in the downloadable ranking.</p>{_table(raw_rows)}</section>
<section><h2>Status</h2><pre>{html.escape(json.dumps(s,indent=2))}</pre></section>
<section><h2>Population-frequency audit</h2><pre>{html.escape(json.dumps(pop,indent=2))}</pre></section>
<section><h2>Downloads</h2>{links}</section></html>"""

    return router

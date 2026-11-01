"""Read-only view of the separate annotation command's progress and exports."""
import json
from pathlib import Path
from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse, HTMLResponse

def annotation_router(root):
    router=APIRouter()
    out=root/'analysis/genome_only'

    @router.get('/api/annotation/status')
    def status():
        p=out/'annotation_status.json'
        return json.loads(p.read_text()) if p.exists() else {'status':'NOT_STARTED','instruction':'Run setup_annotation.sh, then the annotation command.'}

    @router.get('/api/annotation/download/{name}')
    def download(name:str):
        allowed={'annotation_status.json','reference_audit.json','annotation_summary.json','review_candidates.csv.gz','gene_inventory.csv','other_rare_or_unknown.csv.gz'}
        if name not in allowed:raise HTTPException(404)
        state=status()
        if state.get('status')!='ANNOTATION_COMPLETE':raise HTTPException(409,'Run is not complete; exports may be partial.')
        path=out/name
        if not path.exists():raise HTTPException(404)
        return FileResponse(path,filename=name)

    @router.get('/annotation',response_class=HTMLResponse)
    def page():
        return '''<!doctype html><html lang="en"><meta charset="utf-8"><title>Genome annotation</title>
        <style>body{font:16px system-ui;background:#f3f6fa;color:#192c41;max-width:1100px;margin:40px auto;padding:24px}pre{white-space:pre-wrap;background:white;padding:24px;border-radius:12px}a{display:inline-block;margin:8px;color:#135f82}</style>
        <h1>Genome annotation</h1><p>Live status of the local reference-checking, normalization and offline VEP run.</p>
        <p>These outputs are annotations and review candidates. Inheritance assessment and causal ranking are still pending.</p>
        <pre id="status">Loading…</pre><div id="exports"></div>
        <script>async function refresh(){try{const r=await fetch('/api/annotation/status');if(!r.ok)throw Error('Status request failed');const s=await r.json();document.getElementById('status').textContent=JSON.stringify(s,null,2);const box=document.getElementById('exports');box.replaceChildren();if(s.status==='ANNOTATION_COMPLETE'){for(const n of ['annotation_status.json','reference_audit.json','annotation_summary.json','review_candidates.csv.gz','gene_inventory.csv','other_rare_or_unknown.csv.gz']){let a=document.createElement('a');a.href='/api/annotation/download/'+n;a.textContent=n;box.append(a)}}}catch(e){document.getElementById('status').textContent=e.message}}refresh();setInterval(refresh,5000)</script></html>'''
    return router

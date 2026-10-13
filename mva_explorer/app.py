import json
import os
import sqlite3
import uuid
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from fastapi import FastAPI, HTTPException
from fastapi.responses import HTMLResponse, FileResponse, RedirectResponse
from fastapi.middleware.trustedhost import TrustedHostMiddleware
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from .core import qc

ROOT = Path(os.environ.get('MVA_DATA', 'data')).resolve()
ROOT.mkdir(parents=True, exist_ok=True)
RUNS = ROOT/'runs'; RUNS.mkdir(exist_ok=True)
app = FastAPI(title='MVA Genome Explorer — QC prototype')
app.add_middleware(TrustedHostMiddleware, allowed_hosts=['127.0.0.1','localhost','testserver'])
app.add_middleware(CORSMiddleware, allow_origins=['http://127.0.0.1:3000','http://localhost:3000'], allow_credentials=False, allow_methods=['GET','OPTIONS'], allow_headers=['*'])
pool = ThreadPoolExecutor(max_workers=1)
jobs = {}

class Input(BaseModel):
    filename: str
    sample: str | None = None

def run_path(run):
    if not all(c in '0123456789abcdef' for c in run) or len(run)!=32:
        raise HTTPException(404)
    p = RUNS/run
    if not p.exists(): raise HTTPException(404)
    return p

@app.get('/', response_class=HTMLResponse)
def home(): return (Path(__file__).parent/'index.html').read_text()

@app.get('/api/files')
def files():
    return [p.name for p in ROOT.iterdir() if p.is_file() and p.name.endswith(('.vcf','.vcf.gz','.bcf'))]

@app.post('/api/qc', status_code=202)
def start(body: Input):
    p = (ROOT/body.filename).resolve()
    if p.parent != ROOT or not p.is_file(): raise HTTPException(400,'Choose a VCF from the data directory.')
    if any(not f.done() for f in jobs.values()): raise HTTPException(409,'A run is already active.')
    run = uuid.uuid4().hex
    out = RUNS/run; out.mkdir()
    def work():
        try: return qc(p, out, body.sample or None)
        except Exception as e:
            (out/'error.json').write_text(json.dumps({'error':str(e)}))
    jobs[run] = pool.submit(work)
    return {'run':run}

@app.get('/api/runs')
def runs():
    # Oldest -> newest so the UI can reliably select the most recent run.
    dirs=[p for p in RUNS.iterdir() if p.is_dir()]
    return [p.name for p in sorted(dirs, key=lambda x: x.stat().st_mtime)]

@app.get('/api/runs/{run}')
def status(run: str):
    p=run_path(run)
    for name in ('error.json','run_manifest.json'):
        if (p/name).exists(): return json.loads((p/name).read_text())
    return {'status':'QUEUED_OR_HASHING'}

@app.get('/api/runs/{run}/alleles')
def alleles(run: str, offset:int=0, limit:int=20):
    p=run_path(run)
    if offset<0 or limit<1 or limit>200: raise HTTPException(400)
    if not (p/'variants.sqlite').exists(): return []
    with sqlite3.connect(p/'variants.sqlite') as db:
        db.row_factory=sqlite3.Row
        # The UI polls while QC is still creating the database. A SQLite file can
        # therefore exist briefly before the alleles table has been created.
        # Treat that transient state as 'no rows yet' instead of a server error.
        ready = db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='alleles'").fetchone()
        if not ready:
            return []
        return [dict(x) for x in db.execute('SELECT * FROM alleles WHERE carried=1 ORDER BY rowid LIMIT ? OFFSET ?', (limit, offset))]


@app.get('/api/runs/{run}/qc-summary')
def qc_summary(run: str):
    p=run_path(run)/'02_variant_qc_summary.csv'
    if not p.exists(): return []
    import csv
    with p.open(newline='') as f:
        return list(csv.DictReader(f))

@app.get('/api/runs/{run}/download/{name}')
def download(run: str, name: str):
    if name not in ('run_manifest.json','02_variant_qc_summary.csv'): raise HTTPException(404)
    p=run_path(run)/name
    if not p.exists(): raise HTTPException(404)
    return FileResponse(p, filename=name)

from .uploads import upload_router
app.include_router(upload_router(ROOT))

from .annotation_view import annotation_router
app.include_router(annotation_router(ROOT))

from .compatibility_view import compatibility_router
app.include_router(compatibility_router(ROOT))

from .biograph_view import biograph_router
app.include_router(biograph_router())

from .validation_view import validation_router
app.include_router(validation_router(Path(__file__).resolve().parents[1]))

from .acmg_preview import acmg_router
app.include_router(acmg_router(ROOT))

from .genome_statistics import genome_statistics_router
app.include_router(genome_statistics_router(ROOT))

@app.get("/knowledge", include_in_schema=False)
@app.get("/knowledge/", include_in_schema=False)
def knowledge_redirect():
    # Evidence Graph now lives inside the unified MVA application shell.
    return RedirectResponse(url="/?tab=evidence", status_code=307)

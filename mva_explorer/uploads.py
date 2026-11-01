
import asyncio
import os
import re
import uuid
from pathlib import Path
import pysam
from fastapi import APIRouter, HTTPException, Request

def upload_router(root):
    router = APIRouter()
    limit = int(os.environ.get("MVA_MAX_UPLOAD_BYTES", 4 * 1024**3))

    @router.post("/api/upload", status_code=201)
    async def upload(request: Request, filename: str):
        if (not filename or len(filename) > 180 or "/" in filename
                or "\\" in filename or any(ord(c) < 32 for c in filename)):
            raise HTTPException(400, "Invalid filename.")

        suffix = (
            ".vcf.gz" if filename.lower().endswith(".gz")
            else Path(filename).suffix.lower()
        )
        if suffix not in (".vcf", ".vcf.gz", ".bcf"):
            raise HTTPException(
                400, "Select a VCF, compressed VCF, or BCF; not an index or ZIP."
            )

        length = request.headers.get("content-length")
        if length and (not length.isdecimal() or int(length) > limit):
            raise HTTPException(413, "Upload exceeds the configured size limit.")

        safe = re.sub(r"[^A-Za-z0-9_-]", "_", filename)[:80]
        token = uuid.uuid4().hex
        dest = root / (safe + "_" + token + suffix)
        temp = root / (".upload-" + token)
        size = 0

        try:
            with temp.open("xb") as f:
                async for chunk in request.stream():
                    size += len(chunk)
                    if size > limit:
                        raise HTTPException(
                            413, "Upload exceeds the configured size limit."
                        )
                    await asyncio.to_thread(f.write, chunk)

            if not size:
                raise HTTPException(400, "The selected file is empty.")

            def probe():
                with pysam.VariantFile(str(temp)) as vf:
                    samples = list(vf.header.samples)
                    next(vf, None)
                    return samples

            try:
                samples = await asyncio.to_thread(probe)
            except Exception:
                raise HTTPException(
                    400,
                    "Cannot read this as VCF/BCF. Check the file and compression."
                )

            temp.replace(dest)
            return {
                "filename": dest.name,
                "bytes": size,
                "samples": samples,
                "validation":
                    "Header and first record readable; run QC to check the full file."
            }
        except OSError:
            raise HTTPException(
                507, "Cannot save the upload. Check disk space and folder permissions."
            )
        finally:
            temp.unlink(missing_ok=True)

    return router

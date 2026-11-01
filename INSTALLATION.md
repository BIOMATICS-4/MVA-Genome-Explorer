# Installation and reproduction

**Status: partial research-code release. The full application cannot yet be independently installed from this repository alone.**

To obtain current files:

```bash
git clone https://github.com/BIOMATICS-4/MVA-Genome-Explorer.git
cd MVA-Genome-Explorer
```

If the repository is private, access permissions are required.

For preliminary Python setup on Ubuntu / WSL with Python 3.12:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pytest tests/test_core.py tests/test_inheritance.py tests/test_acmg_preview.py
```

These tests have not yet been confirmed to pass from a clean clone, and may depend on additional modules or fixtures. Full application installation requires the remaining backend/UI source, startup scripts, separately developed BioGraph API and React frontend, reference resource acquisition, privacy-safe synthetic example data, and end-to-end validation.

The original local application used `bash install.sh` followed by `bash start.sh` with MVA Explorer usually at http://127.0.0.1:8010. Those commands are **not yet supported by this partial GitHub package**.

Do not publicly distribute patient VCFs or restricted reference databases.

## Updated public-code snapshot (1 November 2026)

The `mva_explorer/` folder now contains the locally installed Model 6 v1.2 backend route modules, evidence/sex-aware/mechanism-ranking modules, versioned API adapters and current `index.html` browser UI. The original local targeted checks passed **43 tests** (1 warning). The public repository still does **not** include restricted patient inputs, frozen model input/output matrices, complete BioGraph service or every runtime reference dependency. Successful local tests are not proof of a clean end-to-end reproduction from the public checkout alone. See [Model 6 release reconciliation](docs/MODEL6_V1_2_RECONCILIATION.md).

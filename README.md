# MVA Genome Explorer — MVA Rare Disease Hackathon

Research software and methods supporting two complementary submissions:

- **Track 1: Variant prioritisation** — genomic annotation, technical QC, inheritance-aware candidate assessment, transparent ranking and ACMG/AMP evidence review.
- **Track 2: Mechanism-informed drug repurposing** — analysis of BUB1B/BUBR1 variants, network-based biological interpretation and literature-derived pharmacological hypotheses.

## Interpretation and limitations

The BUB1B p.Leu737Ter variant is a predicted loss-of-function allele; p.Asn1002Lys is a missense VUS. The two heterozygous variants have unresolved phase. Proposed interventions, including premature-stop readthrough, NAD+-pathway support and 4-PBA, remain hypotheses: rescue efficacy and safety have not been experimentally established. BioGraph is a separate post-ranking interpretation service and does not modify prioritisation scores. The ACMG/AMP-style evidence map highlights criteria requiring expert review and is not a clinical classifier.

## Reproducibility and data safety

Only task-relevant, reviewed code, methods, configuration and tests should be committed. Patient VCF files, identifiable data, giant reference datasets, caches, virtual environments, backup versions and incidental local results are intentionally excluded. Reference resources and authorised example input datasets will be needed to reproduce full-genome runs. No genomic dataset is included in this initial repository documentation.

## AI assistance

Parts of software development and writing involved AI assistance. See each track's methods disclosure for provider, plan/tier and applicable data-handling settings.

## Submission package status

Core research Python modules, preliminary tests, track-specific overview files and dependency requirements have been uploaded. The full original methods workbooks, reference resources, complete BioGraph source and end-to-end reproduction bundle are **not yet uploaded**. This is a partial research-code release, not yet a verified end-to-end reproduction package.

For local tests, create a Python 3.12 virtual environment, install `requirements.txt`, and run `python -m pytest tests/test_core.py tests/test_inheritance.py tests/test_acmg_preview.py`. Full application startup will require further source modules and external resources.

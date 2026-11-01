# Model 6 v1.2.0 — application, submission and hypothesis reconciliation

**Model label is fixed:** Model 6 v1.2.0 (short name Model 6 v1.2). Do not rename the scientific model because of an internal scoring implementation or a release audit.

## Important: two different recorded Top-10 sets

The local application's installed active `model6_active_v12` snapshot (manifest date 18 October 2026; implementation `model6_mechanism_v41.py`; scoring revision `evidence-combined-20261018`) displays in its default mode:

1. BUB1B
2. ZFP57
3. PIK3C2A
4. PKD1
5. ADAMTS1
6. STK36
7. TTN
8. LAMA3
9. GHRHR
10. PEG3

The separate final consolidated six-stage Track 1 report describes its final v1.2 ranking as:

1. BUB1B
2. ADAMTS1
3. PIK3C2A
4. BRINP1
5. ARL9
6. TRIB1
7. NPRL3
8. FOXP1
9. GNRHR
10. ACTB

A saved CSV named `UscBen45_Track1_Model6_v1.2.csv` also lists the first (ZFP57...) set, while the public competition leaderboard reports a later submission filename `model3_UscBen45_Track1_Model6_v1.2.csv`. **The content bytes of that exact portal-uploaded file have not been independently recovered.** Therefore these sources must NOT be merged, renamed as equivalent, or used interchangeably. The app's active default and the final consolidated report are separately labelled until the portal submission is reconciled.

## Genetic conclusion and reporting requirements

- Leading BUB1B hypothesis: p.Leu737Ter (predicted premature stop, ClinVar P/LP evidence) and p.Asn1002Lys (missense VUS).
- Both are heterozygous and technically supported; phase is **UNRESOLVED**, not confirmed in trans.
- In prior final reporting, individual BUB1B p.Leu737Ter AD score 0.979006 and candidate AR pair score 0.882026 were **different constructs**. Do not substitute one for the other or describe these research scores as clinical pathogenicity probabilities.
- Missing population-frequency evidence means UNKNOWN, not rare by default. Never equate unverified QC with QC_PASS.

## Track 2 — submitted research hypotheses

- Niacin: NAD⁺/SIRT2/BUBR1 abundance mechanism — a hypothesis based on indirect precedent, not patient-specific rescue.
- 4-phenylbutyrate (4-PBA): conditional proteostasis/folding or protein-abundance rescue hypothesis.
- Gentamicin: premature-stop readthrough research hypothesis, with uncertainty and safety constraints.
- Amlexanox: premature-stop readthrough/NMD-related hypothesis; pharmacological and clinical feasibility constraints remain.
- NMN: **mechanistic comparator only**, not a fifth primary shortlisted treatment.

None of these has been established as effective or safe for this individual's MVA. Existing genomic, cell and literature work must remain labelled at its actual evidence level.

## Public release status

This repository is a partial source-code release. A full data-free reproduction of the currently installed application, including frontend, runtime assets, reference snapshot and complete model-output checking, has not been independently established by this audit. Never commit patient-level VCF/FASTQ/BAM, clinical documents, private intermediate tables or protected reference datasets.

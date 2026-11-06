# Track 1 — Genomic variant prioritisation

**Objective:** Rank plausible candidate genes and variants from whole-genome variation using variant consequence, population frequency, sample-level quality control, phenotype evidence and inheritance-aware hypotheses.

**Key interpretation:** The BUB1B p.Leu737Ter predicted loss-of-function allele and p.Asn1002Lys missense VUS are both heterozygous. Phase remains unresolved; detection of two variants does not establish a confirmed compound-heterozygous diagnosis.

**Analysis components:** GRCh38 reference and VEP annotation, ClinVar/gnomAD/dbSNP evidence, variant-level QC, inheritance-aware scoring, ReMM × gene linkage and post-ranking BioGraph interpretation. Predictions are not clinical classifications.

## Model 6 v1.2 — reported Top 10

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

This is the single saved ranking shown on the repository homepage and in the [saved Model 6 v1.2 CSV](results/UscBen45_Track1_Model6_v1.2_saved_2026-10-22.csv). The later Track 1 report (`UscBen45_Track1_variant_prioritisation (report).pdf`) includes a Results screenshot supporting the saved order, with BUB1B first and ADAMTS1 fifth. The exact competition-portal CSV is not available here for hash verification; **do not claim that this saved file is a verified copy of the portal upload**.

[Results provenance and supporting records](../../docs/MODEL6_V1_2_PROVENANCE.md).

**Reproduction status:** The repository contains a partial, data-free application snapshot and tests. The full genomic dataset, reference resources and an independently verified end-to-end reproduction are not included.

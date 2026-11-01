# Track 1 — Genomic variant prioritisation

**Objective:** Rank plausible candidate genes and variants from whole-genome variation using variant consequence, population rarity, sample-level quality control, phenotype evidence and inheritance-aware hypotheses.

**Key interpretation:** The BUB1B p.Leu737Ter predicted loss-of-function allele and p.Asn1002Lys missense VUS are both heterozygous. Phase remains unresolved; detection of two variants does not establish a compound-heterozygous diagnosis.

**Analysis components:** GRCh38 alignment and VEP annotation; ClinVar, gnomAD and dbSNP evidence; genomic QC, inheritance modelling and reproducible scoring; ACMG/AMP evidence flags for human review. Predictions are not clinical classifications.

**Reproduction status:** The full reproducible genomic code and dataset-independent verification scripts are not yet included. Patient data and large annotation caches are intentionally excluded. The methods workbook should be added separately after a data-handling review.

## Version and Top 10 provenance

**Model 6 v1.2.0** remains the application model name. The saved [22 October Model 6 v1.2 CSV](results/UscBen45_Track1_Model6_v1.2_saved_2026-10-22.csv) contains the Top 10 **BUB1B, ZFP57, PIK3C2A, PKD1, ADAMTS1, STK36, TTN, LAMA3, GHRHR, PEG3** and corresponds by gene order to the installed active ranking snapshot recorded on 18 October 2026. The six-stage consolidated methods report independently describes **BUB1B, ADAMTS1, PIK3C2A, BRINP1, ARL9, TRIB1, NPRL3, FOXP1, GNRHR, ACTB** as its final table. See [the reconciliation](../../docs/MODEL6_V1_2_RECONCILIATION.md); the exact portal-uploaded final CSV bytes have not been recovered, so neither of these should be relabelled as the verified official final submission.

The BUB1B candidate pair remains **phase UNRESOLVED**. Research scores and quality-control labels are separate from clinical pathogenicity classifications.

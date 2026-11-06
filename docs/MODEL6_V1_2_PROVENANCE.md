# Model 6 v1.2 — reported Track 1 Top 10 provenance

## Single public Top 10

The following ten-gene order is the **saved Model 6 v1.2 research result** used on the repository homepage:

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

## Supporting records

- **Later Track 1 methods/report:** `UscBen45_Track1_variant_prioritisation (report).pdf` (researcher's saved report). Its Results screenshot (Figure 4) shows the leading portion of this ranking, including **BUB1B rank 1** and **ADAMTS1 rank 5**. The report also describes the BUB1B autosomal-recessive candidate pair and unresolved phase.
- **Saved ten-row results CSV:** [`UscBen45_Track1_Model6_v1.2_saved_2026-10-22.csv`](../submissions/track1/results/UscBen45_Track1_Model6_v1.2_saved_2026-10-22.csv). This gives the complete ten-gene order shown above and uses research submission-format confidence values. It is **not independently verified as the exact competition-portal uploaded file**.
- **Application Results screenshot:** [Model 6 v1.2 Results interface](../assets/model6/model6-v12-results-clarified.png). This is visual corroboration of the saved ranking interface.
- **Software snapshot and tests:** [`mva_explorer/`](../mva_explorer/) and [saved-snapshot test](../tests/test_model6_saved_snapshot.py) document the retained implementation context. These cannot replace a byte-level comparison with the original portal submission.

## Interpretation boundaries

BUB1B was ranked first with a **0.882026** combined research prioritisation score for its candidate autosomal-recessive pair, p.Leu737Ter and p.Asn1002Lys. Both variants were heterozygous, and their phase remained **UNRESOLVED**. The separate single-allele value **0.979006** must not be substituted for the pair score. Combined scores are not calibrated clinical pathogenicity probabilities; the ten ranks are not simply sorted by descending combined score.

**Submission verification limit:** The authentic uploaded competition-portal CSV and its hash have not been recovered in this workspace. The order above is supported by the saved CSV and later Track 1 report, but it must not be described as a byte-for-byte authenticated upload. No competing Top 10 is presented as a current result.

This repository is a partial, data-free software release. Patient-level genomic data, protected clinical documentation and private intermediate datasets are intentionally excluded.

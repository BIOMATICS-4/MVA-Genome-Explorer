# Model 6 v1.2.0 — application, submission and hypothesis reconciliation


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
- NMN: **mechanistic comparator only**.



## Public release status

This repository is a partial source-code release. A full data-free reproduction of the currently installed application, including frontend, runtime assets, reference snapshot and complete model-output checking, has not been independently established by this audit. Never commit patient-level VCF/FASTQ/BAM, clinical documents, private intermediate tables or protected reference datasets.


## 4 November 2026 clarification — saved Model 6 v1.2 research results

Two distinct saved candidate lists must **not** be silently substituted for one another:

| Record | Top 10 order |
|---|---|
| Six-stage consolidated comparison dated 8 October 2026 | BUB1B, ADAMTS1, PIK3C2A, BRINP1, ARL9, TRIB1, NPRL3, FOXP1, GNRHR, ACTB |
| Frozen 18 October active Model 6 v1.2 snapshot and 19 October audited candidate CSV | BUB1B, ZFP57, PIK3C2A, PKD1, ADAMTS1, STK36, TTN, LAMA3, GHRHR, PEG3 |

The **later Track 1 methods report**, supplied separately by the researcher, explicitly gives BUB1B an AR-pair score of **0.882026** and ADAMTS1 rank **5** (Figure 2/Table 3), and its Results screenshot shows the second order above. The audited 19 October CSV was validated against that frozen result. The uploaded competition-portal CSV is not available here for hash comparison; the reported 22 October submission filename/leaderboard entry does **not** prove byte identity.

**Score context, not a new model version:** the Results view uses the *selected ranked hypothesis's* combined score (**BUB1B AR pair = 0.882026; phase UNRESOLVED**). An earlier AD-style single-variant BUB1B score (**p.Leu737Ter = 0.979006**) appeared in an ACMG/AMP review screenshot. It is **not** the two-variant score, a separate ACMG clinical rating, or proof of a later model. The interface now labels the two score contexts separately.

**Ranking is not descending combined-score sorting.** For the stored Model 6 v1.2 evidence-led ordering, the predicted effect class and gene-level phenotype evidence precede lower-priority tie-breaks such as QC-review status, established gene–disease support and conservative variant-score lower bound. Consequently ADAMTS1 (combined score ≈0.9228) can rank below BUB1B (AR-pair score ≈0.8820). Those 0–1 scores are research heuristics, **not calibrated disease probabilities**. The candidate confidence (EPCR) column in the audit CSV is a rank-based submission encoding, not this logistic score.

**Provenance safeguard:** the live-result endpoint continues to require matching hashes for both source inputs and frozen outputs. When a source input has changed, a separate explicitly **saved/historical** view is available only if frozen result hashes still verify; it must not be described as a newly recalculated or live-revalidated ranking. Neither the historic CSV nor the frozen ranking files were rewritten to make these records agree.


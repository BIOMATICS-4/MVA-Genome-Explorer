# MVA Genome Explorer — Hackathon submissions

## Model 6 v1.2 — Report-supported Top 10 and ranking scores

| Rank | Gene | Disease / phenotype context | **← Phenotype match (P)** | Disease–HPO support | Variant evidence (V) | Combined score (S) |
|---:|---|---|---:|---:|---:|---:|
| **1** | **BUB1B** | Mosaic variegated aneuploidy syndrome 1 | **← 0.7657** | 0.743605 | 0.799321 | 0.882026 |
| 2 | ZFP57 | Diabetes mellitus, transient neonatal, 1 | **← 0.6721** | 0.216122 | 0.364518 | 0.049533 |
| 3 | PIK3C2A | Oculocerebrodental syndrome | **← 0.6648** | 0.453376 | 0.445816 | 0.092493 |
| 4 | PKD1 | Autosomal dominant polycystic kidney disease | **← 0.6494** | 0.178863 | 0.741125 | 0.566705 |
| 5 | ADAMTS1 | Autosomal dominant prognathism | **← 0.6346** | — | 0.998801 | 0.922801 |
| 6 | STK36 | Primary ciliary dyskinesia | **← 0.6176** | 0.071677 | 0.599528 | 0.203829 |
| 7 | TTN | Tibial muscular dystrophy | **← 0.5847** | 0.132194 | 0.596020 | 0.149730 |
| 8 | LAMA3 | Unresolved disease hypothesis | **← 0.5766** | — | 0.600000 | 0.143764 |
| 9 | GHRHR | Unresolved disease hypothesis | **← 0.5741** | — | 0.593442 | 0.133474 |
| 10 | PEG3 | Unresolved disease hypothesis | **← 0.5695** | — | 0.840055 | 0.585756 |

![Model 6 v1.2 Results interface](assets/model6/model6-v12-results-clarified.png)

*Figure 1. Report-supported Model 6 v1.2 gene ranking. This is the ten-gene order recorded in the saved results CSV and reflected in the later Track 1 report screenshot; exact identity with the competition-portal upload cannot be independently certified. Phenotype-match scores (P, marked ← above) help explain the ranking; the combined scores (S) are separate metrics.*

### How the scores were calculated

**Phenotype match (P):** Imported Exomiser/HiPhive-style gene–phenotype similarity. The arrows in the table identify the score that determined the ordering **within this particular Top 10**, because all ten candidates share the same broad predicted coding/essential-splice effect class.

**Variant evidence (V):** For coding/splice alleles, the higher of the consequence-based support weight and available computational predictor score, multiplied by the population-frequency multiplier. For a two-variant AR hypothesis, V is the mean of the two allele scores. Regulatory evidence uses ReMM × gene-linkage confidence. The table shows the saved V upper-bound value; missing population frequency does not prove rarity.

**Combined hypothesis score (S):** Both P and V enter the existing logistic formula:

`Z = −13.28813 + 10.39451 × P + 9.18381 × V`

`S = 1 / (1 + exp(−Z))`

For **BUB1B**, P = **0.7657** and V ≈ **0.799321** yield **S = 0.882026** for its AR candidate pair (phase **UNRESOLVED**).

**How rank was assigned:** Model 6 v1.2 first prioritises predicted effect class and then P; QC-review status, established gene–disease evidence, conservative variant support and disease–HPO support act as subsequent ordering criteria. **S is a separate combined-evidence score, not the value used to sort these saved ranks.** Hence ADAMTS1's S = 0.922801 can exceed BUB1B's S = 0.882026 while ADAMTS1 remains rank 5.

*The disease-context column lists a stored hypothesis, not confirmed diagnosis. Missing disease–HPO support is shown as — rather than zero. These scores are research prioritisation values, not calibrated pathogenicity probabilities.*

## Explore the submissions

| Track | Methods | Results |
|---|---|---|
| **Track 1 — Variant prioritisation** | — | — |
| **Track 2 — Mechanism-informed drug repurposing** | — | — |

## BioGraph — visual exploration

[![BUB1B disease network](assets/biograph/bub1b-semantic-disease-network.png)](assets/biograph/README.md)

[**See all three BioGraph screenshots (Classic pathways, Semantic diseases, and BUB3 GO terms) →**](assets/biograph/README.md)


## Reproducibility and data sharing

[Installation status](INSTALLATION.md) · [Data availability and restrictions](DATA_AVAILABILITY.md)

The repository currently contains a partial software release, not a validated one-command installation. No open-source license is granted. Do not distribute patient-level genomic files without appropriate authorization.

---

## [Future Prospects — GENE-KNOW-DRUG (GKG)](future-prospects/)

Inspired by the **Hackathon**, I have decided to develop **Drug Cascader**, an interactive pharmacology tool designed to explore and eventually predict what happens to a drug after administration, from its absorption and distribution throughout the body to its interactions with molecular targets, biological pathways, tissues and organs.

The current prototype allows users to interactively explore drug-associated molecular relationships and their connections to broader biological systems. Future development aims to simulate a drug's journey through the body, predicting **which tissues and organs it may reach, how much of the drug may accumulate there, how long it may persist, and the likelihood of beneficial or harmful effects**.

The goal is to move beyond identifying where a drug acts, towards estimating **tissue-specific drug concentrations, exposure trajectories, downstream molecular responses and the probability of different biological outcomes over time**. These predictions will be informed by pharmacokinetic models, molecular interactions, tissue-expression data and established biological evidence, with uncertainty explicitly represented.

Drug Cascader will be integrated with **MVA Genome Explorer** and **BioGraph** to form **GENE-KNOW-DRUG (GKG)**, a broader research software initiative extending beyond MVA, with the ambition of investigating rare genetic diseases and exploring potential therapeutic strategies.

**GKG remains a research-driven ambition, with Drug Cascader currently under development.**

### Drug Cascader — Current Prototype

<!-- Screenshot will be added when the prototype image is available. -->

*Current development-stage interface of Drug Cascader.*

---


Research software and methods supporting two complementary submissions:

- **Track 1: Variant prioritisation** — genomic annotation, technical QC, inheritance-aware candidate assessment, transparent ranking and ACMG/AMP evidence review.
- **Track 2: Mechanism-informed drug repurposing** — analysis of BUB1B/BUBR1 variants, network-based biological interpretation and literature-derived pharmacological hypotheses.

## Interpretation and limitations

The BUB1B p.Leu737Ter variant is a predicted loss-of-function allele; p.Asn1002Lys is a missense VUS. The two heterozygous variants have unresolved phase. Proposed interventions, including premature-stop readthrough, NAD+-pathway support and 4-PBA, remain hypotheses: rescue efficacy and safety have not been experimentally established. BioGraph is a separate post-ranking interpretation service and does not modify prioritisation scores. The ACMG/AMP-style evidence map highlights criteria requiring expert review and is not a clinical classifier.

## Reproducibility and data safety

Only task-relevant, reviewed code, methods, configuration and tests should be committed. Patient VCF files, identifiable data, giant reference datasets, caches, virtual environments, backup versions and incidental local results are intentionally excluded. Reference resources and authorised example input datasets will be needed to reproduce full-genome runs. No genomic dataset is included in this initial repository documentation.

## AI assistance

Parts of software development and writing involved AI assistance. See each track's methods disclosure for provider, plan/tier and applicable data-handling settings.

## Source-code release status

Core research Python modules, preliminary tests, track-specific overview files and dependency requirements have been uploaded. The full original methods workbooks, reference resources, complete BioGraph source and end-to-end reproduction bundle are **not yet uploaded**. This is a partial research-code release, not yet a verified end-to-end reproduction package.

For local tests, create a Python 3.12 virtual environment, install `requirements.txt`, and run `python -m pytest tests/test_core.py tests/test_inheritance.py tests/test_acmg_preview.py`. Full application startup will require further source modules and external resources.

## Model 6 v1.2 — research results provenance

The public `mva_explorer/` directory contains a data-free snapshot of the Model 6 v1.2 software, not a fully verified end-to-end deployment. Confidential genomic input, private intermediates and large reference datasets are not included.

**Top 10 provenance:** The single ranking displayed above follows the saved ten-row [Model 6 v1.2 result CSV](submissions/track1/results/UscBen45_Track1_Model6_v1.2_saved_2026-10-22.csv) and is supported by the Results screenshot (Figure 4) in the researcher's later Track 1 report, `UscBen45_Track1_variant_prioritisation (report).pdf`. The report confirms BUB1B first, ADAMTS1 fifth and the BUB1B AR candidate-pair score of **0.882026** with phase **UNRESOLVED**. The original competition-portal CSV has not been recovered for a file-hash comparison, so this is a **report-supported saved result**, not an assertion of byte-for-byte portal-upload identity.

See [Model 6 v1.2 results provenance](docs/MODEL6_V1_2_PROVENANCE.md) for the supporting records and limitations.

### Track 2 literature and hypotheses

- [Recorded PubMed search strategy (working review notes)](docs/research/PUBMED_SEARCH_STRATEGY.md).
- [Preliminary scoping evidence register — pending completion](docs/research/SCOPING_REVIEW_PENDING_COMPLETION.md). The most populated saved workbook has 87 first-pass evidence entries from a 193-record worklist; **not** 87 fully verified completed studies. The [12-sheet Excel workbook](https://github.com/BIOMATICS-4/MVA-Genome-Explorer/blob/main/scoping%20review/G1_Continuous_Review_193_87_Studies_Audited.xlsx) is now published in the **scoping review** folder (reconstructed from the full worksheet text export).
- The four candidate research hypotheses in the most recent submission materials were **niacin, 4-phenylbutyrate (4-PBA), gentamicin and amlexanox**. **NMN was a mechanistic comparator only.** No patient-specific drug or rescue experiment was performed.

**Browse the preliminary extracted paper records on GitHub:** [87-study workbook text mirror — pending completion](docs/research/G1_CONTINUOUS_REVIEW_87_STUDIES_TEXT_MIRROR_PENDING.txt). The text mirror is a searchable companion to the [Excel workbook](https://github.com/BIOMATICS-4/MVA-Genome-Explorer/blob/main/scoping%20review/G1_Continuous_Review_193_87_Studies_Audited.xlsx); it does not preserve spreadsheet styling.

---

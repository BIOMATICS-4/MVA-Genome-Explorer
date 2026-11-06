# Track 2 — Drug repurposing investigation

**Objective:** Identify testable pharmacological hypotheses for BUB1B/BUBR1-related spindle-assembly checkpoint dysfunction, rather than assume a drug rescues a patient genotype.

**Variant mechanisms:** p.Leu737Ter may trigger nonsense-mediated decay or yield a truncated protein; p.Asn1002Lys is a missense VUS with unproven effects on protein folding or molecular interactions.

**Structural assessment:** Human BUBR1 UniProt O60566 and AlphaFold predicted structure; local residue inspection and solvent accessibility with Biopython; model preparation with PDBFixer; PROPKA side-chain pKa; OpenMM AMBER14 energy minimisation. For Asn1002, AlphaFold pLDDT was 91.06; the estimated solvent-accessible area was 31.50 Å² (WT Asn) versus 43.43 Å² (modelled Lys). Lys1002 predicted pKa was 9.28. No Asp/Glu side-chain atom was identified within 7 Å of the investigated residue in examined models. Raw minimised potential energies cannot establish folding ΔΔG.

**Candidate review:** Literature-derived hypotheses include NAD+-pathway support (NMN evidence does not prove niacin rescue), premature-stop readthrough, and 4-PBA as a conditional protein-stability intervention. DrugBank and PubChem support compound property and mechanism review; PubMed/PMC supports the preliminary literature assessment. Clinical and experimental rescue for the exact genotype remain unproven.

**BioGraph:** A separately developed network explorer integrated as a post-ranking biological interpretation adjunct. It does not feed information into genomic ranking.

**Limitations:** Variant phase, missense mechanism, off-target effects, cytotoxicity, pharmacological rescue and clinical efficacy are unresolved. Candidate selection is hypothesis-generating.

**Reproduction status:** Structural workflow files, complete reproducible source and the methods workbook still need independent file-level audit and upload. No patient data are included.

## Final nominated research hypotheses

The four selected candidates in the final Track 2 presentation were:

| Candidate | Research hypothesis | Evidence boundary |
|---|---|---|
| **Niacin** | Potential NAD⁺–SIRT2–BUBR1 abundance support | Indirect mechanism; the exact patient genotype has not been rescued |
| **4-phenylbutyrate (4-PBA)** | Conditional folding/proteostasis or residual-protein stability hypothesis | Evidence from other genetic diseases is not BUBR1-specific rescue |
| **Gentamicin** | Premature-stop readthrough research hypothesis | Exact stop context, tolerability and activity remain untested |
| **Amlexanox** | Readthrough and NMD-related exploratory hypothesis | Off-target action, feasibility and rescue are unproven |

**NMN was included only as a mechanistic comparator**, not an additional nominated candidate. These are literature- and mechanism-driven research hypotheses, not established therapies or validated efficacy for the proband.

### Scoping-review working records

The [recorded PubMed search strategy](../../docs/research/PUBMED_SEARCH_STRATEGY.md) covers **niacin, NMN and 4-PBA only**; it was not a search covering gentamicin or amlexanox. The [87-entry cumulative evidence register](../../docs/research/SCOPING_REVIEW_PENDING_COMPLETION.md) is a preliminary review pending completion and full-text verification. The [12-sheet Excel literature register](https://github.com/BIOMATICS-4/MVA-Genome-Explorer/blob/main/scoping%20review/G1_Continuous_Review_193_87_Studies_Audited.xlsx) is now available in the public **scoping review** directory as a reconstructed Excel-format copy of the complete text-extracted worksheets.


**Date and submission distinction:** The publicly documented search strategy comes from saved working records. The detailed Track 2 methods text does not specify a calendar search date, so no exact literature-search date is attributed to the submitted Track 2 methods here.

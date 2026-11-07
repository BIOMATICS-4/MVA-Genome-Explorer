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

### Scoping review and evidence sharing

The **current preliminary scoping review** includes [108 publications with examined Methods and Results](../../scoping%20review/README.md): 94 investigating 4-PBA, eight investigating NMN, and six investigating niacin. Its [Data extraction table](../../scoping%20review/Data_Extraction_Table.csv) includes the 108-publication synthesis and separately identifies other citations in the original worklist.

The [saved original 14-study preliminary extraction table](results/G1_Focused_Scoping_Review_14_studies.csv) is also available. This is a transcription of the preserved 14-row research file (eight 4-PBA-focused, four NMN-focused, and two niacin-focused studies); its exact byte identity with a hackathon portal attachment has **not** been independently authenticated. The 14 were used in the original Track 2 write-up, whereas the 108-study synthesis was developed subsequently.

The [recorded PubMed search strategy](../../docs/research/PUBMED_SEARCH_STRATEGY.md) covers **4-PBA, NMN and niacin only**. Gentamicin and amlexanox were separate Track 2 hypotheses and were not included in this saved search.

The earlier [87-entry register](../../docs/research/SCOPING_REVIEW_PENDING_COMPLETION.md), [12-sheet reconstructed workbook](../../scoping%20review/G1_Continuous_Review_193_87_Studies_Audited.xlsx), and [193-publication working bibliography](../../scoping%20review/Data_Extraction_Table.csv) remain available as historical research materials, **not** as the latest synthesised result. Of the original 193 worklist citations, 108 appear in the present synthesis and 85 do not. The older 18-excluded/67-pending classification has not been fully reconciled to individual records.

**Search-date and submission distinction:** The saved search has no independently verified calendar execution date in the submitted Track 2 methods text. The later 108-publication analysis must not be portrayed as part of the original hackathon submission.

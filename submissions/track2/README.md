# Track 2 — Drug repurposing investigation

**Objective:** Identify testable pharmacological hypotheses for BUB1B/BUBR1-related spindle-assembly checkpoint dysfunction, rather than assume a drug rescues a patient genotype.

**Variant mechanisms:** p.Leu737Ter may trigger nonsense-mediated decay or yield a truncated protein; p.Asn1002Lys is a missense VUS with unproven effects on protein folding or molecular interactions.

**Structural assessment:** Human BUBR1 UniProt O60566 and AlphaFold predicted structure; local residue inspection and solvent accessibility with Biopython; model preparation with PDBFixer; PROPKA side-chain pKa; OpenMM AMBER14 energy minimisation. For Asn1002, AlphaFold pLDDT was 91.06; the estimated solvent-accessible area was 31.50 Å² (WT Asn) versus 43.43 Å² (modelled Lys). Lys1002 predicted pKa was 9.28. No Asp/Glu side-chain atom was identified within 7 Å of the investigated residue in examined models. Raw minimised potential energies cannot establish folding ΔΔG.

**Candidate review:** Literature-derived hypotheses include NAD+-pathway support (NMN evidence does not prove niacin rescue), premature-stop readthrough, and 4-PBA as a conditional protein-stability intervention. DrugBank and PubChem support compound property and mechanism review; PubMed/PMC supports the preliminary literature assessment. Clinical and experimental rescue for the exact genotype remain unproven.

**BioGraph:** A separately developed network explorer integrated as a post-ranking biological interpretation adjunct. It does not feed information into genomic ranking.

**Limitations:** Variant phase, missense mechanism, off-target effects, cytotoxicity, pharmacological rescue and clinical efficacy are unresolved. Candidate selection is hypothesis-generating.

**Reproduction status:** Structural workflow files, complete reproducible source and the methods workbook still need independent file-level audit and upload. No patient data are included.

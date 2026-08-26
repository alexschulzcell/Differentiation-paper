> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# Preregistration S10 — R2 with a pinned primary reference (the C8 collapse as an exclusion)

Written and dated **2026-08-19**, **before** any new S10 number. It builds on
`PREREG_S9.md` and on the decision of 2026-08-19 to recompute the R2 impact
question with a pinned reference.

## 1. Scope and binding earlier results

- S8 and S9 remain binding. `GSE337700`, `GSE301458` and `GSE255460` were
  tested under the C8 reference; all three collapse to 93 to 99 %
  `transitioning`.
- **The C8 collapse is now a fixed exclusion criterion:** the six MSigDB C8
  signatures (`SU_HO_FOETAL_FEMUR_C0-C5`) are **not** admissible as a reference
  for adult skeletal tissue; a result with at least 90 % `transitioning` under
  C8 counts as a limit of the reference, not as a biological finding. This is
  an *exclusion*, not a retrospective reinterpretation of a result.

## 2. Two external, version-fixed references (named in advance)

Both from `celldex` **1.16.0** (pinned, with the package version logged), with
marker genes per label following a fixed, reproducible rule (below).

- **Primary reference R_A:** `celldex::HumanPrimaryCellAtlasData()`
  (`label.main`).
- **Second reference R_B:** `celldex::BlueprintEncodeData()` (`label.main`).

**The marker rule (identical for both):** per label, the genes with mean log
expression in the label ≥ 1.0 **and** a difference to the second-best label
≥ 1.0 log2. Labels without markers drop out. This rule was already implemented
before S10 and remains unchanged.

## 3. The assignment rule (identical to S7, S8 and S9)

Per cell and signature, the score is the share of signature genes with a count
above 0; the cell is assigned to the signature with the highest score if that
score is at least 0.10 and at least 2 percentage points above the second-best;
otherwise it is not assigned.

## 4. A coarse ontology for cross-reference stability (in advance)

| coarse class | R_A (HCA label.main) | R_B (Blueprint label.main) |
|---|---|---|
| chondral / osseous | Chondrocytes, Osteoblasts | Chondrocytes |
| endothelial | Endothelial_cells | Endothelial cells |
| stromal / mesenchymal | MSC, Fibroblasts, Tissue_stem_cells | Fibroblasts, Pericytes, Mesangial cells, Adipocytes |
| smooth muscle / muscular | Smooth_muscle_cells | Smooth muscle, Myocytes, Skeletal muscle |
| immune / haematopoietic | B_cell, T_cells, Macrophage, Monocyte, Neutrophils, NK_cell, DC, Myelocyte, Pro-Myelocyte, CMP, GMP, MEP, HSC_-G-CSF, HSC_CD34+, Erythroblast, BM, BM & Prog., Pro-B_cell_CD34+, Pre-B_cell_CD34-, Platelets | B-cells, CD4+ T-cells, CD8+ T-cells, DC, Eosinophils, Macrophages, Monocytes, Neutrophils, NK cells, Erythrocytes, HSC |
| other | Astrocyte, Neurons, Neuroepithelial_cell, Epithelial_cells, Keratinocytes, Hepatocytes, Embryonic_stem_cells, iPS_cells, Gametocytes | Astrocytes, Epithelial cells, Keratinocytes, Melanocytes, Neurons |

**Stability** is the share of cells whose coarse class agrees under both
references (counting only cells assigned under both). The threshold is 0.80.

## 5. Data and quality control

- `GSE255460` (osteoarthritic cartilage, 8 donors, against non-osteoarthritic
  control, 3 donors): the matrix already loaded,
  `GSE255460_sc_counts.txt` (symbol by cell) and `GSE255460_metadata.csv`
  (cell to `ID`, `trait` OA or control, `nFeature_RNA`).
- Quality control: `nFeature_RNA >= 1000` (from the authors' metadata; the
  matrix is already author-filtered; a mitochondrial fraction is not deposited
  in the supplement and is therefore not checked — documented, not fudged).
- The unit of comparison is the **donor** (OA1 to OA8 against C1 to C3), not
  the cell.
- `GSE301458` is **not** recomputed in S10: its control arm failed quality
  control (at most 71 cells per sample), independently of the reference.

## 6. The R2 statistic (sample or donor as the unit)

Per population that passes the gate (§7):
1. **Composition:** the share of the population per donor, OA against control
   (median, range, Wilson). That is the cell-type-shift statement.
2. **Within-population signal (exploratory):** the mean expression of the known
   173-gene module per population and donor, OA against control — as a
   descriptive, direction-free addition (no `dWT`, no `iv`, no new gene-set
   search).

## 7. Gate and stopping criteria

- **The gate:** at least TWO populations reach at least 10 % per donor and at
  least 100 cells across both groups, **and** the cross-reference stability is
  at least 0.80.
- **S10-AB1:** the gate is not met (fewer than 2 populations, or stability
  below 0.80). Report; no cell-type claim.
- **S10-AB2:** a planned analysis needs a threshold, gene list, population or
  reference to be relaxed after a result is known. **Not autonomously; stop the
  session.**
- **S10-AB3:** a robust finding contradicts a statement in the paper. No
  submission decision without a new decision by the authors.

## 8. Step E — the outcomes, in advance

1. **The gate is met and OA against control is robust** (at least one
   population with non-overlapping median Wilson intervals per donor): report
   the result as an exploratory cell-type-specific lesion response; the
   existing bulk null statement is not quietly retained; a supplement section
   and possibly a supplementary figure, no sixth main figure and no claim of a
   nineteenth data set.
2. **The gate is met but no population separates OA from control:** a clean
   null or not-decidable finding.
3. **The gate is not met (S10-AB1):** a limitation finding.

## 9. Technical

- Scripts carry a reference to this preregistration, the date, the role, the
  sources, the package versions, the seeds and the note "no `dWT`/`iv` main
  metric".
- `figure_style/publication_style.R` is not changed. The word "specific" is not
  used for the scissors.

## Addendum 1, 2026-08-19 — extending the data to `GSE337700`

**What was known before this addendum:** under the pinned HCA and Blueprint
references `GSE255460` reaches a cross-reference stability of 88 % but yields
only **one** population (`Chondrocytes`, 82 % of the cells passing quality
control) — cartilage is a single-cell-type tissue. The R2 gate (at least 2
populations) is therefore **not** attainable in GSE255460, independently of the
reference.

**Consequence (prospective, before any new GSE337700 number):** the primary R2
candidate is extended to **`GSE337700`**: fracture non-union against healed
fracture (3 against 3 samples, the authors' final set), a **heterogeneous**
callus tissue that shows several populations under the HCA reference (chondral
about 25 %, endothelial about 24 %, immune about 17 %, smooth muscle about
9 %, stromal about 24 % — from the S8 HCA analysis). The same rules of §2 to §7
apply to the non-union cohort (GSE337700) as well: HCA primary, Blueprint
second, quality control at least 1 000 genes and at most 20 % mitochondrial
reads, the sample as the unit, the gate at 2 populations and stability at least
0.80. The authors' exclusions
NUBF01 and NUBF03 remain visible as a sensitivity row.

## 10. The statement this phase starts from

> We are not changing the question, we are changing the reference. Now that the
> fetal C8 signature has demonstrably collapsed on adult skeletal tissue, we
> test with two pinned, independent cell-atlas references whether a lesion
> response (OA against control) is visible within stably annotated cartilage
> and skeletal populations. If it is not, that is a clean limitation result.

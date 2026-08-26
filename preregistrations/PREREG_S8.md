> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# Preregistration S8 — an orthogonal single-cell test of the claims of the paper

Written and dated **2026-08-19**, before any new S8 number and before any new
download. Changes to this document are made only as a dated addendum at the
end, stating what was known at the time of the change.

This session is **not** a rescue of the S7 extension. S7-AB2, the closed axis
question and `PREREG_S7.md` addenda 1 and 2 remain binding. `GSE196652` is not
attempted again as a replication. The paper stays submittable; the version with
five figures and iScience is the fixed fallback.

---

## 1. The three roles and the data sets

| role | the question put to the independent cells | data sets | fixed role |
|---|---|---|---|
| **R1** | Is the already known `dWT` programme visible as a module (173 genes) in independent mesenchymal or skeletal cell states along a differentiation? | `GSE166824` (BMMSC undifferentiated / D0, D3, D6), `GSE324998` (femoral osteolineage), `GSE255646` (four-donor multiome) | an external module statement; **no** new study convergence of 18 against 19 |
| **R2** | Is a lesion response visible in an independent clinical data set within predefined populations (non-union against healed)? | `GSE337700` (primary), `GSE241505`, `GSE150768` (fixed sensitivity contexts) | non-union against healed per population; **no** `iv`, **no** complete 2 x 2 |
| **R3** | Do cell shares and the within-population signal come apart in independent samples? (Bulk can deceive through cell mixture) | only data sets with at least 2 independent samples per group **and** 80 % external assignment stability | the share component against the within-population component |

These roles are fixed before the result. No candidate is exchanged because of
an interesting pattern.

## 2. The 173 genes: confirmatory, not discovering

`derived_data/reference_tables/S5_konvergente_gene.csv` is read in unchanged
(173 genes, known by name, with their direction from the main study). **No**
new gene discovery and **no** new cut-off is formed from S8 data. The 173 genes
are the **only** foreground in R1 and R3. Their use is reported as
confirmatory, not as discovering.

## 3. External, version-fixed cell-type signatures and the assignment rule

- **The first fixed reference:** the six MSigDB C8 signatures from
  `PREREG_S7.md` §3 (database `2026.1.Hs`, `msigdbr` 26.1.0, retrieved
  2026-08-19): `transitioning`, `resting_chondrocyte`,
  `hypertrophic_chondrocyte`, `proliferating_chondrocyte`, `fibroblast`,
  `endothelial`.
- **The second reference (sensitivity analysis only, named in advance):** the
  Human Cell Atlas cell-type markers, **pinned** through
  `celldex::HumanPrimaryCellAtlasData()` with the `celldex` version number and
  retrieval date, aggregated to skeletal and mesenchymal labels (MSC or
  fibroblast, chondrocyte, osteoblast, endothelial, immune). This reference is
  used **only** as a sensitivity analysis, to test the stability of the
  assignment (the R3 gate).
- **The assignment rule (identical to S7):** per cell and signature the score
  is the share of signature genes with a count above 0 (or expressed);
  assignment to the signature with the highest score if that score is at least
  0.10 and lies at least 2 percentage points above the second-best; otherwise
  unassigned.

**Important:** if the six MSigDB C8 signatures do not discriminate in
`GSE337700` because of the dominant `transitioning` assignment, a better-fitting
hand-made list is **not** built afterwards. The result of that transfer is then
simply "cell type not reliably determined" (eligible for S8-AB2). If the two
references do not agree, the result is likewise "cell type not reliably
determined".

## 4. Quality control, the sample replicate rule and the authors' exclusions

- **Quality-control parameters** (only for cells whose data set provides cell
  attributes at all): where cell attributes are available, the rule of S7
  applies: at least 1 000 detected genes, at most 20 % mitochondrial reads. For
  data sets that supply only an already filtered matrix or annotation, the
  authors' filtering is adopted and logged; only the minimum-unit check then
  applies.
- **The sample replicate rule:** the unit of comparison is **the sample**, not
  the cell. Cells of the same sample do not substitute for an independent
  biological unit. A role is load-bearing if at least **two independent
  biological sample units** are present per comparison group (R1:
  undifferentiated against differentiated; R2: non-union against healed). Fewer
  than two is S8-AB1.
- **The authors' exclusions:** in `GSE337700`, `NUBF01` and `NUBF03` remain
  visible as the authors' exclusions; the main table uses the authors' final
  set of six samples, the sensitivity row all eight. No exclusion is quietly
  discarded.

## 5. Primary and sensitivity analysis per role

### R1 — the known module (173 genes) along the differentiation
- **The primary statistic:** the share of the 173 genes whose expression rises
  in the predefined differentiation comparison (for example D0 or
  undifferentiated against D3 and D6 in `GSE166824`); missing or unmeasurable
  genes stay in the denominator. The direction per gene comes from the mean of
  the assigned cells per state. Every number carries the sample count, an
  uncertainty (a Wilson interval at sample level) and a noise expectation (a
  permuted direction count with the sample structure preserved).
- **R1 succeeds** if at least **two independent data sets** show the
  differentiation rise of the known set named in advance, and the direction is
  not carried by a single data set or by individual cells.
- **Sensitivity:** (a) a computation without the sample deviating most in the
  authors' assessment, (b) a Wilson interval rather than a point estimate,
  (c) no new cut-off.

### R2 — the lesion response within populations
- **Primary** is `GSE337700`. The result is non-union against healed per
  population; **no** `iv`, and **no** equating with a 2 x 2 interaction.
- The unit of comparison is the sample. Per population: median, range and a
  sample-level confidence interval (Wilson). "No difference" only with a
  predefined equivalence test and MDE80; otherwise "not decidable".
- If no population can be annotated robustly, it is **not** replaced by post
  hoc clusters.
- **The sensitivity contexts** `GSE241505` and `GSE150768` are used only if the
  species and arm limits are documented beforehand; they neither confirm nor
  refute the 173-gene or non-union statement.

### R3 — mixture against cell-intrinsic signal
- Only in data sets with at least 2 independent samples per group **and** at
  least 80 % assignment stability between the two external references (C8 and
  HCA).
- Population shares per sample; the 173-gene module per population and sample;
  optionally a bulk comparison that decomposes an observed module difference,
  using the linear mixture of the observed shares fixed in advance, into a
  share component and a within-population component.
- No covariate of the baseline, no new z null, no new gene-set search.
- If R3 is not load-bearing for lack of stable cell types, that is the result.
  The S7 signature bound is not relaxed.

## 6. MDE80 and confidence intervals, noise expectation, seeds

- Every quantitative number receives an MDE80 (the formula from `03_metric.R`,
  through the project function — **no** test statistic written by hand) **or**
  an explicit justification for reporting confidence intervals only. At
  n = 3 against 3, Wilson intervals are admissible; **no** biological null
  verdict is derived from them.
- Every direction count of the 173 genes is set against a predefined random
  gene and direction expectation; in permutations the sample structure is
  preserved (whole sample or data-set columns are flipped in sign, never
  individual cells permuted freely).
- Every search for cell types or markers receives a noise expectation or a
  reference comparison rate.
- Seeds depend on `(role, data set, hypothesis, round)`. Package versions and
  dates are logged in every script header.

## 7. Step E — outcomes fixed before the result

1. **R1 positive, R2 and R3 not decidable:** a supplement section as an
   orthogonal module validation; no sixth main figure and no claim of a
   nineteen-data-set replication.
2. **R2 shows a robust cell-type-specific lesion response:** the existing bulk
   null statement is not quietly retained; it is reported as an exploratory new
   narrative in need of independent confirmation.
3. **R1 negative or R3 not determinable:** the limitation is reported; the five
   figures stand.
4. **S8-AB4:** no submission decision without a new decision by the authors.

## 8. Stopping criteria

- **S8-AB1:** no new data set carries a role with at least two independent
  sample units per group. Write the report.
- **S8-AB2:** the cell-type annotation is not at least 80 % stable between the
  two external references. No cell-type claim.
- **S8-AB3:** the 173-gene module is not interpretable above the predefined
  MDE80 in the differentiation data. A negative or exploratory result, no new
  axis.
- **S8-AB4:** a robust external analysis contradicts a statement of the paper.
  Do not smooth it over and do not redefine; extend the protocol and name the
  consequence for the narrative.
- **S8-AB5:** a planned analysis needs a threshold, gene list, population or
  comparison group to be relaxed after a result is known. **Do not decide
  autonomously; stop the session.**

## 9. Technical rules

- New scripts carry in their header: the preregistration, the date, the role,
  the data sources, the package versions, the seeds and the sentence **"no
  `dWT`/`iv` main metric"**.
- `GSE337700` is not loaded again in full; the files already downloaded are
  used.
- `GSE166824`, `GSE324998`, `GSE255646`, `GSE241505` and `GSE150768` are loaded
  only after this document; archive sizes and SHA-256 are logged.
- Partial files are written per data set and round; existing partial files are
  skipped.
- RDS objects are loaded one at a time and memory is released after every
  sample.
- There is no second implementation of `kern()`.
- `figure_style/publication_style.R` is not changed.
- The word "specific" is not used for the scissors.

## 10. The statement this session starts from

> We are not looking for a new axis and not for a substitute for the missing
> replication. We test whether an already known differentiation module is
> visible in independent cells, and whether the clinical lesion context looks
> different within stably annotated populations from how it looks in bulk. If
> the data do not carry that, it is a clean limitation result and no reason to
> rewrite the narrative.

---

## Addendum 1, 2026-08-19 — the data files of two R1 data sets, named in advance

**What was known:** no S8 number yet and no new download beyond the GEO
supplement. `GSE324998` is downloaded as `GSE324998_RAW.tar`
(319 877 120 bytes); `GSE255646` as
`GSE255646_filtered_feature_bc_matrix.h5` (152 110 328 bytes) together with
`GSE255646_per_barcode_metadata.tsv.gz`; `GSE166824` as
`GSE166824_scRNAseq_expression_matrix_D0_D3_D6.txt.gz` (21 188 112 bytes).
For each of the three data sets the naming of the arms and time points is taken
from the existing metadata; **no** new arm or time classification is
introduced that is not documented beforehand.

---

## Addendum 2, 2026-08-19 — the R2 cell-type question in `GSE337700` (a fixed rule in advance)

For R2 in `GSE337700` the six C8 signatures remain the first reference. It is
known in advance that this transfer classifies about 98 to 99 % of cells as
`transitioning`. The following therefore holds firmly: the R2 population
analysis is reported as load-bearing **only** if at least two of the six C8
populations reach at least 10 % assignment per sample in the authors' final set
and at least 100 cells across both groups, **and** the second reference (HCA)
agrees with that assignment to at least 80 %. Otherwise the R2 result reads
"cell type not reliably determined — not decidable", and R3 does not apply in
this data set. That is a consequence of the assignment rule fixed in advance,
not a relaxation afterwards.

---

## Addendum 3, 2026-08-19 — the operational definition of the R1 statistic, and a structural finding about the R1 data sets

**What was known:** the R1 base data set `GSE166824` had been downloaded (with
its SHA-256 in the download log) and its column structure checked; **no**
direction or module number had been computed.

**The R1 direction statistic (defined firmly, before any number):** for
`GSE166824` the predefined differentiation comparison is `PS(D0)` or
undifferentiated against all `D3` and `D6` cells (2 kPa and 25 kPa, both with
and without induction medium, taken together as "differentiated"). Per gene the
direction is the sign of the difference in mean expression between
differentiated and undifferentiated. The statistic is the **share of the 173
genes with a measurable value** whose external direction reproduces the known
sign `ri` from `S5_konvergente_gene.csv` (ri = 1: a rise expected; ri = −1: a
fall expected); missing genes stay in the denominator. The noise expectation
comes from 500 permutation rounds in which the **condition-wide** sign pattern
is preserved and the assignment of differentiated and undifferentiated is
flipped jointly (genes are not permuted individually); MDE80 as in
`03_metric.R`. R1 "success" requires at least **two independent data sets** with
this differentiation rise named in advance.

**The structural finding (data preparation, not a direction number):**
- `GSE166824`: a single BMMSC preparation, with cells labelled
  `PS(D0)`, `2kPa(D3)`, `25kPa(D3)`, `2kPa(D6)` and `25kPa(D6)`; **no** donor
  or replicate identifier. There is therefore **no** second independent
  biological sample per comparison group (undifferentiated and differentiated
  both come from the same preparation). Under the sample replicate rule (§4)
  the cells are therefore **pseudoreplicates**, not biological units.
- `GSE324998`: 14 osteolineage 10x samples (barcodes, features and matrix per
  sample), **no** undifferentiated or temporal differentiation comparison in
  the data → descriptive module expression only, **no** direction verdict.
- `GSE255646`: one mixed multiome h5 (7 763 barcodes, 4 donors), **no**
  differentiation time comparison → descriptive module expression only, **no**
  direction verdict.

**The consequence (fixed before any direction number):** under the frozen
replicate rule, **none** of the three R1 data sets can on its own deliver a
documented "differentiation rise of the known set" with at least two
independent biological units per group (`GSE166824` is one preparation;
`GSE324998` and `GSE255646` carry no direction verdict). **S8-AB1 for R1 is
therefore satisfied in advance:** no new data set carries R1 with at least two
independent sample units per comparison group. The `GSE166824` direction is
nonetheless computed as a **single-preparation, exploratory** addition,
explicitly without the status of reaching the R1 success threshold.

---

## Addendum 4, 2026-08-19 — R2 and R3: structural conditions of validity

- `GSE337700`: 6 samples in the authors' final set (8 in all); three non-union
  against three healed. **No** time or perturbation axis. R2 is load-bearing
  only if at least two C8 populations can be annotated robustly (addendum 2).
- `GSE241505`: a clinical OPLL and ligament scRNA data set, **no**
  undifferentiated arm and no complete 2 x 2 cell → it serves **only** as a
  fixed sensitivity context for the R2 description and carries neither R1 nor
  R3.
- `GSE150768`: a mixed human and murine lesion model (mouse smooth-muscle
  cells, a `cond` label), **no** human undifferentiated-against-differentiated
  control pair → **only** a fixed sensitivity context with a documented species
  and arm limit: it can neither confirm nor refute the 173-gene direction or
  the non-union statement.
- R3 (mixture against cell-intrinsic signal) is load-bearing only in a data set
  with at least 2 independent samples per group **and** at least 80 %
  assignment stability between C8 and HCA. If that stability is missing
  (addendum 2), that is the R3 result.

> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# Preregistration S7 — the cell-type question and a possible nineteenth data set

Written and dated **2026-08-19**, before any S7 quality-control, aggregation or
metric number. It builds on `PREREG_whole_study.md` addendum 2 and on the
stopping rule of `PREREG_S6.md` §1, which remains in force. The confirmatory
axis question is not reopened.

Changes to this document are made only as a dated addendum at the end, stating
what was known at the time of the change.

---

## 1. What this session tests

`GSE196652` is a patient-paired single-cell data set of control bone and NF1
pseudarthrosis under osteogenic differentiation. S7 tests only:

1. whether the convergent differentiation programme runs within predefined cell
   populations or arises only from a mixture of populations;
2. whether the absence of convergence in `iv` can be explained by cell-type
   averaging in the bulk;
3. whether the data set can carry the existing computation as a nineteenth
   class M data set.

New axes, new gene-set scans and any reopening of the closed confirmatory axis
question are excluded.

## 2. Preparation and replicates

The complete aggregation step stands in `PREREG_whole_study.md` addendum 2. The
frozen rules in brief:

- raw counts; `Gene Symbol` to Ensembl through `org.Hs.eg.db` 3.20.0, unique
  assignments only, with ambiguities and unmappable symbols in the mapping CSV;
- at least 1 000 detected genes per cell and at most 20 % mitochondrial reads;
- one pseudobulk sample per `source x time point x patient_matched_1`, at least
  10 cells passing quality control, no artificial splitting and no merging of
  time points;
- `TP1` = `undiff`, `TP2`, `TP3` and `TP4` = `diff`;
- after that unchanged: `rlog`/VST, z standardisation, `kern()`, 20 classes,
  the VIF correction, `kontrast_f()` and `einzel_f()`.

S7-AB2 applies if, after quality control and aggregation, a cell of the 2 x 2
has fewer than two pseudobulk samples. S7-AB3 applies if fewer than two
populations reach the minimum size of §3. S7-AB4 applies at a
population-level detection limit above 1.0 z. S7-AB1 is the AB4 violation; with
S7-AB5 the reference implementation would be faulty.

## 3. Cell populations — in advance and external

No clusters are formed from these data and no population is named after its
result. The candidates are, in advance, the six external cell-type signatures
from **MSigDB C8 / Cell Type Signature, database version 2026.1.Hs**, loaded
with `msigdbr` 26.1.0 on 2026-08-19:

| population | MSigDB C8 `gs_name` |
|---|---|
| transitioning | `SU_HO_FOETAL_FEMUR_C0_TRANSITIONING_CELL` |
| resting_chondrocyte | `SU_HO_FOETAL_FEMUR_C1_RESTING_CHONDROCYTE` |
| hypertrophic_chondrocyte | `SU_HO_FOETAL_FEMUR_C2_HYPERTROPHIC_CHONDROCYTE` |
| proliferating_chondrocyte | `SU_HO_FOETAL_FEMUR_C3_PROLIFERATING_CHONDROCYTE` |
| fibroblast | `SU_HO_FOETAL_FEMUR_C4_FIBROBLAST` |
| endothelial | `SU_HO_FOETAL_FEMUR_C5_ENDOTHELIAL_CELL` |

For every cell passing quality control and every signature, the score is the
share of signature genes detected in that cell with a count above 0. A cell is
assigned to the signature with the highest score if that score is at least
**10 %** and lies at least **2 percentage points** above the second-best score.
Otherwise it stays unassigned. The signature genes are restricted to the genes
of the data set that map uniquely to Ensembl; the complete signature CSV with
name, version and retrieval date is written.

A population is evaluable for S7 if it has at least **10 assigned cells passing
quality control per source and exact time point**. The pseudobulk minimum size
of addendum 2 applies in addition. Fewer than two populations above this
minimum is S7-AB3.

## 4. Four computations in step C

All four computations run per evaluable population. With every reported number,
MDE80, the cell or pseudobulk count used, and a noise expectation are written.

### 4.1 `dWT` per population

`dWT` is computed within each population with the same reference function. The
173 genes from `derived_data/reference_tables/S5_konvergente_gene.csv` are known
by name. Their use is therefore **confirmatory, not discovering**, and is
reported exactly that way. The per-population sign agreement is computed as the
share of the 173 with a valid value; missing genes are reported in the
denominator.

### 4.2 Mixture against programme

The bulk `dWT` is computed twice: once from the real cells and once with the
population shares per source and time point held at the respective reference
value. The comparison quantity is the sign agreement of the 173 genes.

The alternative explanation "only a shift in mixture" counts as excluded if the
held shares preserve at least **80 %** of the real sign agreement of the 173
genes. That is the threshold fixed in advance, not a criterion adjusted
afterwards.

### 4.3 `iv` per population

The interaction term is computed per population with `kern()` and the unchanged
20-class main null. Consistent genes are counted at a sign consistency of at
least 90 %. The noise expectation comes from **1 000 rounds** in which whole
pseudobulk columns within the population are flipped in sign jointly, which
preserves the correlation between genes. Seeds depend on `(population,
hypothesis, round)`.

The alternative explanation "a cell-type-specific `iv` response is averaged
away in the bulk" counts as excluded if in **no** population the observed
number of consistent `iv` genes exceeds the respective noise expectation by
more than a factor of **two**, and §4.4 is at the same time below 1.0 z. If a
population exceeds that bound, the result is reported as an exploratory
cell-type-specific response and the existing narrative is not quietly retained.

### 4.4 Detection limit

The established procedure is applied unchanged at population level: neutral,
disjoint gene pairs from the external neutral set `GO:0007268` via
`org.Hs.eg.db` 3.20.0, the same 20-class null, additive offsets of
`0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00` z, five rounds per size, and
MDE80. The largest testable sides `15`, `30` and `60` are used, as far as both
sides have at least eight genes. The detection limit is the smallest offset at
which the median of the absolute corrected z reaches the threshold of 2. If it
lies above **1.0 z**, S7-AB4 applies, and a null result is then not reported as
excluding the alternative explanation.

## 5. Step D — a nineteenth data set

Only if step B passes the replicate and pool gates is `GSE196652` computed as
class `M`, arm `osteogenic`, through the same reference implementation. The
expected direction, in advance: the 18-data-set numbers should change only
slightly; a large jump would first have to be checked as an aggregation error.
A number seen before this preregistration would have to be marked in the
protocol as a procedural deviation; that has not happened.

## 6. Outcomes and journal — fixed before the result

- **Outcome 1:** §4.2 and §4.3 are satisfied, §4.4 lies below 1.0 z, and no
  population exceeds the noise expectation by more than a factor of two. The
  paper receives a sixth figure on cell-type resolution, and the main statement
  is extended by "and this holds within individual cell populations as well".
  The target journal is chosen afresh, with reasons, between **Genome Biology**,
  **Cell Reports** and **eLife**, with iScience as the fallback.
- **Outcome 2:** `iv` converges within a population. The paper is rebuilt; the
  new statement is exploratory and needs independent confirmation. The journal
  is chosen afresh.
- **Outcome 3:** S7-AB1, S7-AB2, S7-AB3 or S7-AB4 applies. The five finished
  figures stay unchanged; the single-cell attempt is reported in the supplement
  with the reason it failed, and the paper goes, as accepted in advance, to
  **iScience**.
- **Outcome 4:** S7-AB5. The reference implementation is not quietly repaired;
  the result is clarified technically and the journal decision is then taken
  afresh.

In every outcome: every number with its MDE80, every search with a noise
expectation, no matching to a covariate of the baseline, no new axis, and the
word "specific" not used for the scissors.

## 7. Release

This document is in force before the first S7 computation. The scripts carry
their individual decisions and seeds in their headers and write partial files
per source, time point, population, round and hypothesis; existing partial
files are skipped.

---

## Addendum 1, 2026-08-19 — a new scRNA GEO screen after S7-AB2

**What was known at the time of this addendum:** `GSE196652` has triggered
S7-AB2. No new GEO search has been run and no new data set has been assessed or
downloaded. This addendum is a new, prospective search step; it revokes neither
S7-AB2 nor the decision for iScience and changes no number of the existing
paper.

### Purpose and honesty

The new search is to establish whether the public literature holds a second,
independent single-cell cohort that can carry the cell-type question or the
bulk comparison. A hit is **not confirmation** before the same reference
computation has been run with a preparation fixed in advance. The search may
not be stopped at a data set with a convenient result, nor restricted to one
hit.

### The single new search query

NCBI GEO DataSets (`gds`), E-utilities `esearch`, run on 2026-08-19:

```text
("Homo sapiens"[Organism]) AND
("Expression profiling by high throughput sequencing"[DataSet Type]) AND
("single-cell"[All Fields] OR "scRNA"[All Fields] OR "scRNA-seq"[All Fields] OR "single cell RNA"[All Fields]) AND
(osteogenic[All Fields] OR chondrogenic[All Fields] OR osteoblast*[All Fields] OR chondrocyte*[All Fields] OR mesenchymal[All Fields] OR skeletal[All Fields] OR bone[All Fields]) AND
(knockout[All Fields] OR knockdown[All Fields] OR CRISPR[All Fields] OR mutant[All Fields] OR perturb*[All Fields] OR patient[All Fields] OR pseudarthrosis[All Fields])
```

All hits up to `retmax=1000` are screened, not only the first by relevance. The
complete hit list, the search URL, the date and the return code are stored. GSE
accessions already screened are marked `A8` and are not counted twice.

### Screening rules, in advance

A data set is a **candidate for download** only if all of the following are
documentable from the GEO metadata before the download:

1. a human single-cell RNA sequencing design with an accessible matrix or
   accessible raw counts and sample metadata;
2. a real control arm and a real perturbed or disease-associated arm;
3. an undifferentiated or day-0 starting point and a defined differentiation
   time point;
4. at least two independent biological units in every cell of the 2 x 2, with
   cells of the same sample not substituting for units;
5. no sample overlap with the 18 existing data sets and no design already
   excluded in the screening record.

A data set is not taken up even if the matrix is large, the perturbation sounds
biologically interesting, or a result could support the existing narrative.
Missing raw data, a missing undifferentiated arm, a missing control arm, only
one patient or donor, or unclear sample units are logged with the concrete
reason. `GSE196652` remains the demonstration that many single cells do not
create replication.

### The limit on downloading and analysis

After the complete screen, only candidates with a documented 2 x 2 and sample
replicates are downloaded. The download changes no existing figure. Before any
computation on a new candidate, a further dated addendum is written with the
aggregation step, the quality control, the cell-type definition and the
outcomes. The existing reference implementation stays unchanged; no new
single-cell metric is introduced.

## Addendum 2, 2026-08-19 — the orthogonal context data set `GSE337700`

The complete screen found no new data set with the 2 x 2 of undifferentiated
against differentiated and control against perturbation that `dWT` and `iv`
require. `GSE337700` is therefore not taken up as a main candidate. It is
nonetheless downloaded as an **orthogonal context data set**, because it is an
independent clinical scRNA cohort with three non-union and three healed
fracture controls and can therefore address the cell-type and lesion question,
though not the differentiation metric.

Fixed in advance:

- `GSE337700_RAW.tar` is downloaded with all eight deposited patient-level
  Seurat objects, including the two non-union samples the authors excluded
  because of batch interference;
- the authors' exclusions are not quietly discarded but reported separately in
  the metadata check;
- this data set does **not** enter the reference implementation as a nineteenth
  data set, delivers no `dWT` or `iv` number, and cannot confirm the bulk
  convergence statement;
- a possible context analysis is preregistered separately, descriptively and
  exploratorily, before the Seurat objects are read. It may describe only
  cell-type shares and the distribution of markers named in advance, and may
  construct no new axis and no new metric.

> Translated from the German original of 2026-08-18, with its addenda of
> 2026-08-18 and 2026-08-19. The content, the dates and every number are
> unchanged.

# Preregistration of the whole study — a typology of matrix production failure

Written on 2026-08-18, **before the first new data set was screened**. It builds
on §4.0 of the project brief. The two gate criteria are passed (S1
`PREREG_S1_FN1_induction_null.md`, S2 `PREREG_S2_positive_control.md`); this
preregistration governs everything that follows.

From here on: **changes only as a dated addendum at the end**, with a reason and
with a statement of which numbers were already known at the time of the change.
An addendum that relaxes a rule once a result is known is marked as such, and
the affected result counts as exploratory.

---

## 1. The statement being tested

> Matrix production fails along at least two transcriptionally distinguishable
> routes. The capacity-against-cargo axis separates them, it is measurable
> across systems, and membership of a quadrant says something about the
> underlying human genetics.

That is a statement about **structure in a distribution**, not about a single
data set. It can fail in three places: the points are not distributed in a
structured way (§5), the quadrants have no genetic signature (§6), or the metric
measures something generic (§5.4).

---

## 2. The metric — frozen

Unchanged from the reference script of the FN1 induction-null session. **That
file is from today the reference implementation.** Changes to it are addenda
under the rule above.

1. Count matrix → rlog/VST; where a transformed matrix is publicly available it
   is used and **not** transformed again. Ensembl IDs without a version.
2. z standardisation per gene over the samples of the respective contrast.
3. Per gene: the interaction `iv` = (KO_diff − KO_naive) − (WT_diff −
   WT_naive). The induction `dWT` = WT_diff − WT_naive.
4. The pool = genes with `dWT ≥ 0.5`. **A pool below 1000 genes means the data
   set drops out** (a rule from the generalisation script, already in the code).
5. Gene sets, external and unchanged: `S_MASCHINE` (GO:0006396, 0042254,
   0006412, 0005730, 0006888, 0048193, 0006486), `S_BIOSYN` (the first four),
   `S_DISTAL` (0006888, 0048193, 0006486, 0030968), `S_FRACHT` (GO:0005615).
   Overlaps removed on both sides, and `S_BIOSYN` and `S_DISTAL` kept disjoint.
6. The null: induction-matched, **20 classes** (no longer 4 — S1 showed that the
   refinement carries and costs nothing). The VIF correction as in the reference
   script.
7. **An MDE80 with every number.** A non-rejected null result without an MDE80
   is not reported.

### 2.1 Two metrics per data set

- **apparatus z** = the corrected z of the single set `S_MASCHINE`
- **cargo z** = the corrected z of the single set `S_FRACHT`

The plane of §4.3 of the brief is spanned by these two. The **decoupling**
(cargo minus machine) is reported in addition, but it is **not** the axis of the
typology — it projects the plane onto a line and cannot separate capacity
failure and programme failure from a global failure.

### 2.2 In addition, not as a criterion

`S_DISTAL` and `S_BIOSYN` separately; the blunting index; the T2
residualisation (`iv ~ dWT`, lowess) as a robustness row for **every** data set,
not only for FN1.

---

## 3. The data sets

### 3.1 Inclusion criteria, exhaustive

Mandatory, all four:
- an undifferentiated arm (day 0 or similar) **and** a differentiated arm,
- two groups (perturbation against control),
- a public count matrix or a public transformed matrix at gene level,
- n of at least 2 per cell of the 2 x 2.

After the computation, in addition: a pool of at least 1000 genes (§2.4).

### 3.2 Grounds for exclusion, exhaustively enumerated

`A1` no undifferentiated arm; `A2` no control arm; `A3` no gene-level matrix; `A4` n below
2 in at least one cell; `A5` "baseline" is a treatment control and not a time
control; `A6` pool below 1000; `A7` species not human; `A8` sample overlap with
an already included data set.

A data set may be excluded **on no other ground**. In particular, "the result is
implausible", "the effect is too small" and "it does not fit the pattern" are
**not** grounds for exclusion.

### 3.3 The screening procedure — the decisive protection

For every candidate the following is entered in the screening record **before**
the first computation on it: GSE, design, arm type, n per cell, perturbation
gene, **lesion class (§3.4)**, inclusion or exclusion with a code from §3.2.

Only once that entry stands is the metric computed. The order is binding. A data
set excluded after the computation can be excluded only with `A6`; any other
retrospective exclusion is an addendum under the rule at the head of this
document.

### 3.4 The lesion class — in advance and external

The perturbed gene is assigned to a class **before** the computation, from
external sources, not from the result:

- **`M` matrisome / secretory route:** core matrisome after Naba, or a member of
  `S_MASCHINE` or `S_DISTAL`, or a curated skeletal dysplasia gene (PanelApp 309
  green).
- **`N` other:** everything else (transcription factors without a matrix
  relation, lncRNA, signalling molecules, ion channels, chromatin factors).

Borderline cases are decided **before** the computation and documented with
their source. Class `N` is the negative class for §5.4.

### 3.5 The starting list

**Already computed, entering unchanged** (5 points): `LAMA5` KO osteogenic (M),
`LAMA5` KO chondrogenic (M), `FN1` C123R (M), `FN1` C231W (M), `SERPINA3` KD
(N), `MIR181A1HG` KD (N).
`COL2A1` patients do **not** enter (A1, no undifferentiated arm).
`GSE219011` plasma cell does **not** enter the plane (A2, no perturbation arm) —
it stays a calibration and is marked separately in the figure as a reference
pole.

**To be screened, in this order** (the 11 with priority from §3.3 of the brief):
`GSE218101`, `GSE221128`, `GSE227512`, `GSE196652`, `GSE319266`, `GSE113253`,
`GSE125167`, `GSE331389`, `GSE205432`, `GSE245585`, `GSE112318`.

**Of lower priority, each looked at once:** `GSE302312`, `GSE317531`,
`GSE226565`, `GSE202147`, `GSE198914`, `GSE166824`, `GSE272495`, `GSE37521`,
`GSE188759`, `GSE188760`, `GSE58123`, `GSE56900`, `GSE225446`.

The list is thereby closed. **New data sets may be taken up only if they are
found through the search query documented in §3.2 of the brief and the addendum
is written before the computation.** Searching deliberately for a data set that
fills a particular gap is excluded.

### 3.6 The minimum extent

The quadrant analysis (§5) is carried out only if **N of at least 8** data sets
meet the inclusion criteria, of which **at least 5 of class M** and **at least 3
of class N**. If that is missed, the study is reported as a **case comparison**,
not as a typology, and the statement of §1 is not made.

---

## 4. The arm rule

Osteogenic and chondrogenic are **never** computed against each other or
compared in their magnitudes (a factor of 4 in our own system). In the plane the
arm is carried as a mark. The quadrant assignment is made per data set in
absolute terms, not relative to another arm. Every statement of the form
"stronger than" is admissible only within one arm.

---

## 5. The quadrants

### 5.1 Boundaries

On each axis `|z_corr| < 2` counts as **flat**. That gives nine fields:

| | cargo ↓ (< −2) | cargo flat | cargo ↑ (> +2) |
|---|---|---|---|
| **apparatus ↑ (> +2)** | **programme failure** | **programme failure** | intact / expanding |
| **apparatus flat** | programme failure | **no lesion** | intact / expanding |
| **apparatus ↓ (< −2)** | *global failure* | **capacity failure** | **capacity failure** |

- **capacity failure** — the `LAMA5` quadrant.
- **programme failure** — the `FN1` quadrant.
- **no lesion** — where the negative controls must lie.
- **intact / expanding** — where the plasma cell lies as a reference pole.
- ***global failure*** (both sides collapse) is **not a type**. Such data sets
  are reported but **not** assigned to either of the two lesion types and do not
  enter §6. It is the most likely artefact cell for severely damaged systems and
  is neutralised in advance.

### 5.2 What a typology requires at a minimum

All three:
1. **at least 2 data sets in each of the two lesion quadrants**,
2. from **at least 2 independent studies per quadrant** (two mutations of the
   same GSE count as **one** study — `FN1` C123R and C231W are therefore **one**
   piece of evidence, not two),
3. **at least 1 data set per quadrant that does not come from our own
   laboratory.**

If one of these is missed: **the typology is not demonstrated.** What is
reported is a case comparison with an express statement that the quadrant claim
does not carry.

### 5.3 The structure test

The null hypothesis: the points scatter without structure around the origin.
Tested through the distribution of the angles in the plane (apparatus z,
cargo z), a Rayleigh or Hodges-Ajne test for uniformity, and in addition through
the share of points in the two named lesion quadrants against the share a
uniform distribution over the nine fields would lead one to expect.

**In advance:** this test is weak at N = 8 to 16. It is reported with its power
and is **not** a stopping criterion. What decides are §5.2 and §5.4.

### 5.4 The ubiquity test — this is where the paper dies, if it dies

If *all* systems show the scissors, the scissors are a generic differentiation
correlate.

**The rule:** among the data sets of **class N** (§3.4), the share with
`|decoupling z_corr| > 2` may be **at most 1/3**. If that is exceeded, the
scissors count as generic and **the study is not published as a typology.**

The counterpart: among class **M** the share with `|decoupling z_corr| > 2` must
be **at least 1/2**. If that is missed, the axis is not generally informative
for matrix defects; what is reported is then the subgroups in which it applies,
expressly as exploratory.

Both shares are reported with a confidence interval (Clopper-Pearson).

---

## 6. The human-genetics anchor

The procedure is unchanged from the distal-secretion disease-gene note:

- **PanelApp 309** (skeletal dysplasia), **green only**. Panel 1471 (short
  stature) is **not** used — it is a phenotype list, enriches only the cell-cycle
  contrast set (20/1009) and not distal secretion (2/523).
- **Comparisons between sets only.** Absolute enrichment is not interpretable,
  because panel genes are well-studied genes. The reference comparison stays
  **distal against biosynthetic** (39/523 against 35/2 192, OR 4.97,
  p 7 x 10⁻¹¹).
- **The cell cycle as a negative control** runs through the same test.
- **No enrichment test on the aggregated height loci.** 68.5 % of all genes hit
  a height locus.

**The new question preregistered here:** do the perturbation genes of the data
sets in the **capacity quadrant** lie closer to distal secretion than those in
the **programme quadrant**? The test: the share of perturbation genes per
quadrant that lie in `S_DISTAL` or in PanelApp 309 green, Fisher.

**Established in advance:** at N of at most 16 and two quadrants this test is
almost certainly underpowered. It is reported with its MDE, or minimum
detectable OR, and is **not** a criterion for §1. The carrying statement in
sentence 3 of §1 stays the set-against-set comparison already available (distal
against biosynthetic), which does not depend on N.

---

## 7. Confounders, to be computed per data set

For **every** included data set, not only in the compendium:

| confounder | check |
|---|---|
| library composition | repeat without the 500 most frequent genes |
| induction strength / ceiling | T2 residualisation (`iv ~ dWT`, lowess) |
| gene size | a length-matched null as in the cargo-load script |
| proliferation | repeat without the 1 724 cell-cycle genes |
| stress response | UPR and ERAD reported separately (direction, not only magnitude) |
| positive and negative control | in **every** data set, not only in the compendium |

A data set whose quadrant assignment flips under one of these checks is marked
**unstable** and does not enter §5.2 as evidence. It is nonetheless reported in
full.

---

## 8. Equivalence instead of a failed difference

Every statement of the form "X and Y do not differ" is carried as a **TOST**,
with the equivalence bound derived from the positive control of the same data
set. A non-rejected null result is not a finding and is not phrased as one.

---

## 9. Stopping criteria

The project ends, without rephrasing and without a substitute question:

- **AB1** — §5.4 missed: more than 1/3 of class `N` shows the scissors. They are
  generic.
- **AB2** — §5.2 missed: a lesion quadrant remains occupied by fewer than two
  independent studies. Then a case comparison instead of a typology; the target
  journals of §6 of the brief are thereby void, and plan B `Matrix Biology`
  remains possible as a case comparison.
- **AB3** — §3.6 missed: fewer than 8 eligible data sets.
- **AB4** — the reference implementation proves faulty and the already reported
  numbers of S1 and S2 change qualitatively.

---

## 10. What would refute the typology

Expressly, so that it is not negotiable afterwards:

1. The negative controls of class `N` show the scissors about as often as class
   `M` (§5.4). → A generic differentiation correlate.
2. The points accumulate in **one** quadrant or in the global-failure field
   instead of occupying two lesion quadrants. → One failure mode plus a severity
   axis, not a typology.
3. The quadrant assignment flips for more than half of the data sets under the
   checks of §7. → The axis measures a confounder.
4. `LAMA5` and `FN1` stay the only occupied points of their quadrants. → Two
   individual cases, not a type.

---

## 11. What is expressly not part of this study

From the record of fallen hypotheses, not to be tested again and not to be
mentioned except as a limitation: convergence in vitro against in vivo; "the two
arms are one superordinate programme"; "the scissors generalise to ECM genes";
a signal from outside (laminin → integrin → TOR); a secretory regulator arm as
the mechanism; comparisons across arms; `CDH1` without consulting the documented
trap; the word **"specific"** for the scissors.

---

## 12. The duty to report

Every screened data set appears in the manuscript or in the supplement — with
its inclusion or exclusion, the reason, and where included with all the numbers
of §2.1, §2.2 and §7, regardless of whether it fits the picture. A flow diagram
of the screen (candidates → checked → included, with exclusion codes) is part of
the main figures.

---

## 13. Release

**2026-08-18, before the first new data set was screened.** The two points at
which a milder version would have been possible are expressly confirmed in the
strict version:

- **§5.2** — two mutations of the same GSE count as **one** piece of evidence,
  and at least one study from another laboratory is required per quadrant. The
  consequence, deliberately accepted: at the time of this release the programme
  quadrant has **one** piece of evidence (`FN1`, `GSE251698`) and the capacity
  quadrant **none** from another hand. The screen under §3.5 must deliver in
  **both** quadrants, otherwise AB2 applies.
- **§5.4** — thresholds of N at most 1/3 and M at least 1/2, with a
  Clopper-Pearson interval.

This document is thereby in force. From here the addendum rule at the head
applies.

---

## Addendum 1, 2026-08-18 — exclusion code `A9` (single-cell data)

**The state of knowledge at the time of this addendum:** no number has been
computed from any of the affected data sets. All that is known is the assay
modality from the GEO metadata.

**The occasion — an error in §3.2.** It says there that the grounds for
exclusion are exhaustively enumerated. That is wrong: the list does not cover
the **assay modality**. During the screen it turned out that `GSE196652` and
`GSE319266` are single-cell data sets (time-series scRNA-seq and single-cell
CRISPRi/STING-seq respectively). Both meet the four criteria of §3.1 literally
if one aggregates them to pseudobulk — but the reference implementation of §2 is
defined on bulk group means and contains no aggregation step.

**The new code `A9`:** single-cell data from which a bulk 2 x 2 arises only
through a self-defined aggregation step that is not part of the frozen reference
implementation.

**Why that is not a relaxation:** `A9` **tightens** the inclusion criteria, it
does not relax them. No data set becomes eligible that was not eligible before.
And it is independent of any result: no value from `GSE196652` or `GSE319266` is
known.

**What it costs, and that is considerable.** `GSE196652` (somatic `NF1`
pseudarthrosis) would be class **M** under §3.4 — `NF1` is green in
PanelApp 309 — from another laboratory, osteogenic, with day 0 and a
patient-paired control. That is exactly the type of data set the capacity
quadrant needs under §5.2. The exclusion is methodologically grounded, but it
takes one of the study's most promising pieces of evidence away. **The
alternative — inclusion through pseudobulk as a documented extension of the
reference implementation — stays expressly open and would have to be written as
addendum 2 before anything is computed.**

**The decision of 2026-08-18, before any computation on the five new data
sets:** `A9` stands, and `GSE196652` is **not** taken up through pseudobulk. The
metric stays identical across all data sets; no second implementation is
introduced. The consequence is expressly accepted: if the capacity quadrant
fails for want of external evidence, AB2 applies and the work becomes a case
comparison. `GSE196652` is carried in the screening display of the manuscript
with the reason `A9`, so that it stays visible that it was checked and excluded
on methodological, not on substantive, grounds.

## Addendum 2, 2026-08-19 — a documented pseudobulk extension for `GSE196652`

**The state of knowledge at the time of this addendum:** the data set and its
GEO characteristics are known; no S7 quality control, no aggregation, no metric
and no cell-type figure has been computed or looked at. The decision to open the
route was expressly instructed by the user on 2026-08-19.

### The occasion and the extent

The alternative expressly left open in addendum 1 is opened for exactly one data
set: `GSE196652`. The extension is not a new gene-set test and not a new axis.
It answers only the cell-type question and tests whether the data set can be
prepared as a regular nineteenth point of the existing computation. `GSE319266`
and all other single-cell data sets stay excluded under `A9`.

### The aggregation step, fixed in full

1. The two public raw count matrices of `GSE196652` are used. The six annotation
   columns are not aggregated.
2. `Gene Symbol` is mapped to Ensembl without a version suffix using
   `org.Hs.eg.db` 3.20.0. A symbol is used only where there is exactly one
   Ensembl assignment; unmappable and ambiguous symbols are excluded and logged
   with their status, count and reason in a mapping CSV. After the assignment,
   duplicate Ensembl rows are summed.
3. A cell enters only if it has at least **1 000 detected genes** with a count
   above 0 and at most **20 % mitochondrial reads**. Mitochondrial reads are
   determined through the chromosome annotation `MT`/`M` and, where it is not
   labelled there, through symbols with the prefix `MT-`. There is no
   retrospective selection of cells by result.
4. The raw counts are grouped by **source x exact time point x donor**:
   `Control_bone` or `Pseudarthrosis`, `TP1` to `TP4`, and the donor documented
   in GEO as patient-matched, `patient_matched_1`. Time points are not merged.
   `TP1` is `undiff`; `TP2`, `TP3` and `TP4` are `diff`, following the existing
   project rule for several differentiation time points.
5. Each such group is exactly **one** pseudobulk sample. It must contain at
   least **10 cells passing quality control**; smaller groups are discarded. The
   counts are added per gene. Cells are not split randomly into subgroups: cells
   of the same patient are not independent biological replicates.
6. Only on these pseudobulk counts does the reference implementation continue
   unchanged: `rlog`/VST along the raw-count route, gene-wise z standardisation
   per contrast, `kern()`, the induction-matched 20-class main null, the VIF
   correction, `kontrast_f()` and `einzel_f()`. Aggregation is a preliminary
   step and not a second metric. A pool below 1 000 genes stays `A6`.

### The replicate rule and AB4

After quality control and aggregation, the 2 x 2 of `Control_bone` against
`Pseudarthrosis` and `undiff` against `diff` is checked for whether every cell
contains at least two pseudobulk samples. The number of cells is no substitute
for that condition. On a violation `S7-AB2` applies; the data set drops out and
step C is not interpreted as a null finding.

The **AB4 check** is a condition: **if an already reported number from S1, S2,
S5 or S6 changes qualitatively, the extension is faulty and is withdrawn.** The
existing 18 points are not re-aggregated and not altered; the check is
nonetheless documented against the frozen result files.

`A9` stays fully in force for all other data sets. What is opened is exactly
this one documented route, not the category and not the exclusion rule as such.

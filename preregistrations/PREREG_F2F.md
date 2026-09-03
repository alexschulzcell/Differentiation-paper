> Translated from the German original of 2026-08-24. The content, the dates
> and every number are unchanged.

# PREREG_F2F — the three-way decomposition as a hypothesis

**Written 2026-08-24**, before the decomposition was computed on any further
data set. This document does **not** make the earlier finding confirmatory —
GSE151315 was and remains unplanned and exploratory, and the text says so. It
makes the **replication** confirmatory.

---

## 1 · What already stands, and why it needs a rule fixed in advance

Panel **F2F** (from `06_orthogonal_layers/51_lineage_contrast.py`, run of 2026-08-23)
decomposes the failed calibration of the H3K27ac cohort (GSE151315) into two
steps:

| marker set | on the osteogenic axis | reading |
|---|---|---|
| undifferentiated | z −1.46 … −2.80 | the undifferentiated state **is left** |
| osteogenic | z −1.16 … −1.30 | the lineage is **not reached** |
| adipogenic | z +1.70 … +2.51 | the *wrong* markers rise |
| **module** | **z +3.10 … +7.02** | **runs regardless** |

The calibration passes 0 of 8; the module is above its own limit in 8 of 8.

Three reservations attach to this and remain in force independently of this
document:

1. **unplanned** — the rule fixed in advance for that run governed the lineage
   contrast, not this decomposition;
2. **calibration D underpowered** — MDE80 0.29–0.55 with two replicates per
   state; "failed" there means partly "not reached" and partly "not testable";
3. the module test is **licensed by no axis calibration**; the statement is
   the **pattern**, not the individual z value.

---

## 2 · The hypothesis, before looking at the data

> In a differentiation cohort with an undifferentiated arm,
>
> **(i)** the **undifferentiated markers fall** — z ≤ **−2** against their own
> null;
> **(ii)** the **lineage markers of the target axis do not rise** above
> z = **+2**;
> **(iii)** the **module lies above its own detection limit**.
>
> The replication is **confirmed** when **(i) and (iii)** hold and **(ii)**
> does not.
>
> If **(ii)** holds — the lineage is reached — and **(iii)** as well, that is
> **not a contradiction** but the **other case of decoupling** (as in the
> adipogenic axis of GSE332758); it is then reported as such.
>
> The decomposition is **refuted** when **(i)** holds and **(iii)** does not:
> the module would then run only where the lineage is reached.

The wording for the count: across the cohort the decomposition counts as
**replicated** when, among the data sets satisfying (i), the majority satisfy
(iii) and not (ii). It counts as **refuted** when, among the data sets
satisfying (i), the majority fail (iii). Every data set is reported
individually with its own detection limit, including those that do not satisfy
(i) — those are called **"not decomposable"** and carry **neither** a positive
**nor** a negative finding (project rule 1).

---

## 3 · The material on which it is tested

**All 18 perturbation data sets of the paper have an undifferentiated arm.**
The decomposition is computable there **without any new download**, from the
calibration data that the paper carries anyway. That is the material of this
preregistration.

That also names the level: the 18 are **transcriptome**, GSE151315 is
**chromatin**. This is deliberate — the decomposition needs no chromatin, it
needs an undifferentiated arm. If it reproduces in the transcriptome it is not
tied to one measurement technique; if it does not, that is a finding about the
measurement technique and is reported as one.

**What is not done:** the 18 data sets are **not 18 independent experiments**
in the sense of a meta-analysis — six study units of the donor level are among
them. No arithmetic mean over the 18 z values is therefore formed and no
pooled test statistic is computed over them. What is reported is a **count**
with a Wilson interval, and every data set stands individually in the table.

---

## 4 · The computation, fixed before the run

- **The statistic** is unchanged:
  `00_shared/_module.py::kontrast` — the mean of set A minus
  the mean of set B, against a null that draws two random sets of the same
  sizes from the same measurable background. Nothing is implemented anew.
- **The marker sets** come unchanged from
  `00_shared/_marker.py`, fixed and free of overlap. No set is
  touched for this run.
- **The axis** per data set is `dWT` (differentiated minus undifferentiated),
  unchanged from the frozen tables `20d_gene_*.csv`, the same source as in
  `03_lineage_calibration/10_calibration_18_datasets.py`. **Here `dWT` is a
  measurement axis and a tool, not a finding.**
- **Four quantities per data set**, each against **the same** null:
  1. the undifferentiated markers against all other marker sets,
  2. the lineage markers of the data set's **own** axis (osteogenic or
     chondrogenic) against all others,
  3. the lineage markers of the **other** axis against all others,
  4. the **module** — taken from `ws6_p1p2_modul_je_datensatz.csv`, the run
     that the paper carries anyway; **not** recomputed.
- **Seed** 20260821 from `_module.py`, 20 000 draws, unchanged.
- **Thresholds** as above: z ≤ −2 for (i), z > +2 for (ii), "above its own
  MDE80" for (iii). They stand here and are not moved after the run.

## 5 · What the result does to Figure 2

- **If the replication carries** (section 2), F2F shows **the decomposition
  across the 18 data sets**, coloured by calibration status, and GSE151315
  moves into the supplement as a single case (S9).
- **If it does not carry**, F2F stays as it is and the discussion says that the
  decomposition does not reproduce in the transcriptome data.

## 6 · The expectation, recorded so that it is not reinvented afterwards

Under the model we expect the **failed** data sets to show the same pattern as
GSE151315 (the undifferentiated state left, the lineage not reached, the module
running), and the **passing** ones to show rising lineage markers in addition.
This expectation is explicitly **not** the decision rule — that stands in
section 2 — but is recorded so that any deviation from it stays visible.

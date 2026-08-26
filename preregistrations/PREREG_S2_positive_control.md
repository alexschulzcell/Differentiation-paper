> Translated from the German original of 2026-08-18. The content, the dates
> and every number are unchanged.

# Preregistration S2 — plasma-cell differentiation as a positive control

Written 2026-08-18, **before the first computation**. A calibration criterion
of the project. Deviations are added as a dated addendum at the end.

## 1. The question

The metric claims to be able to measure the build-up of the secretory
machinery. The textbook case of that build-up is the B cell becoming a plasma
cell: the differentiation that multiplies the ER compartment. **If the metric
shows no massive, coordinated rise of the apparatus there, it is not calibrated
and everything that follows is worthless.**

## 2. The data set

`GSE219011` (Tooze laboratory, Leeds). 22 RNA-seq samples, human peripheral
B cells, in vitro differentiation, donors A, B and C:

| time point | cell type | donors |
|---|---|---|
| day 0 | B cell | A, B, C |
| day 0 + 30 min … + 12 h, day 1 | B cell | A, B |
| day 3 | plasmablast | A, B, C |
| day 6 | plasmablast | A, B, C |
| day 10 | plasma cell | A, C |
| day 13 | plasma cell | A, C |

The public VST matrix with Ensembl identifiers
(`GSE219011_RSEM_Gene_count_VST_wGeneSymbols.txt.gz`, 28 370 genes). It is
already transformed and is therefore, like the FN1 matrices, **not**
transformed again.

**Primary contrast (fixed in advance):** day 0 against **day 6**, n = 3 against
3, all three donors in both arms. Chosen because the pairing is complete.

**Confirmatory contrast:** day 0 against **day 13**, n = 2 against 2 (donors A
and C). The more mature state, but with less power.

The intermediate time points (30 min to day 3, and day 10) enter **no** test.
They are carried only as a time course.

## 3. What is measured differently here from S1, and why

In S1 the target quantity is the **interaction** (genotype x condition). Here
there is no perturbation, only one arm. What is tested is therefore the
**induction term itself**: the rise from day 0 to day 6 per gene, averaged over
the gene set, against a matched null.

The matching is on the **baseline expression at day 0** (decile classes of the
mean VST value), not on the induction — that is the target quantity here. A
global shift effect ("everything rises in the plasma cell") cancels out,
because the null is drawn from the same pool with the same baseline expression.

The VIF correction follows the project standard; the mean correlation is
estimated on the residuals of **all 22 samples** against the time-point factor,
because the primary contrast alone would have only 2 degrees of freedom. MDE80
is reported with every null.

The gene sets are unchanged: `S_MASCHINE`, `S_DISTAL`, `S_BIOSYN`, `S_FRACHT0`.

## 4. Control sets

- **Neutral set (negative control within the same data set):** GO:0007268,
  chemical synaptic transmission. A programme that has no business in B cells.
  If it does not stay flat, the null is defective and the test says nothing.
- **Cargo (GO:0005615)** is reported alongside but is **not** a negative
  control here: in the plasma cell the cargo (immunoglobulin) rises as
  expected. Its value is a description, not a criterion.

## 5. Decision rule

**CALIBRATED** — only if all four hold:
1. primary contrast: `S_MASCHINE` z_corrected > **+3**;
2. primary contrast: `S_DISTAL` z_corrected > **+2** (the skeletally relevant
   half has to move too, not only the ribosomal one);
3. primary contrast: neutral set |z_corrected| < 2;
4. confirmatory contrast at day 13: the same sign for `S_MASCHINE` and
   `S_DISTAL`.

The threshold of +3 rather than +2 is deliberately higher than in S1: what is
required is not significance but a **massive** rise in the system that
represents the textbook case. A barely significant machinery in the plasma cell
would be a warning sign, not a pass.

**NOT CALIBRATED** — if (1) or (2) is missed.
→ The metric does not measure what it claims to measure. **The project ends
here**, independently of the outcome of S1.

**NULL DEFECTIVE** — if (3) is missed. No verdict on the metric; the matching
is corrected and the test repeated. This is documented explicitly and is
**not** read as a pass.

**WEAK CONFIRMATION** — (1) to (3) met, (4) not. The calibration then counts as
passed, but the contradiction between day 6 and day 13 is reported in the
manuscript and not left out.

## 6. What this control does NOT show

- It calibrates the measurement of the **apparatus**, not that of the
  **decoupling**. A system without a perturbation arm can say nothing about the
  interaction quantity.
- B cell to plasma cell is a change of cell type, not a differentiation arm of
  a mesenchymal lineage. The control establishes the sensitivity of the metric
  to secretory capacity, not its transferability to the skeletal system.
- n = 3 and n = 2, with three and two donors. A non-significant secondary
  result in this data set is not interpretable without MDE80.

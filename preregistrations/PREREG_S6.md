> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# Preregistration S6 — discovery and validation, in separated halves

Written and dated **2026-08-19**, **before the first metric of the scan**. It
replaces `PREREG_S5.md` as the governing session preregistration; `PREREG_S5.md`
and `PREREG_whole_study.md` remain in force unchanged, except where this
document says otherwise.

It builds on the S6 brief. Changes to this document only as a dated addendum at
the end, stating which numbers were known at the time of the change.

---

## 0. What was known before the date of this document

Disclosed on the same pattern as `PREREG_S5.md` §0. Three things were computed
**before** this dating, and all three are structural quantities, not results:

**0.1 The split itself** (the split script). It uses only arm, lesion class,
null integrity, study, perturbed gene and pool size. **Not one row of the
per-point contrast tables was read by that script.** The split stands by name in
§2 and is not touched again thereafter.

**0.2 The set check** (the set-check script). What was collected is set sizes
and pool intersections at the **nine discovery points**, plus a pure timing
measurement. From it comes the number **K = 105**. That is the same procedure as
the set check before `PREREG_S5.md`, and it is necessary: K has to stand in
advance (§4.3), and K can only be determined once one has counted how many
candidates are testable at all. The timing measurement computed three contrasts
at two discovery points; their `z_corr` values were **not written out and not
looked at** — what was measured was the number of seconds.

**0.3 What is known from S5 and therefore cannot be a discovery.** The contrast
values of the three S5 hypotheses lie open for **all eighteen** points in the
per-point contrast tables, and the three conspicuous H2 points are known by name
(`LAMA5` osteogenic +6.36, `LINC01638` KD +4.31, `RB1` del isogenic +3.99). From
that follows the **honesty rule** in §4.5 (c): the GO terms from which
`S_STRESS`, `S_ABBAU` and `S_ECM_KERN` are built are **removed** from the
candidate list, and a candidate that overlaps strongly with one of the three
sets is not a discovery and does not go into the validation.

**0.4 What is not known.** No contrast of a candidate set against the reference
side has been computed or looked at at any point. The **validation half has not
been opened since the end of S5** and will not be opened until step D.

---

## 1. The stopping rule — word for word from §2.1 of the brief

> **A binding stopping rule: if no hypothesis holds in S6 either, the axis
> question is closed.** It is not reopened afterwards with new gene sets, with
> more points or in a different construction. What is then written is the
> methods and negative-finding paper, which already stands in the S5 report §6
> with a structure, a target journal and a figure plan.

It is not weakened in this document, not made conditional and not altered by an
addendum.

---

## 2. The split — by name, and frozen from here

Drawn with seed **20260819** from the **98** splits that satisfy all the
stratification conditions. Blocks are formed **by perturbed gene** (13 blocks),
so that no point in the validation has a sibling in the discovery — which is
stricter than §5.1 of the brief requires and closes the leakiest place in the
separation. The conditions, hard and in advance: exactly 9 points per half,
exactly 3 chondrogenic, 3 or 4 of class `M`, 4 to 6 with a demonstrably intact
null. It was **not optimised** but drawn from the admissible splits.

### 2.1 The discovery half (step C)

| # | point | arm | class | null | pool |
|---|---|---|---|---|---|
| 3 | `FN1` C123R | ch | M | defective | 11 461 |
| 4 | `FN1` C231W | ch | M | defective | 11 815 |
| 5 | `SERPINA3` KD chondrogenic | ch | N | ok | 1 347 |
| 11 | `RB1` +/− patient line | os | N | defective | 8 369 |
| 13 | `SERPINA3` KD osteogenic | os | N | ok | 1 033 |
| 14 | `ERCC6L2` KD | os | N | ok | 6 502 |
| 15 | `TP53` LFS | os | M | ok | 6 219 |
| 17 | `RB1` mut isogenic | os | N | defective | 5 720 |
| 18 | `RB1` del isogenic | os | N | defective | 5 974 |

n 9 · chondrogenic 3 · class M 3 · null intact 4

### 2.2 The validation half (step D) — **untouched until step D**

| # | point | arm | class | null | pool |
|---|---|---|---|---|---|
| 1 | `LAMA5` KO chondrogenic | ch | M | ok | 9 861 |
| 2 | `LAMA5` KO osteogenic | os | M | ok | 9 063 |
| 6 | `MIR181A1HG` KD | os | N | ok | 12 039 |
| 7 | `LINC01638` KD | os | N | ok | 3 941 |
| 8 | `ARSB` MPS VI | ch | M | ok | 7 505 |
| 9 | `ACVR1` FOP | ch | M | ok | 9 132 |
| 10 | `RNF4` KD | os | N | defective | 8 306 |
| 12 | pseudarthrosis | os | N | defective | 2 251 |
| 16 | `YAP/TAZ` KD | os | N | defective | 2 712 |

n 9 · chondrogenic 3 · class M 4 · null intact 6

**No look, no interim analysis and no figure at these nine points until
step D.** If one is opened, that belongs in the protocol and the point counts as
a discovery point.

### 2.3 What §14 of the brief demands here, recorded

The three post hoc points behind the H2 conspicuity are distributed as follows:
**`RB1` del isogenic (18) lies in the discovery**, **`LAMA5` osteogenic (2) and
`LINC01638` (7) lie in the validation.** That is the draw, and it is **not
corrected**. Two of the three therefore lie in the half that under §14 they
expressly may not determine — which is favourable for the cleanliness of the
session, but was neither sought nor steered.

---

## 3. Step B — the door stays shut, n = 18

§6 of the brief leaves it open whether to open the pseudobulk extension as
**addendum 2** to the whole-study preregistration or to leave it. **It is not
opened.** The full reasoning stands in the pseudobulk session protocol; the
three carrying grounds:

1. An addendum 2 would be a **relaxation** of an existing preregistration (`A9`
   would fall away for a class of data sets). §12.2 of the brief requires for
   exactly this case that the autonomous session **halts** and does not decide
   itself. The session therefore takes the other route, expressly named as
   admissible in §6.
2. **n was not the bottleneck in S5** (§1 of the brief): the jump from 11 to 18
   points changed not a single verdict. Two to four further points would have
   shifted the threshold of the validation half from 9/9 to 10/10 or 11/11
   without improving discrimination.
3. The binding bottleneck of this session is the **sensitivity of the metric**
   (§6), not n. The computational cost of an scRNA aggregation including the AB4
   check would come at the expense of exactly that step.

**n = 18** is thereby the binding number for §2. The door stays open for a later
session, unchanged as described in `PREREG_whole_study.md` addendum 1.

---

## 4. Step C — the scan

### 4.1 Construction

The metric is unchanged: **a contrast of two gene sets against the 20-class main
null**, `kontrast_f()` from the reference implementation, called unchanged. No
second implementation. The VIF correction always, `NVIF` 100, the seed keyed to
(point, candidate, round).

**The reference side `REF` = `S_MASCHINE`**, unchanged from
`PREREG_whole_study.md` §2.5. It therefore stands **in advance** and is not
chosen out of the result. The scan asks a single question, and it is the axis
question of this work: *which external gene set is systematically decoupled from
the biosynthetic capacity scaffold?*

Per candidate `C`: side A = `C ∖ REF`, side B = `REF ∖ C`. Both sides are
thereby disjoint, as in all three S5 hypotheses.

**The direction: two-sided.** A scan has no prediction. The direction is fixed
only in step C, from the majority sign of the discovery half, and then enters
step D **frozen**.

### 4.2 `NB` is **not** reduced

§7.1 of the brief permits `NB` 500 in the scan. That permission is **not taken
up**: at K = 105 instead of the 1 000 assumed there, the scan is computable with
**`NB` = 2000** (measured: 0.65 s per contrast at the smallest and 2.5 s at the
largest point at `NB` 500; at `NB` 2000 roughly 3.5 times that). The
recomputation of §7.1 thereby falls away entirely, and discovery and validation
compute with **the same** number.

**The one exception, expressly:** the empirical noise arm (§4.4) computes with
`NB` = 500. It produces no reported metric but a distribution, and it is not
computable at `NB` 2000. That is the only place in this session where `NB` is
reduced, and it is named here.

### 4.3 The scan extent K — in advance, as a number

**K = 105.**

How it arises, entirely rule-driven (the set-check script):

| step | number |
|---|---|
| GO BP terms with a gene assignment (`org.Hs.eg.db` 3.20.0, `GO2ALLEGS`) | 15 413 |
| of these within the size window of **50 to 500 genes** genome-wide | 2 564 |
| less the **16 GO terms** from which the sets of this work are built | 2 557 |
| of these **testable at all nine discovery points** (at least 15 genes per side in the pool) | **105** |

The size window is set in advance and with reason: below 50 genes the pool
intersection at the smallest point (pool 1 033) regularly falls under the
testability bound, and above 500 genes the sets become generic — exactly the
objection on which `S_ABBAU` (891 genes) hung in S5. The testability bound of 15
genes per side is stricter than the 8 that `kontrast_f` requires, so that a hit
does not hang on a handful of genes.

**K is therefore not chosen but counted.** It stands before the first metric and
is not extended, not even if the scan delivers nothing.

### 4.4 The noise expectation — both routes, mandatory

**(a) Analytically.** An exact sign test, two-sided, n = 9:

| hit threshold | p | expected hits at K = 105 |
|---|---|---|
| **9 / 9** | 0.0039 | **0.41** |
| at least 8 / 9 | 0.0391 | **4.10** |

**(b) Empirically.** A permutation of the group assignment: at each of the nine
discovery points the `genotype` labels are exchanged **within** the respective
`condition` (seed 20260819 plus point plus round), `kern()` is recomputed and
the whole scan is repeated over all K candidates. **`NPERM` = 10 rounds per
point**, `NB` 500.

From that, **1 000 noise rounds** are formed by drawing, per round and per
point, **independently** one of the 10 permutations (10⁹ combinations, seed
20260819). That is admissible because under the null hypothesis the points are
independent, and it preserves the correlation structure **between** the 105
candidates within one point — exactly the quantity the analytical computation
does not know. From this distribution are reported: the mean number of
candidates with 9/9 and with at least 8/9, the distribution of the **best** sign
count, and the empirical p of the actual best hit.

> **A scan hit is not reported without its noise expectation.** If the observed
> number of hits matches the expected one, that is the result of step C,
> **S6-AB2** applies, and the session proceeds to outcome 2.

### 4.5 The control on the baseline expression (§7.3) — and what it excludes

Under §2.4 of the brief, a free scan on this scale sorts by the baseline
expression. Therefore, for **each** of the 105 candidates and at every discovery
point:

- **ΔBaseline** = the mean baseline expression (`basis` from `kern()`, the
  starting point of the WT undifferentiated arm on the gene-wise z scale) of side A
  **minus** that of side B, in each case over the genes in the pool.
- **ΔLevel** = the same for the absolute mean expression `expr` before the z
  standardisation, standardised over the pool.

Both are reported per candidate as the median over the nine points and plotted
in the figure against the hit list. In addition the Spearman correlation between
a candidate's median `z_corr` and its median ΔBaseline is reported over all 105
candidates — it quantifies how strongly the scan actually sorts by the
confounder.

**The exclusion rule, deterministic and in advance:** a candidate is an
**artefact candidate** and does **not** go into the validation if
`|median ΔBaseline|` lies in the **top third** of all 105 candidates. The bound
is thereby defined structurally (a tertile of the candidate distribution) and
not read off a result. Artefact candidates are reported in full, with their
number and rank.

**Not admissible, and it is not done:** matching the null to the baseline
expression. That is shown in the null-diagnostic protocol §4 to be structurally
inadmissible and stands in the record of fallen hypotheses. The baseline
expression is **reported**, in order to exclude candidates, not **computed
away**.

### 4.6 How the scan yields **exactly one** hypothesis

The rule stands here in full and is applied as it stands.

**(a) The ranking criterion, primary:** the **sign count** `v` = max(number of
points with `z_corr` > 0, number of points with `z_corr` < 0) over the nine
discovery points, descending.

**(b) The minimum threshold:** `v` at least 8. Candidates with `v` of at most 7
are not considered (p at least 0.18, which is noise).

**(c) Exclusions, in this order:**
1. **Artefact candidates** under §4.5.
2. **Already known S5 hypotheses.** The 16 constituting GO terms are already
   removed from the candidate list. Additionally excluded is any candidate with
   a **Jaccard similarity of at least 0.5** to `S_STRESS`, `S_ABBAU`,
   `S_ECM_KERN`, `S_DISTAL`, `S_BIOSYN` or `S_FRACHT0`. It is then not a
   discovery but a renaming.
3. **Too few points with an intact null.** A candidate whose sign count `v` is
   carried solely by the five points with a defective null (that is, at most 2
   of the 4 null-intact points in the majority direction) is excluded. That is
   the counterpart of version (c) of `PREREG_S5.md` §5.4, set here in advance as
   an exclusion rather than as a secondary count.

**(d) Tiebreak, deterministic, in this order:**
1. the larger `|median z_corr|` over the nine points;
2. the smaller `|median ΔBaseline|`;
3. the smaller numerical GO identifier.

**No selection by the plausibility of the story.** The order above is complete;
there is no step at which a human decides which narrative sounds better. That is
the rule on which "identity holds, output flips" died, and it is not reopened
here.

**(e) The result:** exactly **one** candidate, with a **frozen direction** (the
majority sign of the discovery half). If no candidate remains after (b) and (c),
**S6-AB2** or **S6-AB3** applies, and the session proceeds to outcome 2 without
step D.

### 4.7 The sensitivity — the detection limit of the metric (§7.4)

Independent of the outcome of the axis search, and belonging in **both**
versions of the paper. Fixed in advance:

At each of the nine discovery points two **disjoint** subsets of `S_NEUTRAL` in
the pool are drawn, at sizes of **15, 30 and 60 genes per side** (capped at half
of `S_NEUTRAL` in the pool). A known offset **δ** is added to `iv` on the genes
of side A, with **δ ∈ {0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0}** in z units.
**5 rounds** per combination, seeds keyed to (point, size, δ, round). `NB` 2000.

**The detection limit** is the smallest δ at which the median `|z_corr|` over
the rounds reaches **at least 2** — the project's flatness threshold. It is
reported per set size and per point, and additionally as the median over the
points.

**What that answers:** whether the negative findings of S4 and S5 mean "no
effect" or "too small for this metric". The answer is reported **however it
turns out** — in particular also if it qualifies the negative findings so far.

**What it does not answer,** and this is not claimed: δ is an additive offset on
`iv` with the noise otherwise unchanged. A biological effect that also changes
the spread is not captured by it.

### 4.8 The candidate from §14 of the brief — it does **not** run along, and why

§14 of the brief leaves it open to enter a reworked K6 ("PanelApp 309 against a
panel of another organ group of similar extent") in advance and run it alongside
the scan. **It is not entered.** The reason is technical, not substantive:

- Exactly **two** PanelApp panels are pinned in this work: `panel_309.json`
  (skeletal dysplasia, v10.3, 440 green genes) and `panel_1471.json` (short
  stature).
- **Panel 1471 is excluded as a comparison side**, by `PREREG_whole_study.md`
  §6: it is a phenotype list that enriches the cell-cycle contrast set (20/1009)
  and not distal secretion (2/523). It is not "another organ group" but the same
  phenotype area.
- A **third, freshly retrieved panel** would be a comparison side chosen by
  hand. Which organ group "of similar extent" one takes is exactly the freedom
  this session's split is built against, and §11 of the brief requires gene sets
  to be external and version-fixed — not external and **picked**.

The decision is a cost-benefit judgement and is named as one: it costs the last
opportunity to test the human-genetics anchor **within** the transcriptome data,
because §1 closes the axis question afterwards. It is nonetheless taken this
way, because a hand-picked comparison side would not carry the finding it could
produce. The anchor stays what it has been since S5: a finding **outside** the
transcriptome data, reported as a boundary condition.

---

## 5. Step D — the validation

### 5.1 What is computed

**Exactly one** hypothesis — the one that emerged from §4.6, with the direction
frozen there — at the **nine** points of §2.2, with **`NB` 2000**, `NVIF` 100,
the 20-class main null, `kontrast_f()` unchanged.

In addition, unchanged from `PREREG_S5.md` §5.1:
- the four §7 checks (a length-matched null, the T2 residualisation, without the
  500 most frequent genes, without the cell-cycle genes),
- the set-size check (20 subsets of the larger side at the size of the smaller
  one, seeds keyed to (point, round)),
- the two single sets descriptively,
- **an MDE80 with every number**,
- the **neutral contrast**: the reference side `S_MASCHINE` lies throughout
  above 105 genes in the pool and therefore outside the range checked in the
  convergence-axes session protocol §4, so it is computed — 10 rounds per point.
  **The refinement as in S5 (decision 3):** the constructible size class is that
  of the **smaller** side, capped at half of `S_NEUTRAL` in the pool;
  `S_NEUTRAL` has 614 genes and cannot represent the size of the reference side.
  The size actually used stands in every row of the CSV.

### 5.2 The threshold — exact, with its p value

A paired sign test, an exact binomial test, two-sided, p = 0.5. **α = 0.05, no
multiplicity correction** — that is the whole point of the split: the
denominator holds **one** hypothesis.

| testable points | required | p |
|---|---|---|
| **9** | **at least 8 of 9** in the frozen direction | 8/9 → 0.0391; 9/9 → 0.0039 |
| 8 | **8 of 8** | 0.0078 (7/8 → 0.070, not enough) |
| 7 | **7 of 7** | 0.0156 |
| at most 6 | **not validatable** → S6-AB4, outcome 2 | — |

Testability means: at least 8 genes per side in the pool, as `kontrast_f`
requires. It is determined on the validation half **only in step D**; points
that are not testable are carried as a row of their own, not omitted.

### 5.3 When the hypothesis counts as **demonstrated** — all five conditions

1. **The sign count** of the table in §5.2 is reached, in the direction frozen
   in step C.
2. **Stability:** under each of the four §7 checks at most **two** of the nine
   signs change, and the sign count stays above the threshold of §5.2 under each
   check individually.
3. **Not size-driven:** in the set-size check the median flips at at most **one**
   point (the rule from `PREREG_S5.md` §5.3, unchanged).
4. **The neutral contrast is flat** (`|median z_corr|` < 2) at all nine points.
5. **Not an artefact candidate** under §4.5 — already ensured in §4.6 (c),
   repeated here for completeness: the ΔBaseline value is reported on the
   validation half as well.

If one of the five is missed, the hypothesis counts as **not held**, **S6-AB4**
applies, and §1 closes the axis question.

### 5.4 Equivalence

As in `PREREG_S5.md` §5.6, unchanged: "no difference" only as a **TOST**,
Δ = 2.0, **separately per arm** (arm rule §4). With 3 chondrogenic points in the
validation half the chondrogenic TOST is practically not feasible; that is
reported with its power and not read as a finding.

---

## 6. Step E — the outcomes, **defined here, before any result**

### Outcome 1 — the axis paper

All five conditions of §5.3 are met. Then:

- The hypothesis is reported as a **demonstrated axis**, with the exact p value
  of the validation, the separation of discovery and validation as a carrying
  procedural feature, and the scan balance of §4.4 as context.
- The methods part of the S5 report §6 moves into this paper as a **methods
  chapter**, together with the detection limit of §4.7 and the neutral contrast
  at 52 + 9 point-size combinations.
- The target journal, structure and figure plan are fixed in the S6 report.

### Outcome 2 — the axis question is closed

Every other case, in particular each of the criteria of §7. Then §1 applies. The
methods and negative-finding paper is written, now additionally carried by

- the **detection limit** of §4.7 — the first positive statement about the
  sensitivity of the metric that this work has;
- the **scan balance** of §4.4 — *so many contrasts searched, so many hits
  predicted by noise, so many found*;
- the **artefact balance** of §4.5 — how strongly a free scan on this scale
  actually sorts by the baseline expression.

In **both** outcomes: every number with its MDE80, the power ladder in the text,
the noise expectation, every fallen hypothesis by name, and **one figure**
through `figure_style/publication_style.R`.

---

## 7. Stopping criteria

- **S6-AB1** — the pseudobulk extension changes a reported number qualitatively.
  **Moot**, because the extension is not made under §3; the criterion stays for
  completeness.
- **S6-AB2** — the scan delivers no more hits than noise predicts (§4.4), or no
  candidate reaches `v` of at least 8: **no hypothesis goes into the validation,
  outcome 2.**
- **S6-AB3** — all candidates with `v` of at least 8 are artefact candidates
  under §4.5 or fall under §4.6 (c): **outcome 2**, and that is a finding in its
  own right, worth reporting.
- **S6-AB4** — the validated hypothesis does not hold (§5.3): **outcome 2**, and
  under §1 the axis question is closed.
- **S6-AB5** — the reference implementation proves faulty.

**Step E is carried out under each of these outcomes.**

---

## 8. What stays excluded in this session

- **`AB1`, `S3-AB1`, `S4-AB2`, `S4-AB3` and `S5-AB2` are not revoked.** In
  particular the scan is **not** a second chance for the typology: it tests
  axes, not quadrants, and the lesion class enters no selection rule.
- **No matching to a covariate of the baseline** (§4.5).
- **No second implementation of the metric.** Call it, do not copy it out.
- **No gene sets by hand.** GO through `org.Hs.eg.db` 3.20.0, with the retrieval
  date in the script header and in the CSV.
- **The arm rule §4:** signs, never magnitudes across arms.
- Nothing from the record of fallen hypotheses without addressing the ground of
  its refutation.
- The word **"specific"** for the scissors.
- **The grant proposal** is not a subject of this session.

---

## 9. Release

**2026-08-19, before the first metric of the scan.** The three points at which a
milder version would have been possible are expressly confirmed in the strict
version:

- **§2** — blocks formed by perturbed gene instead of by point, although that
  lowers the number of admissible splits from several thousand to 98 and makes
  the stratification coarser thereby.
- **§4.2** — `NB` stays at 2000, although §7.1 of the brief expressly permits
  the reduction.
- **§4.6 (c) 3** — the exclusion of candidates whose sign count is carried
  solely by the points with a defective null is set in advance as an
  **exclusion** and not as a secondary count as in S5.

This document is thereby in force. From here the addendum rule at the head
applies.

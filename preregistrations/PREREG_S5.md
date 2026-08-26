> Translated from the German original of 2026-08-18, with its addendum of the
> same date. The content, the dates and every number are unchanged.

# Preregistration S5 — axes of convergence

Written on **2026-08-18**, **before any computation of one of the three
metrics**. It builds on §5 of the S5 brief and fulfils step A of the procedure
binding there.

This preregistration extends `PREREG_whole_study.md`; it does not replace it.
Everything stated there — the reference implementation §2, the inclusion and
exclusion criteria §3.1 and §3.2, the screening procedure §3.3, the arm rule §4,
the confounders §7, the duty to use TOST §8 — continues to apply unchanged.
**AB1, S3-AB1, S4-AB2 and S4-AB3 are not revoked.**

From here the same addendum rule applies as in `PREREG_whole_study.md`: changes
only as a dated addendum at the end, with a reason and with a statement of which
numbers were known at the time of the change.

---

## 0. What was known before the date of this document

Completely, and without gloss. There are three items.

**(1) Set sizes and pool intersections — deliberately collected in advance.**
The set-check script ran **before** this dating and counted genes only: set
sizes genome-wide, overlaps, and the number of genes per side in the pool at all
eleven points. It computed no metric, no `z` and no `iv` mean. That is required:
§5.2 of the S5 brief expressly demands for K1 that the pool size be "checked in
advance", and the testability condition of `kontrast_f` (at least 8 genes per
side) is exactly the condition on which S4-AB3 failed. The same disclosure as in
S4 ("disclosed in advance"). The numbers stand in §2.4 below and in the
set-check pool table.

**(2) Individual values of `S_UPR` and `S_ERAD` lie open.** They stand as
`test = "UPR"` and `test = "ERAD"` in the eleven-point metric table, against the
20-class main null, for all eleven points. **They were not opened before this
dating** — neither the CSV columns nor the corresponding block in the log or the
session protocol. That does not change the fact that they **could have been
inspected at any time**, and that is what counts, not the self-report. The
tightened duty to report of §6.2 below therefore applies to H1: the count is
additionally reported over the **newly added** data sets alone, following the
pattern of `PREREG_S4.md` §4.2. For H2 and H3 this problem does not exist — the
sets `S_ABBAU`, `S_ECM_KERN` and `S_SEKRET_LOESLICH` have never been computed
against anything in this work.

**And for all three:** the **contrasts** themselves — the actual metrics of this
preregistration — exist nowhere. Neither `S_STRESS − S_MASCHINE` nor
`S_ABBAU − S_MASCHINE` nor `S_ECM_KERN − S_SEKRET_LOESLICH` has ever been
computed in this work.

**(3) The state of knowledge from S1 to S4**, as it stands in the S1-S3 and S4
reports. In particular it is known that `ARSB` (MPS VI) stands **against** the
original guiding idea both on the apparatus axis (+6.63) and in the distal
decomposition (distal +2.01, biosynthetic +6.36). H2 is expressly formulated
**out of** that known counter-finding — see §2.2. That is not a hypothesis drawn
from a look at the results of the quantity to be tested, but neither is it
independent of them, and it therefore stands here.

---

## 1. Why three hypotheses and not four

§5.1 of the brief allows up to four and requires the number to be chosen **in
advance**, because the correction bound depends on it. **Three** are chosen.

**Bonferroni, α = 0.05 / 3 = 0.0167.** Each hypothesis is additionally reported
uncorrected; the verdict hangs on the corrected value.

The reason for three instead of four is not convenience but that the fourth
hypothesis has no admissible candidate. The reasoning behind the selection
stands in §2.5. **No replacement:** if one of the three drops out before the
computation, the correction bound stays at 0.0167.

---

## 2. The three hypotheses — frozen

All three are **contrasts of two gene sets named in advance against the 20-class
main null**, computed with `kontrast_f()` from the reference implementation,
called unchanged. Single-set `z` values carry no verdict; they run along
descriptively.

The target quantity in all three cases is the interaction term `iv` from §2 of
the whole-study preregistration. A **positive** contrast means: the
first-mentioned set is induced more strongly in the perturbation relative to the
wild type than the second-mentioned.

### 2.1 H1 — the stress response instead of the capacity scaffold (from K1)

> **The statement.** The lesion sits in the stress response of the secretory
> route (UPR/ERAD), not in the capacity scaffold of the apparatus.

| | |
|---|---|
| **metric** | `kontrast_f(iv, ZIEH[["20"]], S_STRESS, S_MASCH_H1, ...)` |
| **`S_STRESS`** | (GO:0030968 ∪ GO:0036503) **∖** `S_MASCHINE` — **158 genes** |
| **`S_MASCH_H1`** | `S_MASCHINE` **∖** (GO:0030968 ∪ GO:0036503) — **3 167 genes** |
| **the direction predicted** | **positive** (z_corr > 0) at the majority of points |
| **origin** | textbook biology: secretory stress induces the UPR, and the UPR is precisely *not* the capacity it regulates. External, through GO, not from these data |

**The reasoning behind the direction.** If matrix production fails on the cargo,
the accumulation of unfolded material in the ER is the classical trigger of the
UPR: the stress response is raised while the capacity scaffold (ribosome,
translation, nucleolus, RNA processing) does not follow. Hence
`S_STRESS − S_MASCH_H1 > 0`.

**A known burden.** The individual values `UPR` and `ERAD` lie open (§0.2) →
the tightened duty to report of §6.2.

**A known weakness, named in advance.** After the separation, `S_STRESS` still
shares **51 genes with `S_DISTAL`** (GO:0030968 is in both). H1 is therefore
**not independent** of the distal axis left undecided in S4. That is not
circularity in the sense of S4-AB3 — the test is against `S_MASCH_H1`, not
against `S_DISTAL` — but a positive result for H1 may **not** be read as an
independent confirmation of the distal axis, and it will not be so read. It
stands here so that it is not negotiable afterwards.

**Testable at 10 of the 11 existing points.** `SERPINA3` KD has only **4** genes
from `S_STRESS` in the pool (below 8) and drops out for H1. That is known in
advance and is carried in the count as a row of its own, not silently omitted.

### 2.2 H2 — degradation against build-up (from K2)

> **The statement.** What is decoupled is not production against cargo but
> production against **degradation**.

| | |
|---|---|
| **metric** | `kontrast_f(iv, ZIEH[["20"]], S_ABBAU, S_MASCH_H2, ...)` |
| **`S_ABBAU`** | (GO:0006914 autophagy ∪ GO:0005764 lysosome ∪ GO:0000502 proteasome) **∖** (`S_MASCHINE` ∪ `S_FRACHT0`) — **891 genes** |
| **`S_MASCH_H2`** | `S_MASCHINE` **∖** the unseparated degradation set — **3 062 genes** |
| **the direction predicted** | **positive** (z_corr > 0) at the majority of points |
| **origin** | MPS VI (`ARSB`) is a **lysosomal** storage disease; the textbook answer to undegraded substrate is an upregulation of the lysosomal-autophagic programme (the TFEB axis). External, through GO |

**The reasoning behind the direction.** If the lesion strikes the degradation
arm or draws on it, degradation rises relative to the biosynthetic capacity
scaffold. Hence `S_ABBAU − S_MASCH_H2 > 0`.

**Honesty about the origin, expressly.** H2 is formulated out of a **known
counter-finding**: `ARSB` stands twice against the guiding idea. That is the
admissible case of §5.1 of the brief — a finding from S1, S2 or S4 that is
**not the quantity to be tested**; the quantity to be tested (`S_ABBAU` against
`S_MASCHINE`) has never been computed. But it is also not the same as an origin
in pure textbook biology, and it therefore stands here and in §0.3.

**A known weakness, named in advance.** `S_ABBAU` overlaps with **82 cell-cycle
genes** and is large and generic at 891 genes. The §7 check "without the cell
cycle" is therefore **not** routine for H2 but the decisive control: if H2 flips
there, it counts as not held.

**Testable at 11 of the 11 existing points.**

### 2.3 H3 — the secretome by destination (from K4)

> **The statement.** The scissors sit **within** the cargo: the matrix cargo
> fails, the soluble cargo does not.

| | |
|---|---|
| **metric** | `kontrast_f(iv, ZIEH[["20"]], S_ECM_KERN, S_SEKRET_LOESLICH, ...)` |
| **`S_ECM_KERN`** | `S_FRACHT0` **∩** `NABA_CORE_MATRISOME` — **188 genes** |
| **`S_SEKRET_LOESLICH`** | `S_FRACHT0` **∖** `NABA_MATRISOME` (core **and** associated) — **2 898 genes** |
| **the direction predicted** | **negative** (z_corr < 0) at the majority of points |
| **origin** | the human-genetics anchor separates disease genes along secretion (OR 4.97); the scissors have never been tested **within** the cargo. The set is external: MSigDB **2026.1.Hs**, `NABA_CORE_MATRISOME` and `NABA_MATRISOME`, through `msigdbr` 26.1.0, retrieved 2026-08-18 |

**The reasoning behind the direction.** All lesions of this study are lesions of
**matrix** production. If something fails, it is the matrix cargo. Hence
`S_ECM_KERN − S_SEKRET_LOESLICH < 0`.

**A known weakness, named in advance — and it is considerable.** The contrast is
a **partition of `S_FRACHT0`**, not a comparison of two independent sets: both
sides lie in the same cargo set, whose overall value (`cargo z`) is already
known at every point. H3 therefore **cannot** show that the cargo fails — that
is already measured. It can only show whether the failure is unevenly
distributed **within** the cargo. It is formulated and reported that way and no
other. In addition the sizes are strongly unequal (188 against 2 898) — the
set-size check of §5.3 is the decisive control for H3.

**Testable at 11 of the 11 existing points.**

### 2.4 The pool sizes collected in advance

From the set-check script, before this dating. Genes per side **in the pool**;
`kontrast_f` requires at least 8 per side.

| point | pool | H1: stress / machine | H2: degradation / machine | H3: ECM / soluble |
|---|---|---|---|---|
| `LAMA5` KO chondrogenic | 9 861 | 60 / 643 | 299 / 601 | 92 / 621 |
| `LAMA5` KO osteogenic | 9 063 | 76 / 692 | 356 / 638 | 75 / 770 |
| `FN1` C123R | 11 461 | 45 / 838 | 280 / 805 | 82 / 695 |
| `FN1` C231W | 11 815 | 53 / 878 | 288 / 841 | 84 / 746 |
| `SERPINA3` KD | 1 347 | **4** / 135 | 33 / 137 | 11 / 118 |
| `MIR181A1HG` KD | 12 039 | 59 / 798 | 300 / 756 | 59 / 666 |
| `LINC01638` KD | 3 941 | 33 / 576 | 207 / 543 | 37 / 396 |
| `ARSB` MPS VI | 7 505 | 72 / 759 | 360 / 713 | 97 / 823 |
| `ACVR1` FOP | 9 132 | 59 / 784 | 242 / 754 | 86 / 612 |
| `RNF4` KD | 8 306 | 29 / 821 | 203 / 795 | 50 / 485 |
| `RB1` +/− | 8 369 | 67 / 1 542 | 355 / 1 487 | 50 / 705 |
| **testable** | | **10 / 11** | **11 / 11** | **11 / 11** |

### 2.5 Why K3, K5 and K6 are **not** on the list

Reasoned **before** the computation, as §5.2 of the brief requires.

- **K3 (`abstumpfung` as a metric) — struck, on two grounds.** First, it is not
  a contrast of two gene sets and so violates the construction rule of §5.1,
  which S4 step B had just shown to be the viable one. Second, it is a
  two-sample problem of `M` against `N` with 6 against 5 points: a Wilcoxon rank
  sum test at 6 against 5 reaches a smallest two-sided p value of **0.0043** —
  formally below 0.0167, but only at **complete** separation of the two classes,
  and §5.2 of the brief itself records that `ARSB` is the only point with a
  positive blunting. The test would not be underpowered but without
  discriminating power. It is not taken up.
- **K5 (degree of differentiation as an axis) — struck.** The arm rule §4
  forbids throwing the arms together; n thereby halves to 5 and 6. At n = 5,
  even 5/5 gives p = 0.0625, which is not below 0.05, let alone below 0.0167.
  **Not demonstrable, under any outcome.**
- **K6 (PanelApp as a gene set) — struck, and the reason is the most important
  of all.** K6 requires a "size- **and expression-matched** remainder" as its
  comparison set. Exactly that matching is demonstrated in the null-diagnostic
  protocol §4 to be **structurally inadmissible** on this scale: the gene-wise z
  standardisation forces cor(baseline, dWT) = −0.566 at the median, negative at
  all eleven points, and every matching to the baseline also matches part of the
  target quantity. The entry in the record of fallen hypotheses says expressly:
  *"This holds for **every** covariate of the baseline on this scale."* K6
  without that matching is an artefact test (disease genes are well-studied,
  highly expressed genes); K6 with it is a revived, already fallen hypothesis.
  Both are excluded. The anchor stays what it is: a finding **outside** the
  transcriptome data.

---

## 3. The power ladder

A paired sign test is an exact binomial test, two-sided, p = 0.5. Computed in
the power script, in full in the power-ladder table. **α = 0.0167** is the
binding column.

| n | required at α = 0.05 | outliers tolerated | **required at α = 0.0167** | **outliers tolerated** |
|---|---|---|---|---|
| **11** | 10/11 | 1 | **10/11** | **1** |
| 13 | 11/13 | 2 | 12/13 | 1 |
| **15** | 12/15 | 3 | **13/15** | **2** |
| 17 | 13/17 | 4 | 14/17 | 3 |
| 19 | 15/19 | 4 | 16/19 | 3 |
| **20** | 15/20 | 5 | **16/20** | **4** |
| 22 | 17/22 | 5 | 18/22 | 4 |
| 24 | 18/24 | 6 | 19/24 | 5 |

**The expected n.** At the time of this dating, n = 11 is secured. From step B
up to 13 lower-priority GSE plus the hits of **one** preregistered search query
(§4.2) are added. An honest expectation, formed from the inclusion rate of the
priority list (**5 of 11 = 45 %**, and the priority list was the *better* one):
**n between 14 and 18**. The lower-priority list is of lower priority because it
looked weaker; a rate above 45 % would be a surprise, not a basis for planning.

**The binding number is entered below as a dated addendum 1 after step B, before
a single metric of §2 is computed.**

**The consequence, accepted in advance.** At n = 11 to 14 an axis hypothesis is
practically only refutable, not demonstrable (a single deviating point
suffices). That is the reason for the lower bound in §5.5.

---

## 4. Step B — the rules of the screen

### 4.1 Unchanged

The thirteen lower-priority GSE from the screening record are screened by the
**unchanged** criteria of `PREREG_whole_study.md` §3.1 and §3.2 (including
addendum 1, `A9`) and the procedure of §3.3: the entry with design, arm, n per
cell, perturbation gene and **lesion class under §3.4** stands **before** any
computation on that data set. **No criterion is relaxed in order to raise n.**

`A9` (single-cell data) stays a ground for exclusion. The reference
implementation is not extended by an aggregation step.

New data sets go into the **same loading block** in `03_metric.R`, not into a
parallel script (§12 of the brief).

### 4.2 The one permitted search query — recorded in advance

It is executed **unconditionally**, regardless of how many of the thirteen are
eligible. That excludes the possibility that its execution is itself a
result-dependent decision.

| | |
|---|---|
| **database** | NCBI GEO DataSets (`gds`), through E-utilities |
| **date** | 2026-08-18 |
| **search string** | `("Homo sapiens"[Organism]) AND ("expression profiling by high throughput sequencing"[DataSet Type]) AND (osteogenic[All Fields] OR chondrogenic[All Fields] OR osteoblast[All Fields] OR chondrocyte[All Fields]) AND (knockout[All Fields] OR knockdown[All Fields] OR CRISPR[All Fields] OR mutant[All Fields] OR silencing[All Fields]) AND (differentiation[All Fields])` |
| **extent** | the first **100** hits by relevance, **all** screened |
| **no follow-up search** | one query, one date, one extent. No second, no extended, no filtered one |

Hits already screened (the priority list, the lower-priority list, the starting
list) are noted as such and not counted twice (`A8`).

### 4.3 The result

One number: **n**, with which step C computes. It enters §3 as addendum 1,
before anything is computed.

---

## 5. Step C — what is computed and what is decided

### 5.1 Per hypothesis, at every point

1. **The metric**: the contrast against the 20-class main null, with MDE80.
2. **Descriptively**: the two single-set `z` values against the same null.
3. **The four §7 checks** of the contrast: a length-matched null, the T2
   residualisation, without the 500 most frequent genes, without the cell-cycle
   genes.
4. **The set-size check** following the pattern of the distal-axis script §6.3:
   **20** random subsets of the larger set at the pool size of the smaller one,
   with seeds keyed to (point, round).
5. **The neutral contrast** (§5.2).

All through `04_load.R` and the unchanged functions `kontrast_f`, `einzel_f`,
`kern`, `mk_zieh` and `mk_zieh_L`. **No second implementation.**

### 5.2 The neutral contrast — the gate from S4 does not carry over automatically

S4 step B calibrated the contrast construction at **one** size class. For every
new set size it is checked afresh: at every point **10** random halvings of the
neutral set `S_NEUTRAL`, at the side sizes of the respective hypothesis, and in
addition 10 halvings at the respective **smaller** side size. Ten instead of
fifty suffices because the construction is already confirmed — that is fixed
here in advance, not chosen afterwards.

**The condition:** a median |z_corr| below 2 at at least `n − 1` points. If it
is missed for a hypothesis, **S5-AB4 applies and that hypothesis drops out.**
The correction bound stays at 0.0167.

### 5.3 When a hypothesis counts as **demonstrated**

All five conditions together:

1. **The sign count in the predicted direction is at least the bound of §3** for
   the n made binding by addendum 1, at **α = 0.0167**.
2. **Stable under all four §7 checks** — no check may turn the verdict. Counting
   is as in S4: the number of points with a sign change; if the count falls
   below the bound, the hypothesis is not demonstrated.
3. **The same in all three counting versions** (§5.4).
4. **Not size-driven**: in the set-size check the median flips at at most
   **one** point (the rule word for word from `PREREG_S4.md` §6.3).
5. **The neutral contrast is flat at its set sizes** (§5.2).

If 1 is missed, the hypothesis is **not demonstrated**. If 1 is met and one of 2
to 5 is missed, it is **likewise not demonstrated** and is reported as "failed
at the confounder criterion", not as a result.

### 5.4 The three counting versions — in advance

For **every** hypothesis, always all three:

- **(a) all points** at which it is testable;
- **(b) only the points newly added in step B**;
- **(c) only the points with an intact null** (all three external neutral sets
  |z| < 2, the column from the null diagnostic; to be determined afresh for new
  points).

Version (a) carries the verdict. (b) and (c) are the counter-checks; if they
deviate, the verdict is not demonstrated under condition 3 of §5.3.

### 5.5 The lower bound — S5-AB1

If n after step B does **not** reach at least 15, no axis hypothesis is
demonstrable in this session. Step C is then **nonetheless computed and reported
in full**, but solely as a **pilot study with an effect estimate and MDE80** —
not as a test, and without any phrasing that suggests a confirmation. That is
**S5-AB1**, and it directly triggers **outcome 2** (§7).

### 5.6 Equivalence — the bound, in advance

Under `PREREG_whole_study.md` §8, "no difference" is carried only as a **TOST**.

- **The equivalence bound Δ = 2.0 in z_corr.** The reason:
  `PREREG_whole_study.md` §5.1 defines |z_corr| < 2 project-wide as **flat**.
  The same bound as in S4, taken over unchanged.
- **And here a tightening relative to S4, expressly named as such:** the TOST
  uses **magnitudes**, and the arm rule §4 permits magnitudes only **within one
  arm**. S4 computed the TOST over all eleven points. In S5 it is computed
  **separately per arm** (osteogenic and chondrogenic), and a cross-arm TOST is
  **not** reported. That is stricter, not milder, and costs power — deliberately
  accepted.
- The **sign test** of §5.3 stays the primary metric and runs over all points:
  it uses only signs, and those are admissible across arms under §4.

---

## 6. Duties to report

### 6.1 All three are computed and all three are reported

Even if the first already has a result. Even if the second obviously fails.
Every fallen hypothesis appears **by name** in the report, with its number.

### 6.2 The tightened duty for H1

Because the individual values `UPR` and `ERAD` lay open (§0.2), the count for H1
is **additionally** reported over the data sets newly added in step B alone —
that is version (b) of §5.4, but for H1 not merely a counter-check but a
mandatory statement in the main text. The pattern: `PREREG_S4.md` §4.2.

### 6.3 An MDE80 with every number

Without exception. A non-rejected null result without an MDE80 is not reported.

---

## 7. Step E — the outcomes, **defined here, before any result**

There are exactly two, and they are of equal standing.

**Outcome 1 — the axis paper.** At least one of the three hypotheses satisfies
**all five** conditions of §5.3, and n is at least 15. Then: a target journal, a
structure, a figure plan, and an honest statement of which of the three did
**not** carry.

**Outcome 2 — the methods and negative-finding paper.** No hypothesis satisfies
§5.3 — or n is below 15 (S5-AB1) and none is demonstrable at all. Then the paper
listed in §2 of the brief is written: the structural limit of the procedure (z
standardisation, cor −0.566, SMD +0.527), the contrast construction as the
remedy (S4 step B), the preregistered negative findings (AB1, S3-AB1, S4-AB2,
S4-AB3, and whatever S5 adds), the human-genetics anchor (OR 4.97), the complete
screening documentation. **With a structure, a target journal and a figure plan
at the same standard as outcome 1.**

**Outcome 2 is not a failure of the session.** It is not phrased as a fallback
position and not reported as one.

In **both** cases the report contains: every number with its MDE80, all three
counting versions, the power ladder in the text, the multiplicity correction,
every fallen hypothesis by name, and **one figure** if a hypothesis has a
result.

---

## 8. Stopping criteria

| | |
|---|---|
| **S5-AB1** | n stays below 15 after step B → no hypothesis demonstrable, step C runs as a pilot study (§5.5), outcome 2 |
| **S5-AB2** | no hypothesis satisfies §5.3 at α = 0.0167 → outcome 2 |
| **S5-AB3** | a hypothesis is not testable (set size, overlap, no external set) → struck, **not replaced**, the correction bound stays 0.0167 |
| **S5-AB4** | the neutral contrast responds at the set sizes of a hypothesis (§5.2) → that hypothesis drops out |
| **S5-AB5** | the reference implementation proves faulty |

**Step E is carried out under each of these outcomes.**

---

## 9. What stays excluded in this session

From the record of fallen hypotheses and §5.2 of the brief, not to be touched
without addressing the ground of their refutation: "the scissors generalise to
ECM genes", cross-arm comparisons, convergence in vitro against in vivo,
cartilage convergence, the contact programme, any statement about the middle of
the pseudotime, **any revival of the null correction rejected in S3 and any
matching to a covariate of the baseline** (see §2.5, K6).

The word **"specific"** is not used for the scissors.

**"Publishable" steers no decision of this session.**

---

## 10. Release

**2026-08-18, before any computation of one of the three metrics in §2.**

The points at which a milder version would have been possible are expressly
confirmed in the strict version:

- **three** hypotheses instead of four, with α = 0.0167 instead of 0.0125 — the
  stricter bound, because the fourth has no admissible candidate (§2.5);
- **TOST separately per arm** instead of over all points, stricter than S4
  (§5.6);
- **a lower bound of n at least 15** for any demonstrability, with the expressly
  accepted result that the session then ends in outcome 2 (§5.5);
- **K6 struck**, although it was the most appealing question in substance,
  because its comparison set rests on a matching already shown to be
  inadmissible (§2.5).

This document is thereby in force.

---

## Addendum 1, 2026-08-18 — the binding n from step B

**The state of knowledge at the time of this addendum:** no metric of §2 has
been computed from any of the seven new data sets. All that is known is the
design information, the occupancy of the cells of the 2 x 2 and the pool size —
the quantities that `PREREG_whole_study.md` §3.2 demands as the only check still
relevant after the computation (`A6`). No contrast, no `z`, no single set.

Step B is complete.

| | |
|---|---|
| so far | 11 |
| newly included | **+7** |
| **n for step C, binding** | **18** |

**S5-AB1 therefore does not apply** (§5.5 requires at least 15).

### The binding row of the power ladder

At **n = 18** and **α = 0.0167** (three hypotheses, Bonferroni):

> **The decision bound is 15 of 18.** Exactly: 15/18 → p = 0.0075; 14/18 →
> p = 0.0309 (above α); 16/18 → p = 0.0013. **Three points may deviate.**

Uncorrected (α = 0.05) the bound would lie at 14/18; that value is additionally
reported under §5.1 but carries no verdict.

### What version (b) of §5.4 now means concretely

The second counting version is the one over the **seven new** points alone. At
n = 7, even 7/7 is only p = 0.0156 — just below α = 0.0167 — and 6/7 at
p = 0.125 is well above it. **Version (b) can practically only refute a
hypothesis, not confirm it**, and is read solely as a counter-check in the sense
of condition 3 of §5.3, not as a test of its own. That applies to the tightened
duty for H1 under §6.2 as well.

### What changes in §5.6

The TOST runs per arm: **osteogenic n = 12, chondrogenic n = 6.** The
equivalence bound stays unchanged at Δ = 2.0.

### What does not change

The three hypotheses of §2, their metrics, their directional predictions, their
conditions of refutation, the multiplicity bound α = 0.0167 and the definition
of the outcomes in §7. **None of it has been touched since the new n became
known.**

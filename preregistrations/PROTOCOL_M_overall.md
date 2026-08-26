> Translated from the German original of 2026-08-21. The content, the dates
> and every number are unchanged. This is a dated internal decision record of
> the study; the journal considered at the time is not the journal the paper was
> finally submitted to.

# Protocol M-OVERALL — the decision gate

Set on **2026-08-21**, following §5.1 of the medical-extension plan of
2026-08-21, after phases A, B and C were complete.

The individual protocols: [`PROTOCOL_M_A.md`](PROTOCOL_M_A.md),
[`PROTOCOL_M_B.md`](PROTOCOL_M_B.md), [`PROTOCOL_M_C.md`](PROTOCOL_M_C.md).
The preregistrations: [`PREREG_M_A.md`](PREREG_M_A.md),
[`PREREG_M_B.md`](PREREG_M_B.md) — both dated before the respective first
download and before the first figure.

---

## 1. Which phase passed its positive control?

| phase | positive control | result |
|---|---|---|
| **A** — the human-genetics anchor | (a) lineage markers in the dysplasia panel: OR 17.1, z +6.84, p 0.0004; (b) the anchor distal against biosynthetic: OR 3.74 (`NOSO`) to 5.75 (`PA309`), p down to 2.9 x 10⁻¹¹ | **PASSED, both parts** |
| **B** — patient concordance | tissue identity per cohort, its own marker set against the others | **PASSED in 4 of 7 cohorts** — gate B met (at least 2) |
| **C** — the diagnosis 2 x 2 | not applicable: phase C is a **screen**, not a measurement. Its result is a design finding and needs no calibration | — |

**At least one level passes its calibration. The rebuild under §5.3 takes
place.** The stopping condition §5.4 does not apply.

---

## 2. Which detection limit applies per phase?

| phase | metric | detection limit |
|---|---|---|
| A | enrichment against the expression- and length-matched null | **OR 1.59** (`GWAS`, 5654 panel genes); OR 1.85 (`KLEIN`, 1130); OR 2.13 (`NOSO_BREIT`, 539); OR 2.37 (`PA309`, 355) |
| B | concordance across patients, a baseline-stratified null | **C 0.576-0.689** per cohort (MDE80); **none** of the four calibrated cohorts reaches its own |
| C | — | the screen: 1424 series, 127 hand candidates, 50 with a diagnosis axis |
| the main part (existing) | a two-set contrast | 0.35 z at 60 genes per side |

---

## 3. What is confirmatory and what exploratory?

| statement | status | basis |
|---|---|---|
| The convergent genes are **not** enriched for dysplasia, short stature or body-height genes, at a detection limit down to **OR 1.59** | **confirmatory** | `PREREG_M_A.md`, dated before the download; gate A passed; 0 of 14 comparisons above the Bonferroni threshold |
| The anchor "distal against biosynthetic secretion" separates dysplasia genes (OR 3.7-5.8) | **confirmatory** (a positive control, reproduced) | `PROTOCOL_M_A.md` §4 |
| In patient cohorts **no** cohort reaches its own MDE80 | **confirmatory** | valid under every omission |
| The reported mean z (−0.69 / −0.94) | **provisional, withdrawn as a core sentence** — a simple mean instead of a study synthesis, carried by one cohort with 2 controls; without it +0.39 / +0.13 | `PROTOCOL_M_B.md` addendum 1 |
| The lesion response lies at chance level in patient cohorts (mean z −0.94) | **confirmatory** | the same computation |
| The one cohort with a directed finding (GSE292600, acromelic dysplasia, z +2.38, p 0.022) | **exploratory** — below its own MDE80, one of four | `PROTOCOL_M_B.md` §4 |
| GSE186141 lies below its null (z −3.95) | **exploratory** — 2 control samples, the highest detection limit of the four; **not a counter-finding**, since the metric does not separate "disagreeing" from "agreeing in the opposite direction" | `PROTOCOL_M_B.md` §5 and addendum 1 |
| 92 % of the candidates with a diagnosis axis have no undifferentiated arm; the two complete 2 x 2 designs are already analysed | **confirmatory** (screening, fully logged) | `PROTOCOL_M_C.md` |
| Three of seven findable patient cohorts measure no skeletal tissue and cannot be calibrated | **confirmatory** (screening) | `PROTOCOL_M_B.md` §3 |

---

## 4. Which core applies?

**Core 2** by §5.2 of the plan — and the condition is met: the positive controls
are passed and the detection limits are measured.

> **What converges in model systems does not converge in patient cells.**
>
> Eighteen perturbation experiments show a shared differentiation programme
> (173 genes against 7.9 expected) and no shared lesion response (7 against
> 8.0). In four calibrated patient cohorts of skeletal disease, **neither** the
> one **nor** the other is above the null: a mean z of −0.69 for the programme,
> −0.94 for the lesion response, and no cohort above its own detection limit. At
> the human-genetics level the programme is not enriched, at **OR 1.00** against
> a detection limit of OR 1.59. And the question can be put at all with the
> available data only in this way, because 92 % of the diagnosis cohorts do not
> co-sequence the undifferentiated reference state.

That is the sharper of the two versions: a statement about the
**transferability of differentiation models to the clinic**, with a measured
detection limit at every level.

**The clinical consequence**, stated expressly in both directions:

- A lesion-spanning downstream biomarker or point of attack is supported at none
  of the tested levels — not in model systems (0.35 z), not in patient cells
  (C 0.58-0.69), not in human genetics (OR 1.59).
- **What that does not say:** that it does not exist. It says that it would have
  to lie below these limits — and it names the design that could find it
  (isogenic series of several lesions in one laboratory, with a co-sequenced
  undifferentiated arm).

---

## 5. What the rebuild has to carry

Following §5.3 of the plan, in this order:

1. the study narrative — the thread under §5.3, target journal Genome Medicine
2. the figures — the new intended order, with the old figures archived
3. the captions, main and supplement
4. the reference list — the medical layer
5. the cover letter
6. the README

Binding, unchanged: confirmatory and exploratory parts never stand in the same
figure; `figure_style/publication_style.R` is not altered; the criticism of the
scale is **not** shortened and **not** softened, it only moves out of the main
thread into S1; the word "specific" is not used for the scissors;
`f4_krankheitsanreicherung.csv` stays unchanged.

---

## Addendum — 2026-08-22: gate B is withdrawn

This addendum changes **no number** of this protocol. It changes the **status**
of phase B, after `PROTOCOL_M_B.md` addendum 2 reassessed the calibration
(§0.5 b and §5 of the donor-resolved plan of 2026-08-22).

**§1 is corrected:** phase B did **not** pass its positive control in a form
that carries a module test. Two of the four "calibrated" cohorts (GSE292600,
GSE77758) were calibrated with `NAIV` — which demonstrates tissue identity, not
the ability of those cells to run the differentiation axis on which the module
is defined. What remains is the OI osteoblast cohort (GSE186141), with two
control samples and the highest detection limit, itself carried as exploratory
in §5, and the enchondromatosis cohort (GSE22855), on the null in both sets. **Gate B counts as not met.**

**§3 is corrected**, in three rows:

| statement | new status |
|---|---|
| "In patient cohorts no cohort reaches its own MDE80" | **stays confirmatory** — a statement about the detection limit, independent of the calibration question, valid under every omission |
| "The lesion response lies at chance level in patient cohorts (z −0.94)" | **exploratory** — the mean is withdrawn under addendum 1, and the calibration under addendum 2 does not carry a module test |
| "GSE292600, z +2.38, p 0.022" | **no longer a reportable module finding** — calibrated with `NAIV` |

What remains confirmatory from phase B is therefore only the **screening and
detection-limit findings**: the complete screen, that three of seven cohorts
measure no skeletal tissue, and that no cohort reaches its own MDE80.

**§4 (core 2) stands unchanged in content**, but the patient sentence within it
is to be read as **exploratory** and no longer as a calibrated level. The
confirmatory weight of the core now rests on phase A (OR 1.00 at a detection
limit of OR 1.59, gate A passed) and on the design finding of phase C (92 %
without an undifferentiated arm). Whether the core keeps its headline is decided by the gate
of **phase D** (`PROTOCOL_M_overall_v2.md`); until then this protocol is the
state of 2026-08-21 with that qualification.

**§5 is replaced** by §9 of the donor-resolved plan of 2026-08-22. The rebuild
carried out on 2026-08-21 stands; the intended order of the figures changes
(phase B moves to S4, phase D becomes Figure 3). The target journal is decided
**only after the gate of phase D** (§9.4 of the plan) — the commitment to
Genome Medicine in §5 point 1 is thereby lifted.

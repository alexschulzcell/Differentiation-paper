> Translated from the German original of 2026-08-22. The content, the dates
> and every number are unchanged. This is a dated internal decision record of
> the study; the journal considered at the time is not the journal the paper was
> finally submitted to.

# Protocol M-OVERALL v2 — the decision gate after phase D

Set on **2026-08-22**, following §8 of the donor-resolved plan of 2026-08-22,
after phases D, E and F were complete.

This protocol **does not replace `PROTOCOL_M_overall.md`**. That document stands
unchanged as the state of 2026-08-21, with its addendum of 2026-08-22. What
stands here is the state **after** the donor-resolved computation.

The individual protocols: [`PROTOCOL_M_A.md`](PROTOCOL_M_A.md),
[`PROTOCOL_M_B.md`](PROTOCOL_M_B.md) (addenda 1 and 2),
[`PROTOCOL_M_C.md`](PROTOCOL_M_C.md), [`PROTOCOL_M_D.md`](PROTOCOL_M_D.md).
The preregistrations: [`PREREG_M_A.md`](PREREG_M_A.md),
[`PREREG_M_B.md`](PREREG_M_B.md), [`PREREG_M_D.md`](PREREG_M_D.md).

---

## 1. Which cells passed their built-in calibration?

`PROTOCOL_M_D.md` §4. The positive control is `dWT` itself: does a data set's
own differentiation contrast find the lineage markers of its own axis?

| | cells | passed |
|---|---|---|
| in total | 14 | **7** |
| of these with a real patient lesion (E2) | 6 | **1** |
| of these an engineering intervention (not E2) | 8 | 6 |

**Half of the published differentiation experiments do not pass their own
positive control** — including the osteogenic arm of our own LAMA5 series. That
is a finding in its own right and belongs in the main text.

---

## 2. Which detection limit applies per level and per statistic?

| level | metric | detection limit (MDE80) | observed |
|---|---|---|---|
| the main part, model systems | a two-set contrast | **0.35 z** at 60 genes per side | convergence 173 against 7.9 expected |
| **D — the programme, pooled** | **S1** (primary) | **0.344** | **0.349** (z +3.00, p 0.0028) |
| D — the programme, per cell | S1 | 0.288 to 0.455 | 2 of 7 above it |
| **D — the programme, study synthesis** | S1 | **+1.151** — **unreachable** on a scale bounded by 1.0 | +0.683 (z +1.66, p 0.249) |
| D — the programme, synthesis | S2 | **the null degenerates** (S2 is invariant under donor flips) | no number |
| D — the lesion response | S1 | **not determinable** — one calibrated cell | **no number** |
| D — the engineering response | S1 | 0.122 | 0.089 (z +1.44, p 0.153) |
| A — human genetics | enrichment, matched on expression and length | **OR 1.59** (GWAS, 5654 genes) | OR 1.00 |
| B — patient against control | concordance, stratified on the baseline | C 0.576 to 0.689 per cohort | no cohort above it |
| C — the diagnosis 2 x 2 | screening | — | 1424 series, 92 % without an undifferentiated arm |
| D — screening by design | screening | — | 22 series, 2 with a complete 2 x 2 and a skeletal axis |

---

## 3. Which level of the plan has been reached?

**None.** `PROTOCOL_M_D.md` §7: the stopping criterion of the preregistration
applies (it demands at least 6 calibrated lesion cells from at least 3 studies;
1 from 1 was reached), and none of the four levels of the table fits, because
the lesion half cannot be tested at all.

The applicable finding:

> **H1 is confirmed exploratorily**, **H2 and H3 cannot be tested.** The
> scissors have **not been measured** in a donor-resolved way — which is
> something different from "not found".

By §9.1 of the plan, the thread for **"weak / not found"** therefore applies:
the fallback core §7.1 (the map of levels) and §7.2 (the map of detection
limits) carry the main figures.

---

## 4. What is confirmatory and what exploratory — the complete state

| statement | status | basis |
|---|---|---|
| The convergent genes are **not** enriched for dysplasia, short stature or body-height genes, at a detection limit down to **OR 1.59** | **confirmatory** | `PROTOCOL_M_A.md`, gate A passed |
| The anchor "distal against biosynthetic secretion" separates dysplasia genes (OR 3.7-5.8) | **confirmatory** (positive control) | `PROTOCOL_M_A.md` §4 |
| **92 % of the diagnosis cohorts have no undifferentiated arm** (1424 series) | **confirmatory** (screening) | `PROTOCOL_M_C.md` |
| **Searching by design instead of by entity leaves 22 series — two with a complete 2 x 2 and a skeletal axis** | **confirmatory** (screening, fully logged) | `PROTOCOL_M_D.md` §3 |
| **7 of 14 donor-resolved cells do not pass their own positive control** | **confirmatory** (a measurement with a fixed rule, before the main computation) | `PROTOCOL_M_D.md` §4 |
| In patient cohorts **no** cohort reaches its own MDE80 (C 0.576-0.689) | **confirmatory** — a statement about the detection limit, valid under every omission | `PROTOCOL_M_B.md` §4 |
| Three of seven findable patient cohorts measure no skeletal tissue | **confirmatory** (screening) | `PROTOCOL_M_B.md` §3 |
| **The differentiation programme runs in the same direction across donors** (S1 0.349, z +3.00, p 0.0028; without the module-forming cells z +4.51) | **exploratory** — the stopping criterion applies, the synthesis is not above the null, and only 5 carrying cells from one donor cohort | `PROTOCOL_M_D.md` §5.1, §6 |
| The engineering response is **not** above its null on the primary statistic (z +1.44) | **exploratory** | `PROTOCOL_M_D.md` §5.3 |
| **The lesion response across donors** | **not testable — no number reported** | `PROTOCOL_M_D.md` §5.2 |
| Phase M-B as a whole (patient against control) | **exploratory**, gate B withdrawn | `PROTOCOL_M_B.md` addendum 2 |
| Chromatin carries the programme (z +4.88 above MDE80), the promoter methylome does not | **exploratory** | `PROTOCOL_orthogonal_layers_2026-08-20.md` |
| Convergence in model systems (173 against 7.9) | **exploratory** (as before) | Figure 2 |

---

## 5. The core, restated

**Core 2 in the version of 2026-08-21 no longer holds.** It rested on "in four
calibrated patient cohorts neither the one nor the other is above the null" —
the calibration of those four has been withdrawn (`PROTOCOL_M_B.md`
addendum 2), and the donor-resolved test could not examine the lesion side.

**The core that carries after phase D** — stated positively, with a measured
limit at every level:

> **A differentiation programme can be recovered across measurement levels and
> across donors. The lesion response can be recovered at no level — and at the
> level that would count, it cannot even be tested with the available data.**
>
> The fixed 173-gene module is recoverable in chromatin (z +4.88 above MDE80),
> between the donors of a calibrated differentiation experiment (S1 z +3.00,
> without the module-forming cells z +4.51) and in the external triangulation
> (synthesis p 0.0002). It is **not** recoverable in the promoter methylome, not
> in human genetics (OR 1.00 at a detection limit of OR 1.59) and not in the
> patient against control contrast (no MDE80 reached). The lesion response is
> above its null at no tested level. And the one level at which the question
> would be decided — the same donors, the same axis, a real diagnosis, a
> co-sequenced undifferentiated arm — exists in the published literature in **one**
> calibrated cell.

**The clinical consequence**, in both directions:

- A lesion-spanning downstream biomarker or point of attack is supported at no
  tested level — not in model systems (0.35 z), not in patient cells
  (C 0.58-0.69), not in human genetics (OR 1.59).
- **What that does not say:** that it does not exist. For the decisive level the
  statement is sharper and more uncomfortable: **the data that could find it
  have not been collected.**
- The corresponding design has stood since phase C and is made precise by
  phase D: isogenic series of several patient lesions in one laboratory, with a
  co-sequenced undifferentiated arm, **and with a calibration per donor** — half of the
  existing experiments fail at the last of these.

---

## 6. The target journal — decided here

By §9.4 of the plan the decision is taken only at this point. What has been
reached is **neither strong nor middling**.

- **Genome Medicine: no.** The commitment made in `PROTOCOL_M_overall.md` §5 is
  lifted.
- **First choice: Disease Models & Mechanisms.** Model validity is the core
  subject there, and that is exactly the core here: what a model transfers, what
  it does not, and what one measures that against. The finding "half of the
  published differentiation experiments do not pass their own positive control"
  belongs in that journal.
- **Second choice: eLife (Reviewed Preprint)** — it carries a methodological
  contribution with measured detection limits well, even without a positive
  headline.
- **Safe harbour: BMC Genomics**, as before.
- PLOS Computational Biology remains possible but fits less well, because the
  contribution is not algorithmic.

---

## 7. What the rebuild has to carry

Following §9 of the plan, in this order: the study narrative, then the figures,
then the captions, then the reference list, then the cover letter, then the
README, then the guards.

Binding, unchanged: confirmatory and exploratory parts never stand in the same
figure; old figures are archived, nothing is overwritten;
`figure_style/publication_style.R` is not altered, and `PUB_DIR` is overwritten
after the `source()`; the criticism of the scale is **not** shortened and
**not** softened; the word "specific" is not used for the scissors;
`f4_krankheitsanreicherung.csv` stays unchanged; phases M-A, M-B and M-C are not
deleted but only extended.

**And the rule from §7.3 of the plan:** at least **two** of the five main
figures carry a positively stated result. The map of levels (§7.1) and the map
of detection limits (§7.2) do that.

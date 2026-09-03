> Translated from the German original of 2026-08-22. The content, the dates
> and every number are unchanged.

# Preregistration M-E — are the scissors a finding or a difference in size?

**Dated 2026-08-22, before the first statistic of this phase.** The occasion is
the objection, raised the same day, that `dWT` and `iv` are not comparable and
that the narrative may therefore be selling a truism as a finding.

---

## 0. The objection that triggers this phase

The paper claims: **the differentiation programme converges, the lesion
response does not.** What is being compared are two quantities that differ in
at least four respects, none of which has anything to do with "shared against
individual":

1. **Order.** `dWT` is a first-order effect, `iv` a difference of differences —
   second order, with about twice the variance from its construction alone.
2. **Expectation.** Osteogenic differentiation is the same programme in
   eighteen laboratories; eighteen different lesions have no reason to be.
3. **Ceiling.** Both gene sets were formed by the same consistency procedure,
   but the attainable upper bound is lower on the `iv` side.
4. **Noise.** Part of the difference is measurement error, and **nowhere** has
   anything been computed against that so far.

**The sharpest expression of the problem:** in this work `dWT` is
simultaneously the **positive control** (phase M-D calibrates every cell on it)
and the **finding** (the scissors). The two together are not tenable.

---

## 1. What is known at the time of this preregistration

Stated explicitly, so that nothing is presented as blind:

- **every** number of phases M-A to M-D, including the protocols;
- a **diagnosis of the orders of magnitude**, computed on 2026-08-22 from
  `derived_data/R_intern/R_interne_genkarte.csv` and reproduced here in full:

  | | `dWT` | `iv` |
  |---|---|---|
  | median absolute effect, all genes | 0.394 z | 0.194 z |
  | median absolute effect, module genes | 1.135 z | 0.200 z |
  | share of genes above 0.35 z | 54.7 % | 23.8 % |
  | genes with absolute consistency ≥ 0.9 | 0.17 % | 0.01 % |

**Not known** is any result of the three tests below. Neither a pseudo-`iv`,
nor an SNR-matched `dWT` convergence, nor a cross-control has ever been
computed.

---

## 2. What this phase does **not** do

- **No new axis and no new gene set.** The fixed 173-gene module and the
  size-matched lesion set stay unchanged. `PREREG_S6.md` §1 applies.
- **No second implementation.** Loading uses the reference loader
  `11_load_18_datasets.R`; the contrasts have the same algebra as `kern()` in
  `12_metric_reference.R`; the between-donor statistic comes from
  `00_shared/_module.py` (`leiter`).
- **No matching to a covariate of the baseline.** The "not a flat null" guard
  applies unchanged.
- **No change to phases M-A to M-D.** This phase changes the
  **interpretation**, not the numbers.

---

## 3. The objects

**The primary object:** the **eighteen perturbation data sets**, loaded with
`11_load_18_datasets.R`, matrix `Z` (per-gene z per data set) and `meta`
(`genotype` in {WT, KO}, `condition` in {undiff, diff}) — unchanged.

**The gene universe**, fixed in advance: the genes measurable in **at least 16
of 18** data sets (the universe `U1` from `31_derive_matrix_programme.py`, 10 177 genes
there). No pool filter and no filtering by effect size.

**The convergence rule**, fixed in advance and unchanged from the main part: a
gene converges if its sign is the same in **at least 90 %** of the data sets in
which it is measurable.

**The secondary object for test 3:** the seven calibrated cells of phase M-D
(`derived_data/M_donoren/zellen.pkl`).

**Seed 20260823. B = 200 draws.**

---

## 4. Test 1 (primary) — the noise floor for quantities of `iv` shape

**The question.** How much convergence does a quantity of the algebraic form of
an interaction term produce when **no lesion at all** is contained in it?

**The construction.** Per data set *d* and per draw *b*, **individual samples
are drawn without replacement** from their cells — so both quantities have an
identical sample count and an identical noise structure:

    iv_1x1(d,b)        = [KO_diff(s1) − KO_undiff(s2)]
                         − [WT_diff(s3) − WT_undiff(s4)]

    pseudo_iv_1x1(d,b) = [WT_diff(a1) − WT_undiff(a2)]
                         − [WT_diff(b1) − WT_undiff(b2)],   a1≠b1, a2≠b2

    dWT_1x1(d,b)       = WT_diff(s3) − WT_undiff(s4)        (the reference)

`pseudo_iv_1x1` has the same algebra, the same order, the same sample count and
the same degrees of freedom as `iv_1x1` — it simply crosses **no genotype
boundary**. It is therefore the most honest null these data allow: better than
any analytical one, because it brings the real measurement noise of the
respective data set with it.

**The precondition** is at least 2 WT samples per condition. The inclusion rule
`A4` (n ≥ 2 per cell) guarantees that; any data set that unexpectedly fails it
is removed from **both** arms of the test and named in the protocol.

**The statistic.** Per draw *b*: the number of convergent genes across the 18
data sets, separately for `iv_1x1` and `pseudo_iv_1x1`. Reported are the mean
and the 2.5th and 97.5th percentiles over the 200 draws.

### Decision rule for test 1

| result | conclusion |
|---|---|
| mean(`iv_1x1`) **≤** 97.5th percentile(`pseudo_iv_1x1`) | The lesion response has been measured **at the noise floor of its own construction**. The sentence "the lesion response does not converge" is **not a biological statement**. → **narrative B**, independently of test 2 |
| mean(`iv_1x1`) **>** 97.5th percentile(`pseudo_iv_1x1`) | `iv` carries structure that is not noise. → continue with test 2 |

---

## 5. Test 2 — the SNR matching

**The question.** Does `dWT` still converge when its **shared component** is
pulled down to the signal-to-noise ratio of `iv`?

**Why not simply rescale:** the convergence rule is sign-based and therefore
**scale-invariant**. Mere rescaling changes nothing. What counts is the ratio
of the shared component to the spread.

**The construction**, on the existing per-data-set tables `20d_gene_*.csv` (per
gene and data set: `dWT`, `iv`):

    m_g      = mean over data sets
    s_g      = SD over data sets
    SNR(X)   = median_g( |m_g(X)| / s_g(X) )
    k        = SNR(iv) / SNR(dWT)

    dWT*_d,g = k · m_g(dWT) + ( dWT_d,g − m_g(dWT) )

The **shared** component is compressed to the level of `iv`, while the
**individual** spread is left untouched. Then the same convergence count.

**The null** is the sign-flip null of the main part: the sign is flipped per
data set jointly (which preserves the correlation between genes), 2 000 rounds,
the same seed.

### Decision rule for test 2

| result | conclusion |
|---|---|
| count(`dWT*`) **>** the 97.5th percentile of its flip null **and** count(`dWT*`) **≥ 2 x** count(`iv`) | The asymmetry survives the matching. → **narrative A** |
| otherwise | The asymmetry cannot be separated from size and noise. → **narrative B** |

---

## 6. Test 3 — cross-control of set against quantity (descriptive)

Separates "it is the gene set" from "it is the quantity measured". It
**decides nothing** but is reported in full.

- On the 18 data sets: the convergence count of `dWT` restricted to the genes
  of the **lesion set**, and of `iv` restricted to the genes of the
  **programme set**.
- On the seven calibrated cells of phase M-D (`_module.leiter`, unchanged): S1
  of the `dWT` vectors on the **lesion genes**, and S1 of the engineering `iv`
  vectors on the **programme genes**.

---

## 7. The two narratives, written out in advance

So that afterwards there is no negotiating over what the result "means":

**Narrative A** — *"A differentiation programme reproduces across data sets,
measurement levels and donors. A lesion response of **comparable size** does
not."* The qualifier "of comparable size" is admissible only if test 2 carries
it. The scissors remain a results section, **not** the title.

**Narrative B** — *"A shared downstream lesion response is not detectable down
to effect size X. Here is the measuring instrument that works at that size on
six levels — and here is the reason why the data for more do not exist."*
`dWT` is then **throughout** the positive control, on every level, consistent
with phase M-D. The scissors disappear as a claim and become the **study
design** (a test plus a built-in positive control).

---

## 8. Rules of honesty for this phase

- **Both narratives are admissible results.** Neither is preferred.
- **All three tests are reported**, even if test 1 already forces the decision.
  No selective reporting.
- **No threshold is lowered after the numbers are known**, and no test is
  declared primary afterwards. The primary test is **test 1**.
- If the decision falls to **narrative B**, the manuscript, the figure legends,
  the cover letter and the README are rewritten accordingly, and the guard
  entry on the scissors is **weakened** — not deleted.
- The numbers of phases M-A to M-D stand. This phase changes their
  **interpretation**, not their content.

---

## 9. Output

`derived_data/M_kalibrierung/`: `test1_rauschboden.csv`,
`test1_ziehungen.csv`, `test2_snr.csv`, `test3_kreuz.csv`, `55_log.txt`.
Protocol: `preregistrations/PROTOCOL_M_E.md`.

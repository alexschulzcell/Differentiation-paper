> Translated from the German original of 2026-08-22. The content, the dates
> and every number are unchanged.

# Protocol M-E — the scissors are a difference in size

Computed **2026-08-22**, following `PREREG_M_E.md` (dated before the first
figure of this phase). Scripts: `08_disease_gene_orthogonality/40_noise_floor_contrasts.R`,
`08_disease_gene_orthogonality/41_noise_floor_tests.py`. Seed 20260823, B = 200 draws,
2 000 flip rounds. All numbers are in `derived_data/M_kalibrierung/`.

**The result in two sentences:**

> **Both decision tests fall to thread B.** The lesion response has been
> measured at the noise floor of its own construction, and the convergence of
> the differentiation programme disappears as soon as its signal-to-noise ratio
> is compressed to that of the lesion response.
>
> **The sentence "the programme converges, the lesion response does not" is
> therefore not a biological statement about shared against individual
> responses. It is a statement about effect size and noise.**

---

## 1. Environment note — what this run required

The reference scripts were run against the data layout of the analysis
sessions, which at the time of this run lay in a separate archive. They were
mounted at their expected location for the run, so that the reference
implementation itself stayed literally unchanged, and unmounted again
afterwards.

**That is a legacy of the layout, not a solution.** The absolute paths in
`10_load_reference_metric.R`, `14_geo_matrices_s5_format.R` and `11_load_18_datasets.R` are a trap for anyone who
receives the repository and have to be moved to relative paths before
submission — together with the Zenodo deposit, not before it, because the
reference implementation would otherwise differ between two runs.

---

## 2. The data situation

`40_noise_floor_contrasts.R` loaded all **18** data sets through the reference loader
`11_load_18_datasets.R` and formed 200 draws of the three single-sample contrasts per
data set (`test1_datensaetze.csv`).

**All 18 are usable**: each has at least 2 wild-type samples per condition, as
the inclusion rule `A4` guarantees. No data set had to be excluded.

The universe, fixed in advance: **10 177 genes** measurable in at least 16 of
18 data sets (`U1` from `31_derive_matrix_programme.py`). The convergence rule is
unchanged: the same sign in at least 90 % of the data sets in which the gene is
measurable.

---

## 3. Test 1 (primary) — the noise null floor

`test1_rauschboden.csv`, `test1_ziehungen.csv`.

Per draw, three quantities were formed from **single samples** — the same
algebra, the same sample count, the same degrees of freedom, except that the
pseudo quantity crosses no genotype boundary:

| quantity | convergent genes, mean | SD | 2.5 % | 97.5 % |
|---|---|---|---|---|
| `dWT_1x1` (the reference) | **66.8** | 30.1 | 26.0 | 143.0 |
| `iv_1x1` (the lesion response) | **6.0** | 5.9 | 0.0 | 22.1 |
| `pseudo_iv_1x1` (**no lesion**) | **4.3** | 4.4 | 0.0 | 15.0 |

**The decision rule:** mean(`iv_1x1`) = 6.0 is at most the 97.5th percentile of
`pseudo_iv_1x1` = 15.0 → **thread B.**

The ratio of `iv` to `pseudo_iv` is **1.40** on average. The lesion response
therefore produces about forty per cent more convergence than a contrast
**without any lesion** — and that lies entirely within the spread of the noise
floor.

**What this test also shows: the instrument works.** At the same resolution,
with the same single samples, `dWT_1x1` finds **66.8** convergent genes on
average — more than fifteen times the noise floor. The null finding on the `iv`
side is therefore not a failure of the procedure but a statement about the size
of what would be there to measure.

> **The sentence that follows from test 1:** the reported non-convergence of the
> lesion response (7 genes against 8.0 expected) does not differ from what a
> contrast of the same construction **without a lesion** produces. It is not an
> observation about lesions.

---

## 4. Test 2 — the SNR adjustment

`test2_snr.csv`. Computed on the existing per-data-set tables
`20d_gene_*.csv`, with a flip null of 2 000 rounds (signs turned jointly per
data set, which preserves the correlation between genes).

The signal-to-noise ratio, median over the genes:

    SNR(dWT) = 0.3002      SNR(iv) = 0.1941      k = 0.6464

The **shared** part of `dWT` was compressed by `k`; the **individual** spread
was left untouched.

| quantity | convergent | flip null | 97.5 % | verdict |
|---|---|---|---|---|
| `dWT` observed | **120** | 4.8 ± 11.0 | 27.0 | above the null |
| `iv` observed | **4** | 4.8 ± 6.4 | 23.0 | on the null |
| **`dWT*` SNR-adjusted** | **11** | 4.7 ± 9.6 | 28.0 | **on the null** |

**The decision rule:** `dWT*` = 11 lies neither above the 97.5th percentile of
its null (28.0) nor at at least twice the `iv` count (8) → **thread B.**

**The real lesson lies in the steepness.** It is enough to lower the
signal-to-noise ratio by **35 %** for the convergence to collapse from 120 to
11 — below the null threshold. The convergence count is therefore a **highly
non-linear function of the SNR**, practically a step function. Comparing two
quantities with different SNR through that count says almost nothing about
"shared against individual" and almost everything about their effect size.

**One deviation deserves naming.** The main part reports **173** convergent
`dWT` genes against 7.9 expected and **7** for `iv` against 8.0. This
recomputation finds **120** and **4**. The difference is the universe: the main
part counts over all genes with a valid value, whereas here the universe `U1`
fixed in advance is used (at least 16 of 18 data sets, 10 177 genes). The
direction and the decision are untouched by that; the numbers are **not** to be
compared one to one with those of the main part, and that is stated here
instead of being smoothed away.

---

## 5. Test 3 — a cross-check of set against quantity (descriptive)

`test3_kreuz.csv`. Named in advance as **not decisive**.

| gene set | quantity | genes in the universe | convergent | null (97.5 %) |
|---|---|---|---|---|
| the programme set (173) | `dWT` | 115 | **96** | 0.0 |
| the programme set (173) | `iv` | 115 | **0** | 1.0 |
| the lesion set (173) | `dWT` | 162 | **1** | 1.0 |
| the lesion set (173) | `iv` | 162 | **4** | 1.0 |

The two sets behave entirely separately: the programme set converges on `dWT`
(96 of 115) and not at all on `iv` (0), and the lesion set on neither quantity
to any notable degree.

**These numbers are largely tautological** and are reported as such: the
programme set **was** selected for `dWT` consistency and the lesion set for `iv`
consistency, both from the same 18 data sets. That the selection reappears in
the count is no confirmation. The one non-trivial point: that the lesion set
delivers only 4 genes on `iv` as well, although it was selected for exactly
that — which is the same statement as test 1, at the level of the gene set.

---

## 6. What follows from it

**Thread B applies.** Spelled out as in `PREREG_M_E.md` §7:

> A shared downstream lesion response is not demonstrable down to a measurable
> effect size. `dWT` is **throughout the positive control** — at every level,
> consistent with phase M-D, where every cell is calibrated against it. The
> scissors disappear as a claim and become the **study design**: a test plus a
> built-in positive control.

**What stands, unchanged:**

- the criticism of the scale (Figure S1) — untouched by this phase;
- the human genetics (OR 1.00 at a detection limit of OR 1.59) — there **two
  gene sets are compared against the same matched background**, and the
  comparability is clean;
- the design findings (92 % without an undifferentiated arm; 22 series along the design);
- the calibration balance of phase M-D (7 of 14 cells do not pass);
- every measured detection limit.

**What falls:**

- "the programme converges, the lesion response does not" as a **biological**
  statement and as a title;
- the reading of the internal core (173 against 7.9 beside 7 against 8.0) as
  evidence of an asymmetry between a shared and an individual response;
- the claim that `dWT` is finding and positive control at once.

**What is newly added and can be stated positively:**

- **A noise floor this field does not use.** Test 1 delivers a procedure with
  which, for any interaction term, one can say how much of its convergence
  arises without any perturbation — built from the replicates of the control
  arm, without additional data.
- **The demonstration that the convergence count is almost a step function of
  the SNR** (120 → 11 at a 35 % loss of SNR). That is a warning to any work
  that places the convergence counts of two effects of different size side by
  side.

---

## 7. What was expressly not done

- **No new axis, no new gene set.** Module and lesion set unchanged.
- **No second implementation.** Loaded with `11_load_18_datasets.R`, computed with the
  algebra of `kern()`, with the between-donor statistics from `_module.py`.
- **No threshold shifted once the numbers were known**, and no test declared
  primary after the fact. Test 1 was and is the primary one.
- **No test suppressed**: tests 2 and 3 are reported in full, although test 1
  had already forced the decision.
- **No number of phases M-A to M-D changed.** This phase changes the reading,
  not the content.

---

## 8. Output

`derived_data/M_kalibrierung/`: `test1_datensaetze.csv`,
`test1_rauschboden.csv`, `test1_ziehungen.csv`, `test2_snr.csv`,
`test3_kreuz.csv`, `55_log.txt`, and `kontraste/` with 54 sign tables (three
quantities x 18 data sets, 200 draws each).

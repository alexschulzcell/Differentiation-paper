> Translated from the German original of 2026-08-18. The content, the dates
> and every number are unchanged.

# Preregistration S4 — the distal axis and the secretory transcription-factor arm

**Dated 2026-08-18, written before any computation on a target quantity of S4.**

`PREREG_whole_study.md` (§4 the arm rule, §7 the confounders, §8 TOST, §12 the
duty to report) and all stopping decisions remain in force. **AB1 and S3-AB1
are not revoked.** S4 is not a second chance for the typology but two new
questions on a different level of analysis.

---

## 0. The two questions

- **A** — does the distal half of the secretory apparatus carry the
  information, and is the biosynthetic half the one that merely follows?
- **B** — does the secretory transcription-factor arm (XBP1, ATF6, CREB3L1,
  MIA3) behave like the apparatus or like the cargo?

---

## 1. What is known at the time of this preregistration

**A mandatory statement. Honesty rather than a fiction of purity.**

### 1.1 Numbers that touch the target quantity of question A

- **`ARSB` / `GSE218101`: distal +1.99, biosynthetic +6.29.** A
  **counter-example** to question A.
- **S2 / plasma cell:** distal +5.91, biosynthetic +6.87 (day 0 to 6); over
  the time course distal rises monotonically while biosynthetic falls back
  after day 3. S2 is **not** one of the eleven data sets.
- **The human-genetics anchor:** distal 39/523 against biosynthetic 35/2 192,
  OR 4.97, p 6.73e-11.
- Under the correction rejected in S3, S2 was measured as distal +3.93 against
  biosynthetic −4.57. That null is **not** used in S4.
- **The remaining ten data sets are unread.** The columns `nur_distal` and
  `nur_biosyn` from `03_metrik_elf.csv` have not been opened up to the date of
  this file, neither as a CSV nor through a script nor as a summary.

**The consequence, fixed in advance:** the result of A is counted **twice** —
once over all eleven data sets, once over the **ten unread** ones (excluding
`ARSB`). If both verdicts agree, the reservation is void. If they differ, **the
version over the ten applies**.

### 1.2 Numbers that touch the target quantity of question B

None. No regulon has yet been computed against any data set.

**What was checked before this preregistration and belongs stated honestly:**
the availability and the size of the regulon databases. Measured on
2026-08-18, before the date of this file:

| source | version | edges of the four factors |
|---|---|---|
| **CollecTRI** via `decoupleR` 2.12.0 / the `OmnipathR` cache | retrieved 2026-08-18 | **126** |
| **DoRothEA** `dorothea` 1.18.0, `dorothea_hs`, confidence A to C | package version | 47 (XBP1 16, ATF6 10, CREB3L1 21, MIA3 0) |

Those are **set sizes, not results** — no z value, no sign, no contrast. They
are recorded here nonetheless, because the choice of the primary database was
made after this was known. The choice is justified in §5.1 and is **binding
from now on**.

### 1.3 A runtime measurement

Also measured before the date of this file, on a **halving of the synaptic
neutral set** at data set 1 (`LAMA5` chondrogenic), z_corrected −1.00: one call
of `kontrast_f` with NB 2000 and NVIF 100 takes 6.8 s; loading all eleven data
sets takes 46 s. This number is a planning quantity and not a target quantity;
it feeds into §3.4 (the scope of B1).

---

## 2. What does **not** change in the machinery

- **No second implementation.** The entry point of every computation is the
  reference loader; the computation uses only `kern()`, `mk_zieh()`,
  `mk_zieh_L()`, `einzel_f()` and `kontrast_f()` from `03_metric.R`, unchanged.
- **The main null stays the main null:** `ZIEH[["20"]]`, 20 induction classes.
  The correction rejected in S3 (level x baseline) is **not** used and not
  revived.
- **Seed 20260818, NB 2000, NVIF 100, THR 0.5, pool rule A6** unchanged.
- **VIF correction always. MDE80 with every number.**
- **The arm rule §4:** what is compared are **signs**, never magnitudes across
  arms.
- **Gene sets unchanged** from `03_metric.R`: `S_DISTAL` 545 genes,
  `S_BIOSYN` 2 711, `S_FRACHT0` 3 621, `S_MASCHINE` 3 205, `S_NEUTRAL` 614,
  `S_ZYKLUS` 1 907. Mutually disjoint, as made there.

---

## 3. Step B — the calibration of the contrast. **The gate.**

The load-bearing claim of S4 is: *a contrast of two sets against the same null
is more robust than the two individual values, because common bias cancels
out.* **That is a claim, not a fact**, and it is tested before any result is
interpreted.

### 3.1 B1 — the neutral contrast

For each of the three external neutral sets — `GO:0007268` (synaptic),
`GO:0007608` (smell), `GO:0007601` (vision), each made disjoint from the metric
sets, the cell cycle, ossification and cartilage development — the set is
randomly halved in the pool at each of the eleven data sets and
`kontrast_f(iv, ZIEH[["20"]], H1, H2, ...)` is computed. The expectation is
z_corrected close to 0 everywhere.

**Scope: 50 halvings per set and data set**, seed 20260818. NB is **not**
reduced. The three sets run as three separate background processes.

A data set counts for a set only if, after `kontrast_f`, both halves have at
least 8 genes in the pool; otherwise it counts there as **not evaluable** and
is not treated as a violation (it is missing from both numerator and
denominator for that set). A known risk: `GO:0007608` in the `LINC01638`
knockdown, `GO:0007601` in the `SERPINA3` knockdown.

### 3.2 B2 — the decisive comparison

At the **four data sets with a defective null** (`FN1` C231W, `FN1` C123R,
`RB1`, `RNF4`) the individual set shows an excursion. If the neutral contrast
shows an excursion there as well, the disturbance does **not** cancel out.
These four data sets are explicitly included in gate condition 1 (§3.5).

### 3.3 B3 — the confounders from §7

Applied to the contrast distal minus biosynthetic, at all eleven data sets: a
length-matched null (`mk_zieh_L`), without the 500 most frequent genes, without
cell-cycle genes, and T2 residualisation. These run along in step C (§5).

### 3.4 The comparison quantity for gate condition 2

"The individual-set z at the same data sets" means:

- **(i)** the z_corrected of the **full** neutral set against the same null,
  `einzel_f(iv, ZIEH[["20"]], S, ...)`, per data set — one number per data set;
- **(ii)** the z_corrected of the **individual halves** as individual sets,
  `einzel_f(iv, ZIEH[["20"]], H1, ...)` and `H2`, that is 100 values per data
  set, so that contrast and individual value are compared at **equal set
  size**.

**(ii) is the decisive comparison quantity**, because it holds the set size
constant; (i) is reported in addition. Both are fixed in advance, so that the
more favourable one is not chosen afterwards.

**Spread** means the standard deviation of z_corrected **over the eleven data
sets**, computed per halving and then averaged over the 50 halvings. For (ii)
likewise, with the 100 half-sets per halving round.

### 3.5 The gate — passed if **both** hold

1. **The neutral contrast lies at |z_corrected| < 2 in at least 10 of 11 data
   sets**, including the four with a defective null. What counts is the
   **median of the 50 halvings** per data set, for the **synaptic** set (the
   only one evaluable at all eleven). For smell and vision the same count is
   reported **in addition** over the data sets evaluable there; a violation
   there is a warning sign that appears in the protocol but does not close the
   gate on its own.
2. **The mean spread of the neutral contrast over the 50 halvings is smaller
   than that of the individual-set z at the same data sets**, per §3.4 (ii). A
   tie counts as **not** passed.

If this is missed, **S4-AB1** applies: the contrast is not more robust than the
individual values, the justification for S4 falls away, and only step E
remains. **This outcome is not circumvented.** In particular, no other neutral
set and no other number of halvings is then sought until the gate falls.

---

## 4. Question A — design

- **The statistic:** `kontrast_f(iv, ZIEH[["20"]], S_DISTAL, S_BIOSYN, pool,
  rho_bar, ...)` — the contrast **distal minus biosynthetic** against the
  20-class main null, with VIF correction and MDE80. One new row in the
  existing machinery.
- **The individual values** `nur_distal` and `nur_biosyn` run along
  **descriptively**; they are **not** the statistic and carry no verdict.
- **The test quantity:** the **sign** of the contrast, paired over the eleven
  data sets, two-sided binomial test against p = 0.5.
- **The confounders of §7** at every data set (§3.3), plus the set-size check
  from §6.

### 4.1 The power ladder, n = 11, two-sided

| same sign | p |
|---|---|
| 11/11 | **0.0010** |
| 10/11 | **0.0117** |
| 9/11 | 0.0654 |
| 8/11 | 0.227 |

> **The decision boundary is 10 of 11.** At 9/11 the result is **undecided**
> and is reported as such — not as a trend and not as "tending towards". At
> n = 10 (the version without `ARSB`) the boundary is correspondingly 10/10
> (p 0.0020) or 9/10 (p 0.0215); 8/10 (p 0.109) is undecided.

### 4.2 Rules of interpretation

- **A counts as established** if the contrast has the same sign in **at least
  10 of 11** data sets **and** the gate of §3.5 is passed.
- **A counts as refuted** if at least 10 of 11 show the **opposite** sign.
- In between: **undecided** → **S4-AB2**.
- **Without a passed gate, A is not interpreted**, whatever the count. The
  numbers are then reported and explicitly marked as not load-bearing.

### 4.3 Three counting versions, all three reported

1. **all eleven** data sets;
2. **the ten unread** ones (without `ARSB`) — decisive in case of a divergence
   (§1.1);
3. **only the seven with an intact null** (without `FN1` C231W, `FN1` C123R,
   `RB1`, `RNF4`).

In addition **always**: the number of data sets with **|z_corrected| > 2** (not
only the sign), and MDE80 with every number.

### 4.4 "No difference" only as a TOST

If A comes out undecided and one wants to say that distal and biosynthetic
behave alike, that is carried as a **TOST**. **The equivalence bound, fixed in
advance: delta = 2.0 in units of the z_corrected of the contrast.** The reason:
2.0 is the same threshold with which the project separates "excursion" from
"flat" throughout. A difference smaller than what the project treats as noise
everywhere else is practically zero. A null result that is not rejected is
**not** a finding.

---

## 5. Question B — design, and the danger of circularity

**The problem, named in advance:** XBP1 targets *are* ER genes. A regulon of
these factors overlaps with `S_DISTAL` by its nature. Without a hard rule, B
measures itself.

### 5.1 The target set — external and version-fixed

- **The primary source: CollecTRI**, obtained through
  `decoupleR::get_collectri(organism = "human", split_complexes = FALSE)`,
  `decoupleR` 2.12.0, the `OmnipathR` cache, **retrieved 2026-08-18**. Name,
  version and retrieval date stand in the script header **and** in every CSV
  produced, as with PanelApp 309.
- **The justification for the choice** (made after the sizes in §1.2 were
  known, and therefore disclosed here): CollecTRI is the curated successor
  resource, covers all four factors, and with 126 edges gives the size gate a
  chance at all. DoRothEA A to C has **zero** targets for MIA3 and would be too
  small for three of the four factors on its own.
- **The secondary source, as a sensitivity check only:** DoRothEA A to C
  (`dorothea` 1.18.0). It **cannot** change the verdict; it is reported.
- **`S_TFZIEL`** = the union of the target genes of the four factors XBP1,
  ATF6, CREB3L1 and MIA3, excluding the four factors themselves, mapped onto
  Ensembl identifiers through `org.Hs.eg.db` (the same mapping as `GOr()`).
- **No gene list by hand.** A list written down from memory is not an external
  set.

### 5.2 The disjointness rule

`S_TFZIEL` is made disjoint from `S_MASCHINE`, `S_FRACHT0`, `S_DISTAL`,
`S_BIOSYN` and `S_ZYKLUS` — the same treatment as the neutral sets. The
subtraction happens **before** any computation and is logged with the sizes
before and after.

### 5.3 The overlap gate, in advance

If fewer than **30 genes** remain in the pool after the subtraction, or if that
is the case at **fewer than 8 of the 11** data sets, B is **not testable** →
**S4-AB3**. That is then reported and B is dropped. This is an admissible
outcome. **No** substitute set is contrived and no confidence level is relaxed
afterwards in order to pass the gate.

### 5.4 The statistics, both fixed in advance, none chosen afterwards

- **B-I:** `kontrast_f(iv, ZIEH[["20"]], S_TFZIEL, S_BIOSYN, ...)`
- **B-II:** `kontrast_f(iv, ZIEH[["20"]], S_TFZIEL, S_FRACHT0, ...)`

That gives the question "like the apparatus or like the cargo" an answer at
all: if the transcription-factor arm behaves like the apparatus, B-I is flat
and B-II shows an excursion; if it behaves like the cargo, the reverse. The
counting and the power ladder are as in §4.1, the three versions as in §4.3.

**An outcome in which both contrasts show an excursion, or both are flat, is
admissible and is reported as "neither"**, not counted towards either side.

### 5.5 The four factors themselves

XBP1, ATF6, CREB3L1 and MIA3 are reported per data set as **individual genes,
descriptively**, with `iv` and `dWT`. **Four genes are not a set; no test is
made from them.**

---

## 6. What would refute S4

Stated explicitly, so that it is not negotiable afterwards. The finding on A
does **not** carry if any one of these occurs:

1. **The contrast depends on gene length** — it turns under the length-matched
   null.
2. **It turns under T2 residualisation.**
3. **It follows set size.** `S_DISTAL` at 545 against `S_BIOSYN` at 2 711 genes
   is a real difference. **The check, fixed in advance:** at every data set the
   contrast is additionally computed with `S_BIOSYN` **reduced** to the size of
   `S_DISTAL` — 20 random subsets of `S_BIOSYN` at the pool size of `S_DISTAL`,
   seed 20260818, reported as the median and range of z_corrected. If the sign
   turns at more than one data set, the finding is size-driven.
4. **It appears in only one of the two arms** — the count is therefore reported
   separately for osteogenic (4 data sets) and chondrogenic (7 data sets)
   in addition. **Without comparing magnitudes across arms** (§4, the arm
   rule).
5. **It disappears without the cell-cycle genes.**
6. **The neutral contrast shows an excursion too** (the gate, §3.5).

---

## 7. Stopping criteria for S4

- **S4-AB1** — the gate of §3.5 is not passed. A then does not carry either.
- **S4-AB2** — A is undecided (9/11 or fewer).
- **S4-AB3** — B is not testable (overlap, set size, or no database).
- **S4-AB4** — the reference implementation turns out to be faulty.

**In every one of these outcomes, step E — the report — is carried out
nonetheless.**

---

## 8. The duty to report

Every folder receives a log **and** a CSV as well as a protocol. The report
contains every number with its MDE80, all three counting versions, the power
ladder of §4.1 in the text, and **the counter-example `ARSB` prominently, not
in a footnote**. In case of a stop: which criterion, with which number, and
what follows from it for the metric.

The word **"specific"** is not used for the scissors.

---

## 9. Release

Dated and in force **2026-08-18**, before the first call of `kontrast_f` with
`S_DISTAL`, `S_BIOSYN` or `S_TFZIEL` and before the first reading of the
columns `nur_distal` and `nur_biosyn`.

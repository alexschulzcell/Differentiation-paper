> Translated from the German original of 2026-08-22. The content, the dates
> and every number are unchanged.

# Preregistration M-D — the donor-resolved scissors

**Dated 2026-08-22, before the first download** of a phase D cohort and before
the first statistic. It replaces phase M-B as the main computation;
`PREREG_M_B.md` remains valid unchanged for what was preregistered there.

What was **known** at the time of this preregistration is stated explicitly, so
that nothing is presented as blind that is not:

- every number of phases M-A, M-B and M-C, including the leave-one-out table in
  `PROTOCOL_M_B.md` addendum 1;
- the existence and rough design of the candidates named in §5 from the screen
  of 2026-08-21 (entity, approximate sample count, presence of an
  undifferentiated arm) — but **not** any module read-out on them. On none of
  the candidate cohorts has `dWT`, `iv`, a concordance or a correlation ever
  been computed. The single exception is **GSE247491/GSE247528** (SERPINA3),
  which has already been analysed study-wise as Figure S3C — **never
  donor-resolved**; it is dealt with separately in §11.

---

## 1. Question and prediction

**The question.** When a donor runs their differentiation programme while
carrying a lesion: is the programme the same across donors and the lesion
response not?

**H1 (the programme is shared).** The per-donor `dWT` vectors agree across
donors on the 173 module genes more strongly than a baseline-stratified
background draw would lead one to expect.

**H2 (the lesion response is individual).** The per-donor `iv` vectors do
**not** do so on the size-matched lesion set.

**H3 (the scissors).** H1 and H2 hold simultaneously, and the difference
between them reaches the threshold in §9.

This is **not a new axis**. It is the fixed 173-gene statistic of the main part
applied to the unit it was meant for. `PREREG_S6.md` §1 remains binding: no new
gene set, no readjustment of the module, no search for a better axis.

---

## 2. The unit, and what a cell is

**The unit is the donor (or patient).** Clones of one line are **one** donor.
Wells, passages, technical and biological replicates never count.

A **cell** of phase D is a triple

> donor x differentiation axis (osteogenic / chondrogenic) x study

with a **complete 2 x 2**:

|  | undifferentiated / day 0 | differentiated |
|---|---|---|
| **control / corrected** | yes | yes |
| **lesion / patient** | yes | yes |

One donor can contribute several cells (osteogenic and chondrogenic). They
count as separate cells, **not** as separate donors; the null treats them
together, as §8 says.

---

## 3. The two quantities per cell

Within **one** study, at symbol level, on a log2 scale (counts: CPM,
`log2(x+1)`; arrays: the deposited normalised value), z-standardised per gene
over the included samples of **the same** study — identical to
`reference_implementations/manuscript/methods/03_metric.R`, function `kern`:

    dWT_p = z(control_p, differentiated) − z(control_p, undifferentiated)
    iv_p  = [z(lesion_p, differentiated) − z(lesion_p, undifferentiated)]
            − [z(control_p, differentiated) − z(control_p, undifferentiated)]

Where several samples exist per arm, they are **averaged** before the
difference is formed. **Nothing** is matched to a covariate of the baseline
(the "not a flat null" guard).

The **lesion set** is the size-matched 173-gene set of highest `iv` consistency
from `derived_data/M_patienten/laesionssatz_173.csv`, deterministically sorted,
as in phase M-B. It is not rebuilt.

---

## 4. Inclusion and exclusion rules

- **E1** human; a complete 2 x 2 within **one** study. No assembling from two
  studies.
- **E2** the lesion is a **diagnosis or a patient mutation**, not a pure
  engineering intervention in a healthy line. **Isogenic correction of a
  patient line is explicitly allowed and preferred.**
- **E3** **donor resolution**: the assignment sample → donor is documentable
  from the GSM metadata.
- **E4** at least 60 of the 173 module genes measurable.
- **E5** **no minimum n per study.** The minimum applies at the level of the
  synthesis: **at least 6 cells from at least 3 independent studies.** (A
  deliberate correction to M-B, where "n ≥ 5 patients per cohort" excluded
  exactly the data sets that carry the right contrast.)

**Exclusion codes:** `A1` batch confounded with the lesion, or assembled from
two studies · `A2` no control group · `A3` no gene-level expression matrix ·
`A7` not human material · `A9` single-cell format · `M2` immortalised line ·
`M3` perturbation unrelated to a patient lesion · `M4` no skeletal axis · `M6`
calibration not passed · **`D1` no donor resolution despite a 2 x 2** (new) ·
**`L1` a sample treated with LiCl or lithium** (§10).

---

## 5. The candidate pool — named in advance

Every candidate is **verified against `data_raw/_meta/<GSE>_proben.csv`**, not
taken from memory. Named in advance:

| accession | expected design | previous exclusion reason |
|---|---|---|
| GSE218101 | MPS VI, 4 patient lines x empty vector / gene-corrected x day 0 / day 14 | `M1` (n = 4) — void under E5 |
| GSE221128 | FOP, iMSC, FOP / resFOP x day 0 / day 6 | `M1` (n = 3) — void under E5 |
| GSE247491 / GSE247528 | SERPINA3, 3 donors, chondrogenic and osteogenic, days 0/3/7 | already used as Figure S3C, never donor-resolved |
| GSE244375 | Gly1170Ser, WT / het / hom x early / late, iPSC cartilage | `M3` (clones) — as an isogenic lesion series, clones = 1 donor |
| GSE148728 | COL10A1 / MATN3 mutants plus isogenic controls, iPSC cartilage | `M3` — likewise |
| LAMA5-USC (our own data) | WT / KO x undifferentiated / chondrogenic / osteogenic | already two of the eighteen; **one** donor (§10) |

**A new search**, logged with its date, search string and hit count. The search
is for the **design**, not for entities — that was the gap in phase M-B:

- "patient-derived" together with "day 0" or "undifferentiated"
- "isogenic control" together with "differentiation"
- "gene-corrected" together with "iPSC"

each crossed with osteogenic or chondrogenic differentiation, organism Homo
sapiens, series format. The search is logged **completely**, including the hits
that are excluded.

**Double use.** GSE218101 and GSE221128 are already among the eighteen. They
may additionally be carried as phase D cells **if** the double use is named in
both figures and is **not** presented as independent confirmation. This is
fixed in advance and justified in the protocol.

---

## 6. Calibration — built in, not appended

**The positive control is `dWT` itself.** A cell carries a finding only if its
**own** differentiation contrast finds the canonical lineage markers of its
axis:

    _module.kontrast(dWT_p, set_a, all other marker sets)
    set_a = OSTEOGEN or CHONDROGEN, directional, p < 0.05, contrast > 0

`NAIV` is **not admissible** as the calibration of a module test
(`PROTOCOL_M_B.md` addendum 2). The marker sets in `_marker.py` are **not**
changed.

If a cell fails: code `M6`, a log entry, **no number reported**. The
calibration runs **before** the main computation.

---

## 7. The statistic ladder — fixed in advance, all three reported

The order stands here and is **not** changed after looking at the numbers. The
best one is **not** selected.

### S1 (primary) — between-donor correlation, continuous

Per gene set and per pair of cells (*p*, *q*) **from different donors**: the
Spearman rho of the vectors restricted to the set genes. The statistic is the
**mean rho over all pairs**. Computed separately for `dWT` (programme) and `iv`
(lesion response).

### S2 — variance decomposition: shared against individual

Per gene set: the share of variance explained by the **first principal
component** of the donor-by-gene matrix. The statistic is the **share of shared
variance**.

### S3 — sign concordance, in both readings

(a) **directional** against `ri` — the statistic from M-B, for continuity: the
mean share of donors with `sign = ri`.
(b) **direction-free** — the mean degree of agreement,
`mean_g |mean_p sign(x_pg)|`.

---

## 8. Null model, seed, detection limit, aggregation

**The null (the same for all three statistics).**
A baseline-stratified background draw as in `_module.py`: gene sets of equal
size from the background measurable in that cell, drawn from **the same
deciles of baseline expression** (the undifferentiated control arm), carrying
**the same set of signs** as the set under test. The statistic is computed on
the drawn set exactly as on the real one.

**20 000 draws, seed 20260822.**

**A donor-flip null** for the synthesis, as in Figure S3B. **Cells of the same
donor are permuted together**, never separately.

**No matching to a covariate of the baseline.**

**Detection limit.** Per cell and per statistic

    MDE80 = null mean + 2.8 x null SD .

**Without this number no finding of this phase is reported** — neither a
positive one nor a negative one.

**Aggregation. Only the study synthesis of Figure S3B, never an arithmetic
mean.** Always together with:

- a complete **leave-one-out computation** (each cell removed once),
- a separate analysis **with and without** cells having fewer than 3 control
  samples,
- a statement of how many cells reach their **own** MDE80.

**Multiple comparison.** Bonferroni over the three statistics,
alpha = 0.0167 for the primary statement. S1 is named as primary in advance.

---

## 9. What "found" means — the thresholds, in advance

| level | condition | status |
|---|---|---|
| **strong** | S1(`dWT`) above the null (synthesis p < 0.01) **and** S1(`iv`) not **and** at least 3 individual cells above their own MDE80 | confirmatory, carries the main title |
| **medium** | S1(`dWT`) above the null, S1(`iv`) not, but no individual cell above MDE80 | the direction reproduces, not confirmed at threshold |
| **weak** | only S2 or S3 separates, S1 does not | exploratory, supplement |
| **not found** | none of the three statistics separates `dWT` from `iv` | a null result with a measured limit |

**Stopping criterion.** If the screen does **not** reach at least 6 cells from
at least 3 independent studies with a passed calibration, the main computation
is **not reported as confirmatory**, whatever the numbers turn out to be; it
then appears as an exploratory computation with a measured limit.

**There is no searching on until it fits.** The search axes in §5 are final. No
threshold is lowered after the numbers are known; the primary statistic is not
switched.

---

## 10. Fixed exclusions

- **No sample treated with LiCl or lithium enters any computation** — not a
  positive control, not a background, not a figure. The sample list of every
  cell is checked against `LiCl` and `lithium` before the computation, **with a
  note in the protocol even where nothing was found** (code `L1`). A LiCl arm
  is no substitute for an undifferentiated arm.
- Our own COL2A1 cohort and its single-cell data are closed questions and are
  not used.
- **COL2A1 and LAMA5 are not assembled into one 2 x 2** (E1).
- Clones are not donors. The LAMA5-USC series (WT1-3, KO9/46/75) is **one**
  donor and serves as an **isogenic lesion series**, not as several donors.
- The word **"specific"** is not used for the scissors.

---

## 11. GSE247491 / GSE247528 — the one cohort already touched

These series have been analysed study-wise as Figure S3C. **Per-donor `dWT`
and `iv` vectors have never been formed from them**, and none of the three
statistics in §7 has ever been computed on them. They are therefore carried as
regular phase D cells, but marked in the protocol and in the figure as
**already touched**, and the leave-one-out computation reports their
contribution separately. If the synthesis carries **only** with them, the
result is **not** confirmatory.

---

## 12. Output

`derived_data/M_donoren/`: `suchlauf.csv`, `zellen_sichtung.csv`,
`eichung.csv`, `statistik.csv`, `nullen.csv`, `synthese.csv`,
`auslassung.csv`, `54_log.txt`.
Protocol: `preregistrations/PROTOCOL_M_D.md`.

New statistics are **added to `reference_implementations/_module.py`** and
documented there, not built beside it. A self-test against the known null rate
runs before the first real computation.

---

## Addendum 1 — 2026-08-22, before the first statistic

Four decisions that became necessary while writing the computation down. They
are dated **before any statistic** of this phase; what was known at that point
were only the **screening results** (sample titles, characteristic fields,
study designs and the search log
`derived_data/M_donoren/suchlauf.csv`, search of 2026-08-22), **no** `dWT`,
**no** `iv`, **no** correlation and **no** calibration number.

### (a) S2 is **not** centred over the donors

§7 had said "genes centred, donors as rows". That is the wrong way round:
centring the gene columns over the donors removes exactly the **shared**
component that S2 is meant to measure, and the statistic could then no longer
see the shared programme component.

**Decision:** S2 is the share of the first principal component in the **total
sum of squares** of the uncentred donor-by-gene matrix. The vectors are already
differences on a z scale and are therefore formed around a meaningful zero
point. The null draws background sets and computes the same uncentred
decomposition, so the comparison is fair.

### (b) E2 applies to the **lesion arm**, not to the differentiation arm

A cell whose "lesion" is a pure engineering intervention in healthy cells (an
siRNA knockdown, a knockout of a healthy line) does **not** satisfy E2. Its
`dWT` is untouched by that, however: the control arm is an ordinary
differentiation of an ordinary donor.

**Decision:** such cells are carried as **`dWT`-bearing but not `iv`-bearing**
(marked `E2-partial`). Their `dWT` enters the programme computation, their `iv`
does **not** enter the lesion computation; it is reported separately as an
**engineering response** and is never mixed with the lesion response in one
number. After the screen this affects the SERPINA3 siRNA pair (GSE247491,
GSE247528) and our own LAMA5-USC series.

The reason for allowing this at all: H1 ("the programme is shared") is a
statement about **donors**, not about lesions. Excluding a cell whose
differentiation arm is sound would weaken H1 artificially. The separation is
carried in every figure and every table, and the leave-one-out computation
reports the contribution of these cells separately.

### (c) Pairs: across studies as primary, within study as a sensitivity

S1 asks about the agreement **between donors**. Restricting the pairs to donors
of the same study would mean that studies with only **one** donor (GSE221128,
GSE245585, LAMA5-USC) could contribute nothing at all, and the statistic would
come down to two studies.

**Decision:** primary are **all pairs of different donors**, across studies as
well. That is the **conservative** direction: pairs spanning studies carry
additional platform, protocol and laboratory noise and can only lower rho. The
null draws its background sets over the same pairs, so that batch structure hits
the module and the null alike — the same argument as in `_module.konkordanz`.

The within-study version is computed and reported in addition (as a
sensitivity). Pairs of the same donor are excluded in both versions; cells of
the same donor on two axes are never correlated against each other.

### (d) Which statistic gets an MDE80 per cell

**S1** and **S3a** are defined per cell (the mean over the partners of that
cell, and the sign agreement of that cell respectively) and receive a null, a z
and an MDE80 per cell — those are the numbers against which §9 measures "at
least 3 individual cells above their own MDE80".

**S2** and **S3b** are statistics of a **set** of donors; they receive an MDE80
at the level of the study or of the pool, not per cell. For studies with only
one donor they are undefined and are reported as such (`—`), not filled with a
substitute value.

### (e) What these four points do not change

The primary statistic remains **S1**. The thresholds in §9, the seed
(20260822), the 20 000 draws, the baseline-stratified null, the donor-flip
synthesis, the Bonferroni correction and the stopping criterion stand word for
word.

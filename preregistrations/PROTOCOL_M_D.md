> Translated from the German original of 2026-08-22. The content, the dates
> and every number are unchanged.

# Protocol M-D — the scissors sought donor-resolved

Computed **2026-08-22**, following `PREREG_M_D.md` (dated before the first
download) and its addendum 1 (dated before the first figure).
Scripts: `data_acquisition/20_donor_search.py`, `30_donor_statistics_self_test.py`,
`31_donor_cells_build_calibrate.py`, `32_donor_statistic_ladder.py`, `33_donor_circularity_control.py`,
`21_donor_manual_screen.py`. Seed 20260822, 20 000 draws. All numbers in
`derived_data/M_donoren/`.

**The result in three sentences, up front:**

1. The differentiation programme runs **in the same direction** between donors —
   S1 (primary) 0.349 against a null of 0.273 ± 0.025, **z +3.00, p 0.0028**,
   above the Bonferroni threshold, robust under every omission, and **stronger**
   (z +4.51) when the cells that helped define the module are removed.
2. The **lesion response cannot be tested**: of the cells with a real patient
   lesion, exactly **one** passes the built-in calibration. One cell allows no
   between-donor metric. **No number** is reported for it — neither a positive
   nor a negative one.
3. The scissors have therefore **not been measured** in a donor-resolved way.
   The stopping criterion of the preregistration (§9) applies; the phase is
   **exploratory**, not confirmatory. The fallback core of §7.1 of the plan
   carries the main figures.

---

## 1. Order and dating — expressly

The preregistration was written **before** the search and **before** any cell
was built. The dates in the files diverge, and that deserves naming:
`suchlauf.csv` carries the system time **2026-08-21**, the preregistration the
date **2026-08-22**. Both are the same session, in the night from the 21st to
the 22nd; the preregistration came first. The order can be reconstructed from
the files and is recorded here as it is, instead of being smoothed over.

Also up front: **nothing was downloaded.** All six study units already lay in
the project or in the backups of its predecessor. The search only fetched
accessions and metadata.

---

## 2. A self-test of the new statistics — before the first real computation

`30_donor_statistics_self_test.py`, `selbsttest.csv`. Cell vectors drawn from the background,
that is **without** any shared structure:

| metric | mean z | sd | share z > 2 | share \|z\| > 2 |
|---|---|---|---|---|
| S1 | +0.05 | 1.03 | 3.5 % | 5.0 % |
| S2 | −0.05 | 0.93 | 3.5 % | 4.5 % |
| S3a | +0.16 | 0.98 | 2.5 % | 3.0 % |
| S3b | +0.07 | 1.07 | 3.5 % | 5.0 % |

**The null is calibrated** (nominally 2.3 % and 5 %). With a shared programme of
strength 0.35 z underlaid — exactly the detection limit of the main part — S1
(z +4.44, 100 %), S3a (z +8.73, 100 %) and S3b (z +3.82, 97.5 %) respond;
**S2 is the least sensitive** of the four (z +0.84, 22.5 %). That is measured
before the first real number and is carried along when S2 is read.

---

## 3. Search and screening

### 3.1 The search by design

`20_donor_search.py`, `suchlauf.csv`, `treffer_roh.csv`. Database GEO `gds`,
filter `"Homo sapiens"[Organism] AND "gse"[Filter]`, crossed with osteogenic and
chondrogenic differentiation, retmax 1000 (no axis ran against the upper bound).

| axis | search term (core) | hits |
|---|---|---|
| `ANLAGE_NAIV` | "patient-derived" plus "day 0" / "undifferentiated" / "baseline" | 17 |
| `ANLAGE_ISOGEN` | "isogenic control/pair/line" plus "differentiation" | 7 |
| `ANLAGE_KORRIGIERT` | "gene-corrected" plus iPSC | 7 |
| **unique series** | | **22** |

The search finds GSE218101 on **two** axes — it is thereby checked against a
known positive example, without that example having defined the search.

**The real finding of the search is how small the number is.** Phase M-C had
screened 1424 series along entities and found 92 % without an undifferentiated arm. If one
searches instead by **design**, **22** series remain in the whole database, and
of those **18** are muscle, tumour or developmental studies without a skeletal
differentiation axis. The bottleneck is not the search strategy. The bottleneck
is the data.

### 3.2 Screening by hand

`21_donor_manual_screen.py`, `sichtung_hand.csv`. All 22 hits plus the six
candidates named in advance, each checked against
`data_raw/_meta/<GSE>_proben.csv` or against the GSM metadata — not against the
series title.

| code | reason | count |
|---|---|---|
| `M4` | no osteogenic or chondrogenic axis (muscle, tumour, development) | 13 |
| `M3` | perturbation unrelated to a patient lesion | 3 |
| `A1` | **no undifferentiated arm** | 2 |
| `A2` | no control group in the 2 x 2 sense | 1 |
| `A9` | single-cell format | 1 |
| **IN** | | **7 study units** |

Notable, because they nearly carried:

| GSE | reason |
|---|---|
| **GSE244375** | Gly1170Ser, iPSC cartilage. **No undifferentiated arm**: day 34 against day 44, both differentiated. The `M3` exclusion of phase M-B was the wrong reason — the right one is `A1`. In addition WT, het and hom are clones of **one** line |
| **GSE148728** | COL10A1/MATN3 with isogenic controls, but **only differentiated cartilage pellets**. "Biological replicates" are not donors |
| **GSE183525** | FOP patient iPSC — but monocytes and an inflammation signature, no differentiation axis |

### 3.3 What was included

| study | design | cells | axis | E2 |
|---|---|---|---|---|
| **GSE218101** | MPS VI (ARSB), 4 patient lines x empty vector / gene-corrected x iPS / day 14 | 4 | chondrogenic | **yes** |
| **GSE221128** | FOP (ACVR1), iMSC, FOP / resFOP x day 0 / day 6 | 1 | chondrogenic | **yes** |
| **GSE245585** | RB1 +/−, the patient's own MSC, day 0 / 7 / 14 / 21 | 1 | osteogenic | **yes** |
| GSE247491 | SERPINA3 siRNA, 3 donors, day 0 / 3 / 7 | 3 | chondrogenic | no |
| GSE247528 | SERPINA3 siRNA, 3 donors, day 0 / 3 / 7 | 3 | osteogenic | no |
| LAMA5-USC | LAMA5 KO, our own line, undifferentiated / chondrogenic / osteogenic | 2 | both | no |
| | | **14 cells** | | |

**The donor count, strictly by rule 0.4:** GSE221128 delivers **one** cell, not
three — `ex1` to `ex3` are experiments on the same line. The LAMA5-USC series
delivers **one** donor, not six — WT1-3 and KO9/46/75 are clones of one line; it
is an **isogenic lesion series**. That is stated expressly here, because it is
tempting to do otherwise.

**E2 (lesion = diagnosis or patient mutation)** is met by GSE218101, GSE221128
and GSE245585. SERPINA3 (siRNA in healthy cells) and LAMA5 KO (a knockout in a
healthy line) do **not** meet it; under addendum 1 (b) they are
**`dWT`-carrying, `iv`-non-carrying**.

### 3.4 Double use — disclosed

**Four of the six study units belong to the eighteen points from which the
173-gene module was formed**: GSE218101, GSE221128, GSE245585 and the LAMA5-USC
series (two points). The SERPINA3 series took no part in forming the module but
has already been analysed study-wise as Figure S3C.

**Phase D therefore contains not a single untouched study.** What is new is not
the material but the **question**: donor-resolved, between donors, with a
built-in calibration per cell. None of it may be presented as an independent
confirmation of the eighteen points. The circularity control in §6 draws the
consequence.

### 3.5 LiCl

`licl_pruefung.csv`. All six study units, **115 samples** in total, checked
against `LiCl` and `lithium`: **no hit**. The note stands here also because the
preregistration expressly demands it **even when** nothing was found.

---

## 4. The built-in calibration — and what it costs

`31_donor_cells_build_calibrate.py`, `eichung.csv`. The positive control is `dWT` itself: does the
cell's **own** differentiation contrast recover the canonical lineage markers of
**its** axis (`_module.kontrast`, directed, p < 0.05)?

| cell | axis | set A | contrast | z | p | verdict |
|---|---|---|---|---|---|---|
| GSE218101 Line #1 | chondrogenic | `CHONDROGEN` | −0.539 | −1.27 | 1.00 | failed |
| GSE218101 Line #2 | chondrogenic | `CHONDROGEN` | +0.484 | +1.03 | 0.302 | failed |
| GSE218101 Line #3 | chondrogenic | `CHONDROGEN` | +0.160 | +0.37 | 0.721 | failed |
| **GSE218101 Line #4** | chondrogenic | `CHONDROGEN` | +1.264 | +2.33 | **0.016** | **passed** |
| GSE221128 FOP | chondrogenic | `CHONDROGEN` | +0.606 | +1.90 | 0.057 | failed |
| GSE245585 WT1 | osteogenic | `OSTEOGEN` | −0.408 | −1.03 | 1.00 | failed |
| **SERPINA3 D1** | chondrogenic | `CHONDROGEN` | +1.412 | +3.57 | 0.0003 | **passed** |
| **SERPINA3 D2** | chondrogenic | `CHONDROGEN` | +2.234 | +5.60 | 0.0001 | **passed** |
| **SERPINA3 D3** | chondrogenic | `CHONDROGEN` | +1.711 | +4.19 | 0.0001 | **passed** |
| **SERPINA3 D1** | osteogenic | `OSTEOGEN` | +0.787 | +2.12 | 0.035 | **passed** |
| **SERPINA3 D2** | osteogenic | `OSTEOGEN` | +0.683 | +2.01 | 0.040 | **passed** |
| SERPINA3 D3 | osteogenic | `OSTEOGEN` | −0.125 | −0.37 | 1.00 | failed |
| **LAMA5-USC** | chondrogenic | `CHONDROGEN` | +1.765 | +2.67 | 0.005 | **passed** |
| LAMA5-USC | osteogenic | `OSTEOGEN` | −0.100 | −0.21 | 1.00 | failed |

**7 of 14 cells pass. Of the cells with a real patient lesion, exactly one
passes.**

That is the most expensive sentence of this protocol, and it is the most
important. The calibration does not ask whether the cell is the right tissue —
that is what phase M-B asked, with `NAIV`, and that is why it was wrong. It asks
whether the cell **runs the axis on which the module is defined**. And for half
of the published differentiation experiments the answer is: not measurable
against the textbook markers of its own lineage.

Two examples, so that this does not stay abstract:

- **GSE245585** (osteogenic, day 0 → 21): `RUNX2` +1.61 and `ALPL` +1.50 go up,
  `SOX9` −0.23 and `ACAN` −1.69 go down — but the set as a whole does **not**
  move against the other lineages, at −0.41. The cell differentiates
  recognisably, but not recognisably **osteogenically** against the
  alternatives.
- **LAMA5-USC osteogenic**: the osteogenic set rises (+0.35), and the
  chondrogenic set rises **more** in the same sample (+1.30). Against the other
  lineages no directed contrast remains. The same series passes the chondrogenic
  calibration easily (z +2.67).

**That applies to our own data as well.** There is no special treatment: the
osteogenic arm of our own LAMA5 series fails and carries no number.

---

## 5. The ladder of statistics

`32_donor_statistic_ladder.py`, `statistik.csv`, `je_zelle.csv`. Calibrated cells only, a
baseline-stratified null, 20 000 draws, seed 20260822.

### 5.1 The programme (`dWT`, 173 module genes) — 7 cells, 5 donors, 4 studies, 19 donor pairs

| statistic | observed | null | z | p | MDE80 | above MDE80 |
|---|---|---|---|---|---|---|
| **S1 (primary)** | **0.3489** | 0.2731 ± 0.0252 | **+3.00** | **0.0028** | 0.3438 | **yes** |
| S2 | 0.4596 | 0.3843 ± 0.0203 | +3.72 | 0.0001 | 0.4410 | yes |
| S3a | 0.7350 | 0.5151 ± 0.0207 | +10.61 | 0.0001 | 0.5731 | yes |
| S3b | 0.5723 | 0.4850 ± 0.0209 | +4.19 | 0.0001 | 0.5434 | yes |

**All four metrics lie above their null and above their own detection limit**,
and S1 also above the Bonferroni threshold α = 0.0167.

**Per cell**, S1 against the **cell's own** MDE80:

| cell | S1 | null | z | MDE80 | |
|---|---|---|---|---|---|
| GSE218101 Line #4 | +0.150 | 0.173 ± 0.041 | −0.56 | 0.288 | below |
| SERPINA3 D1 chondro | +0.335 | 0.269 ± 0.040 | +1.66 | 0.381 | below |
| SERPINA3 D2 chondro | +0.434 | 0.355 ± 0.036 | +2.23 | 0.455 | below |
| **SERPINA3 D3 chondro** | +0.384 | 0.264 ± 0.040 | +2.98 | 0.377 | **above** |
| **SERPINA3 D1 osteo** | +0.457 | 0.340 ± 0.035 | +3.29 | 0.439 | **above** |
| SERPINA3 D2 osteo | +0.375 | 0.279 ± 0.038 | +2.55 | 0.384 | below |
| LAMA5-USC chondro | +0.340 | 0.257 ± 0.036 | +2.33 | 0.357 | below |

**2 of 7** cells reach their own S1 MDE80; for S3a it is **7 of 7** (z +3.21 to
+10.26) — but S3a is precisely the metric affected by the double use (§6).

**Leave-one-out** (`auslassung.csv`, 2 000 draws per run):

| without | n | S1 | z | p |
|---|---|---|---|---|
| — (all) | 7 | +0.349 | +2.98 | 0.0040 |
| GSE218101 Line #4 | 6 | +0.441 | **+3.82** | 0.0010 |
| SERPINA3 D1 chondro | 6 | +0.354 | +2.95 | 0.0030 |
| SERPINA3 D2 chondro | 6 | +0.318 | +2.87 | 0.0070 |
| SERPINA3 D3 chondro | 6 | +0.333 | +2.14 | 0.0360 |
| SERPINA3 D1 osteo | 6 | +0.310 | +2.18 | 0.0290 |
| SERPINA3 D2 osteo | 6 | +0.339 | +2.47 | 0.0080 |
| LAMA5-USC chondro | 6 | +0.353 | +2.43 | 0.0150 |

**No single cell carries the result**, and none reverses it. The range is z
+2.14 to +3.82. There are no cells with fewer than three control samples — the
separate analysis is therefore not applicable and is noted as such.

### 5.2 The lesion response (`iv`) — **no number**

Of the three studies with a real patient lesion, **one cell** passes the
calibration (GSE218101 Line #4). One cell has no donor partner; S1, S2 and S3b
are not defined.

**No value is reported for the lesion response.** Not z, not p, not "lies on the
null" — nothing. That is not a formality: a null finding without a detection
limit is not a finding in this work, and a detection limit needs at least two
donors.

### 5.3 The engineering response (`iv` of the cells without E2) — reported separately

6 cells, 4 donors, 13 pairs. **Never mixed with the lesion response.**

| statistic | observed | null | z | p | MDE80 | above MDE80 |
|---|---|---|---|---|---|---|
| S1 | 0.0887 | 0.0536 ± 0.0243 | +1.44 | 0.153 | 0.1218 | no |
| S2 | 0.3198 | 0.2639 ± 0.0175 | +3.20 | 0.0041 | 0.3129 | yes |
| S3a | 0.6163 | 0.4889 ± 0.0172 | +7.42 | 0.0001 | 0.5370 | yes |
| S3b | 0.3808 | 0.3484 ± 0.0210 | +1.54 | 0.130 | 0.4073 | no |

On the **primary** statistic the engineering response does **not** lie above its
null (z +1.44, p 0.15), and its S1, at 0.089 against 0.349 for the programme, is
smaller by a factor of four. That is the direction the scissors predict — **but
it is not the scissors.** An siRNA and a laboratory knockout are not diagnoses.
The sentence this work is looking for is about patients.

### 5.4 The study synthesis — and why it can decide nothing here

`synthese.csv`, `synthese_je_studie.csv`. The procedure as in Figure S3B: the
metric study-wise against a joint **donor-flip null** in which the sign of
**whole donors** is turned; cells of the same donor always together. 10 000
rounds.

| metric | observed | null | z | p | MDE80 | studies with a defined metric |
|---|---|---|---|---|---|---|
| S1 | +0.683 | −0.003 ± 0.412 | +1.66 | 0.249 | **+1.151** | **2 of 4** |
| S2 | +0.804 | **degenerate** | — | — | — | — |
| S3a | +0.769 | +0.501 ± 0.134 | +2.00 | 0.017 | +0.875 | 4 of 4 |
| S3b | +0.902 | +0.749 ± 0.092 | +1.67 | 0.249 | +1.006 | 4 of 4 |

**Three things deserve to be named honestly here:**

1. **On S1 the synthesis can demonstrate nothing in principle.** Its MDE80 lies
   at **+1.151** — on a scale whose maximum is 1.0. There is no value this
   synthesis could recognise as a finding. The reason is not the computation but
   the data situation: two of the four studies have only **one** donor, where a
   between-donor metric is not defined, and the flip null has a spread of 0.41 at
   two to three donors.
2. **The flip null is mathematically degenerate for S2.** If one turns row
   signs, the new matrix is D·M with D orthogonal-diagonal; the singular values,
   and therefore the PC1 share, do **not** change. The spread of the null is
   exactly zero. **No** number is therefore reported instead of an apparent one.
   The code sets a `degenerate` flag for this and returns `NaN` — the degeneracy
   is documented in `_module.synthese_flip`.
3. **The synthesis does not contradict the pooled result, it measures less.**
   The pool uses 19 donor pairs, the synthesis effectively two studies with two
   and three donors. Both stand side by side, and the weaker number is the one
   the preregistration named as the primary aggregation. **That is why the
   result is not confirmatory.**

---

## 6. The circularity control — after the fact, not preregistered

`33_donor_circularity_control.py`, `zirkularitaet.csv`. **This section arose after the
numbers were known** and is presented as what it is.

The occasion stands in §3.4: four of the six study units helped define the
module. A high directed concordance (S3a) is partly built in there. S1 and S2,
by contrast, ask about the agreement of the donors **with one another** and do
not use `ri` at all.

| subset | cells | S1 | z | S2 z | S3a z | S3b z |
|---|---|---|---|---|---|---|
| all calibrated | 7 | 0.349 | **+3.00** | +3.72 | +10.61 | +4.19 |
| **non-module-forming only** (SERPINA3) | 5 | 0.521 | **+4.51** | +5.36 | +10.29 | +5.02 |
| module-forming only | 2 | 0.309 | **+0.88** | +1.41 | **+9.90** | +0.98 |

The pattern is exactly the one predicted:

- On the **module-forming** cells S3a is high (z +9.90) — and S1 lies **on the
  null** at z +0.88. The directed concordance is partly built in there, the
  between-donor correlation is not.
- On the **non-module-forming** cells S1 is **stronger** than in the whole pool,
  at z +4.51.

**The programme result is not produced by circularity** — on the contrary, the
module-forming cells dilute it. And it confirms that S1 is rightly the primary
statistic: it is the only one of the three that does not profit from the double
use.

**The limitation stands**: the five carrying cells come from **one** donor
cohort (three donors, two axes, one laboratory, one protocol) and are moreover
an engineering intervention on the `iv` side. That is no substitute for three
independent studies.

---

## 7. The gate of the preregistration

**The stopping criterion (§9) applies.** What was required was **at least 6
cells from at least 3 independent studies with a passed calibration**; what was
reached on the lesion side is **1 cell from 1 study**.

Against the level table of the preregistration:

| level | condition | reached? |
|---|---|---|
| strong | S1(`dWT`) in the **synthesis** p < 0.01, S1(`iv`) not, at least 3 individual cells above MDE80 | **no** (synthesis p 0.249; 2 cells above MDE80; `iv` not testable) |
| middling | S1(`dWT`) above the null, `iv` not, no individual cell above MDE80 | **no**, not cleanly: `iv` is not testable, and 2 cells are above MDE80 |
| weak | only S2 or S3 separate, S1 does not | no — S1 separates in the pool very clearly |
| not found | no statistic separates | no |

**None of the four rows fits**, and that is not healed by choosing one. The
applicable finding reads:

> **H1 is confirmed exploratorily** — the differentiation programme runs in the
> same direction between donors (S1 z +3.00, p 0.0028; without the
> module-forming cells z +4.51), robust under every omission, with a measured
> detection limit per cell.
>
> **H2 and H3 cannot be tested.** The lesion response has exactly one calibrated
> cell. The scissors have **not been measured** donor-resolved — which is
> something other than "not found".

Phase D is therefore **exploratory**. It carries no headline. The fallback core
of §7.1 (the map of levels) and §7.2 (the map of detection limits) takes over
the main figures.

---

## 8. What the phase yielded in addition

Three findings that need no lesion data and that belong in the main text:

1. **Searching by design instead of by entity gives 22 series — and 2 of them
   have an undifferentiated arm with a skeletal axis and a patient lesion.** Phase M-C had
   shown the same along entities (92 % without an undifferentiated arm, 1424 series). Two
   independent search strategies, one bottleneck.
2. **Half of the published differentiation experiments do not pass their own
   positive control** (7 of 14 cells). As far as we can see, that has nowhere
   been checked in this literature — and it expressly applies to our own data as
   well.
3. **For the study synthesis the detection limit can be stated, and it is
   unreachable** (MDE80 +1.15 on a scale bounded by 1.0). A statement that would
   otherwise be reported as "no effect" is documented here as "not measurable".

---

## 9. What was expressly not done

- **No new convergence axis**, no gene set by hand, no readjustment of the
  module. `PREREG_S6.md` §1 untouched.
- **No second implementation of the metric.** The new statistics stand **in**
  `00_shared/_module.py` (`leiter`, `synthese_flip`) and use the
  same baseline-stratified null as `konkordanz`. `_marker.py` unchanged.
- **No arithmetic mean over studies** — nowhere.
- **No adjustment to a covariate of the baseline**; the control expression
  appears solely as a stratification variable.
- **No follow-up search after the first result**, no extension of the search
  axes, no lowering of a threshold, no change of the primary statistic.
- **No special treatment of our own data**: LAMA5 osteogenic fails and carries
  no number; LAMA5 counts as **one** donor.
- The questions closed in §6.2 of the plan (our own COL2A1 cohort, the COL2A1
  scRNA-seq, assembling from two studies) were **not** reopened.
- **No LiCl sample in any computation** — checked, 115 samples, no hit.
- The word **"specific"** was not used for the scissors.
- **No number** was reported for the lesion response, although one cell exists
  and although it would have been tempting to present it as "on the null".

---

## 10. Output

`derived_data/M_donoren/`: `suchlauf.csv`, `treffer_roh.csv`,
`sichtung_hand.csv`, `zellen_sichtung.csv`, `licl_pruefung.csv`,
`eichung.csv`, `selbsttest.csv`, `statistik.csv`, `je_zelle.csv`,
`synthese.csv`, `synthese_je_studie.csv`, `auslassung.csv`,
`zirkularitaet.csv`, `gene_level.csv.gz`, `zellen.pkl`, `54b_log.txt`,
`54c_log.txt`.

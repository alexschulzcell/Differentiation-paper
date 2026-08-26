> Translated from the German original of 2026-08-21, with its addenda of
> 2026-08-22. The content, the dates and every number are unchanged.

# Protocol M-B — individuality measured patient by patient

Computed **2026-08-21**, following `PREREG_M_B.md` (dated before the first
download) and its addendum 1 (dated before the first figure).
Scripts: `reference_implementations/50_cohort_search.py`, `50b_screening.py`,
`50c_check_candidates.py`, `51a_fetch.py`, `51_patient_variability.py`.
Seed 20260821, 20 000 draws, all numbers in `derived_data/M_patienten/`.

---

## 1. The question

Does the fixed 173-gene programme run in the same direction between
**patients**, while the lesion response stays individual? That is the medical
prediction from the internal finding (173 against 7.9 expected for `dWT`, 7
against 8.0 for `iv`).

---

## 2. The search — logged in full

`derived_data/M_patienten/suchlauf.csv`, search date **2026-08-21**, database
GEO `gds`, filter `"Homo sapiens"[Organism] AND "gse"[Filter]`, retmax 1000 (no
axis ran against the upper bound).

| axis | hits |
|---|---|
| osteogenesis imperfecta (COL1A1/COL1A2) | 144 |
| mucopolysaccharidoses (ARSB/IDUA/IDS) | 621 |
| pseudoachondroplasia / MED (COMP, MATN3) | 79 |
| achondroplasia / the FGFR3 spectrum | 95 |
| SHOX, Léri-Weill | 17 |
| cleidocranial / campomelic dysplasia (RUNX2, SOX9) | 381 |
| FOP (ACVR1) | 33 |
| free text "skeletal dysplasia / chondrodysplasia" | 9 |
| free text "short stature / dwarfism" | 54 |
| free text "growth plate / epiphyseal" | 43 |
| **unique series** | **1424** |

**The mechanical pre-screen** (`sichtung_mechanisch.csv`, all 1424 series with
their verdict): `A3` no gene-level expression matrix 566; `M4` not a skeletal
entity 430; `M1` fewer than 7 samples 237; `A9` single-cell format 34; `M3`
perturbation without a patient assignment 22; `M2` immortalised line 8;
**127 into the screen by hand**.

**The screen by hand**: each of the 127 accessions checked against its **GSM
metadata** (`data_raw/_meta/<GSE>_proben.csv`), not against the series title.
Notable exclusions, because they nearly carried:

| GSE | reason |
|---|---|
| GSE157587 | 10 OI patients, MSC, patient-resolved — but **no control group** (`A2`) |
| GSE244375 | COL2A1 p.Gly1170Ser — iPSC clones, not patients (`M3`) |
| GSE218101 | MPS VI, isogenic pairs — only **4** patient lines (`M1`) |
| GSE148728 | chondrodysplasia mutant iPSC, biological replicates (`M3`) |
| GSE120558 | microcephalic short stature — 4 patients, mixed platforms (`M1`) |
| GSE156466 / GSE221128 | FOP, n = 3 (`M1`) |
| GSE246390 | growth plate, 4 patients, a loading experiment (`M1`) |
| GSE16464 | OA chondrocytes, 3 against 3 donors (`M1`) |
| GSE185333 | OA knee, 4 patients, no healthy control (`A2`) |

**Seven cohorts** satisfy E1 to E5 (`kohorten_sichtung.csv`):

| GSE | entity | tissue | patients | controls | module genes |
|---|---|---|---|---|---|
| GSE186141 | osteogenesis imperfecta | primary osteoblasts | 6 | 2 | 167/173 |
| GSE22855 | enchondromatosis (Ollier) | cartilage / enchondroma | 7 | 6 | 142/173 |
| GSE292600 | acromelic dysplasia (ADAMTSL2/FBN1) | dermal fibroblasts | 8 | 3 | 171/173 |
| GSE77758 | EDS-HT / JHS | dermal fibroblasts | 5 | 6 | 171/173 |
| GSE160207 | osteogenesis imperfecta | whole blood | 7 | 5 | 153/173 |
| GSE228522 | FOP (ACVR1) | CD14⁺ monocytes | 6 | 6 | 173/173 |
| GSE58435 | Turner syndrome | amniotic fluid, cell-free mRNA | 5 | 5 | 166/173 |

In all seven the unit is the **patient**. Immortalised lines: none. For
GSE228522 only the **untreated** arm was used — the activin A arm is an
intervention, not a state.

---

## 3. Calibration — gate B

The tissue-identity control under addendum 1: `_module.kontrast` on the mean
control profile, its own marker set against all the others (`eichung.csv`).

| GSE | set A | contrast | z | p | verdict |
|---|---|---|---|---|---|
| GSE186141 | `OSTEOGEN` | +1.41 | +2.38 | 0.027 | **passed** |
| GSE22855 | `CHONDROGEN` | +1.72 | +4.77 | 0.0003 | **passed** |
| GSE292600 | `NAIV` | +3.96 | +4.62 | 0.0002 | **passed** |
| GSE77758 | `NAIV` | +1.97 | +2.73 | 0.0075 | **passed** |
| GSE160207 | — | — | — | — | **not calibratable** (`M6`) |
| GSE228522 | — | — | — | — | **not calibratable** (`M6`) |
| GSE58435 | — | — | — | — | **not calibratable** (`M6`) |

**Gate B passed: four cohorts with a passed calibration** (rule `AB1` does not
apply, and the finding is not single-cohort exploratory).

**Three cohorts cannot be calibrated**, because no canonical lineage marker set
exists in `_marker.py` for whole blood, monocytes and cell-free amniotic fluid
mRNA. Their numbers stand in `streuung.csv` with `eichung_bestanden = False` and
**carry no finding**. That is not a formalism: it is the same sentence as in
Figure S4A — a level without a calibration says nothing.

*That is at the same time a result in its own right: three of the seven patient
cohorts of skeletal disease that can be found at all measure no skeletal tissue.*

---

## 4. The result — the four calibrated cohorts

`streuung.csv`, `streuung_null.csv`. `U` is the preregistered agreement in
direction, `C` the concordance against the baseline-stratified null, and MDE80
the cohort's own detection limit.

**The programme (173 genes)**

| GSE | n | U | C | null | z | p | MDE80 |
|---|---|---|---|---|---|---|---|
| GSE186141 | 147 | 0.429 | 0.429 | 0.581 | **−3.95** | 1.00 | 0.689 |
| GSE22855 | 142 | 0.459 | 0.437 | 0.490 | −1.27 | 1.00 | 0.606 |
| GSE292600 | 159 | 0.570 | 0.560 | 0.466 | **+2.38** | **0.022** | 0.576 |
| GSE77758 | 171 | 0.506 | 0.497 | 0.494 | +0.07 | 1.00 | 0.601 |

**The lesion response (173 size-matched genes, `iv` consistency)**

| GSE | n | U | C | null | z | p | MDE80 |
|---|---|---|---|---|---|---|---|
| GSE186141 | 156 | 0.384 | 0.327 | 0.484 | **−4.14** | 1.00 | 0.590 |
| GSE22855 | 129 | 0.540 | 0.543 | 0.503 | +0.90 | 0.42 | 0.626 |
| GSE292600 | 160 | 0.498 | 0.500 | 0.512 | −0.30 | 1.00 | 0.620 |
| GSE77758 | 170 | 0.518 | 0.494 | 0.502 | −0.21 | 1.00 | 0.610 |

**A study synthesis over the four calibrated cohorts, descriptive**
(`synthese.csv`; **not an inferential test** — that would need a joint cohort
null model, the same restraint as for Figure 4C):

| set | mean z | range | mean U | above its own MDE80 |
|---|---|---|---|---|
| the programme | **−0.69** | −3.95 to +2.38 | 0.491 | **0 of 4** |
| the lesion response | **−0.94** | −4.14 to +0.90 | 0.485 | **0 of 4** |

---

## 5. What follows from it

**H1 (the programme is concordant) is not confirmed.** Of four calibrated
cohorts, **one** shows the predicted direction (GSE292600, acromelic dysplasia,
8 patients: z +2.38, p 0.022; two-set contrast z +4.56, p 0.0001), **one lies
clearly below its null** (GSE186141, OI osteoblasts: z −3.95; on the reading see
addendum 1), and two lie on the null. **Not a single cohort reaches its own
MDE80.** The mean z over the four cohorts is **−0.69** — internally, on the same
genes, the finding stood at 22 times the chance expectation.

**H2 (the lesion response is individual) is confirmed.** In no calibrated cohort
does the lesion response lie above the null; mean z −0.94, 0 of 4 above MDE80.
That is the continuation of the internal (7 against 8.0) and external
(Figure 4C) finding at the patient level.

**H3 (the scissors) is not confirmed.** The two sets do not behave differently
in patient cohorts: −0.69 against −0.94, both on or below the null. The
scissors that separate eighteen perturbation experiments **cannot be
recovered** between patients.

**This is core 2 of the plan (§5.2)** — and its condition is met: the positive
controls are passed and the detection limit is measured. The sentence reads:

> What converges in model systems does not converge in patient cells. Across
> four calibrated patient cohorts of skeletal disease, the concordance of the
> differentiation programme lies at a mean z of −0.69, and no cohort reaches its
> own detection limit. The lesion response likewise lies at chance level — as it
> does internally and externally.

**What the finding expressly does not say.** It does not say that the programme
does not exist in patients. It says: at these cohort sizes (5-8 patients) and at
this detection limit (C of roughly 0.58-0.69 required), no shared direction
across patients is measurable. Four readings stand side by side and are not
decided between:

(a) the programme really is not shared between patients;
(b) the effect lies below the cohort's own detection limit — which is supported
    by the limits being high at 5-8 patients;
(c) the patient against control contrast is something other than the
    differentiation contrast the module comes from — the most obvious reading,
    and it is not an excuse but a statement about model transfer;
(d) the cohorts are so heterogeneous in entity, tissue, platform and laboratory
    that a shared direction is lost.

**Against (d) stands the fact that each cohort was tested individually** and
none reaches its own MDE80 — it is not a pooling problem.

**GSE186141 deserves a sentence of its own.** The OI osteoblasts run clearly
below the null against the module in both sets (z −3.95 / −4.14; **not** to be
read as "opposed", see addendum 1). The control group there consists of **two**
samples; the control median is correspondingly uncertain, and the cohort has the
highest detection limit of all four (MDE80 0.689). The finding is reported but
not read as a counter-finding in its own right.

---

## 6. What was expressly not done

- No new convergence axis was sought, no gene set was built by hand, and the
  module was not readjusted.
- **No second implementation**: `_module.konkordanz` and `_module.kontrast`
  unchanged, `_marker.py` unchanged.
- No adjustment to a covariate of the baseline; the control expression appears
  solely as the stratification variable of the null.
- No follow-up search after the first result, no extension of the search axes.
- No change of the calibration assignment once the numbers were known.
- The three non-calibratable cohorts were **not** removed from the protocol,
  although two of them (GSE228522, GSE160207) show conspicuous values for the
  lesion response (z +1.77 / +1.64). They carry nothing.
- The existing negative finding `f4_krankheitsanreicherung.csv` is untouched.
- The word "specific" has not been used.

---

## 7. Output

`derived_data/M_patienten/`: `suchlauf.csv`, `treffer_roh.csv`,
`sichtung_mechanisch.csv`, `kandidaten.csv`, `kohorten_sichtung.csv`,
`eichung.csv`, `streuung.csv`, `streuung_null.csv`, `synthese.csv`,
`laesionssatz_173.csv`, `51_log.txt`. The metadata per series are in
`data_raw/_meta/`.

---

## Addendum 1 — 2026-08-22: two corrections to the report, not to the numbers

**The individual computations in §3 and §4 are valid unchanged.** What is
corrected is the **reading** in §4 and §5 and the **aggregation**. Both were
re-checked against the result files after the session of 2026-08-21.

### (1) The mean z is not an admissible aggregation

§4 reports a **simple mean** over four cohorts with a range from −3.95 to +2.38.
For exactly this heterogeneity the main part of this work uses the **study
synthesis** (Figure S3B, study-wise z against a joint donor-flip null). Here an
average was taken — a break with our own standard.

The complete leave-one-out, recomputed:

| without | the programme | the lesion response |
|---|---|---|
| — (all four) | −0.69 | −0.94 |
| **GSE186141** | **+0.39** | **+0.13** |
| GSE22855 | −0.50 | −1.55 |
| GSE292600 | −1.72 | −1.15 |
| GSE77758 | −0.95 | −1.18 |

**The mean reverses its sign as soon as the cohort with two control samples is
dropped** — the same cohort that §5 itself carries as exploratory and that has
the highest detection limit of the four.

**Consequence:** the numbers −0.69 and −0.94 are **withdrawn** as a core
sentence. They remain as a descriptive statement, always with this table beside
them. The core sentence of this phase is instead the statement that holds under
**every** omission:

> **None of the four calibrated cohorts reaches its own detection limit**
> (MDE80 C 0.576-0.689), for neither of the two gene sets.

The study synthesis is made up for in phase D (§3.6 of the donor-resolved plan
of 2026-08-22).

### (2) GSE186141 is not a counter-finding

§5 read z −3.95 as "runs against the prediction". With this metric that is **not
decidable**. What is measured is the agreement with the direction `ri` predicted
by the module. A value below the null arises both when the patients are
**disagreeing** and when they **agree** and run jointly in the **opposite**
direction. The metric does not separate the two cases.

A direction-free diagnostic (the mean absolute share of agreeing patients, with
an **unstratified** null — hence **not** a reportable test, only an indication)
points to the second case: GSE186141 would reach z +2.79 there instead of −3.95.
That becomes defensible only with a baseline-stratified null; it is preregistered
as statistic S3b in phase D.

**Consequence:** everywhere "lies below its null" instead of "runs against it".
Corrected in the study narrative, in the captions, in `PROTOCOL_M_overall.md`
and here.

### (3) What of §5 is untouched

- The calibration table in §3 and every individual number in §4.
- "No cohort above its own MDE80" — valid under every omission.
- The finding that three of seven cohorts cannot be calibrated.
- The complete screen in §2.

### (4) A third objection, only noted here

The calibration of GSE292600 and GSE77758 ran through `NAIV` and therefore
demonstrates **tissue identity**, not the ability of those cells to run the
osteogenic or chondrogenic programme. For a **module test** that is the wrong
calibration. The reassessment takes place in phase D (§0.5 b and §5 of the
donor-resolved plan of 2026-08-22) and not here, because it concerns gate B and
not only the report.

---

## Addendum 2 — 2026-08-22: the calibration, the gate and the classification

This addendum draws the consequence of the objection that addendum 1 (4) only
**noted**. It changes **no computed number**; it changes which cohort may carry
a module finding, and thereby the status of gate B. The basis: §0.5 (b) and §5
of the donor-resolved plan of 2026-08-22. At the time of this addendum **all**
the numbers of §3 and §4 as well as the leave-one-out table of addendum 1 were
known — which is expressly noted, because the reassessment is therefore made
**after the fact** and may not be presented as preregistered.

### (1) `NAIV` is not a calibration of a module test

The project rule reads: *a level without a calibration carries no finding* —
and, sharpened: **the calibration must fit the level being tested.**

The 173-gene module is defined on `dWT`, the contrast **undifferentiated against
differentiated**. A calibration with `NAIV` shows that the measured cells are
undifferentiated cells of their tissue. It does **not** show that these cells
run the osteogenic or chondrogenic programme at all — that is, not that the axis
on which the module is defined exists in that cohort. For a module test that is
the wrong positive control.

| GSE | set A | calibration | new assessment |
|---|---|---|---|
| GSE186141 | `OSTEOGEN` | passed (z +2.38) | **module-carrying** |
| GSE22855 | `CHONDROGEN` | passed (z +4.77) | **module-carrying** |
| GSE292600 | `NAIV` | passed (z +4.62) | **not module-carrying** — tissue identity, not the differentiation axis |
| GSE77758 | `NAIV` | passed (z +2.73) | **not module-carrying**; in addition EDS-HT / JHS is not a skeletal dysplasia |

The numbers of GSE292600 and GSE77758 are **not deleted**. They stay in
`streuung.csv` and in §4 and will in future be marked as *not module-carrying*.
That affects in particular the only directed individual finding of the phase
(GSE292600, z +2.38, p 0.022): it loses its calibration and is therefore **no
longer a reportable module finding** but an observation without a positive
control.

`_marker.py` stays unchanged — `NAIV` remains valid as a marker set, only no
longer as the calibration of a module test.

### (2) Gate B is not met

What remains is **two** module-carrying cohorts, GSE186141 and GSE22855. The
threshold "at least 2" would formally be touched, but:

- **GSE186141** has **two** control samples, the highest detection limit of all
  four (MDE80 0.689), and is already carried as exploratory in §5 and
  addendum 1. It is exactly the cohort that carried the mean on its own.
- **GSE22855** lies on the null in both sets (z −1.27 / +0.90).

A cohort that is itself carried as exploratory cannot carry a two-cohort
threshold. **Gate B therefore counts as not met.** The patient against control
finding of phase M-B is **exploratory** in its entirety — including the part
carried as confirmatory until now.

**What stays confirmatory is only what needs no module calibration:**

- the complete screen (§2) and its design finding,
- that three of seven findable patient cohorts measure no skeletal tissue (§3),
- that **none** of the four cohorts reaches its own MDE80 (§4) — a statement
  about the **detection limit**, not about the module, valid under every
  omission and independent of the calibration question.

### (3) The aggregation, with the procedure of the main part

Addendum 1 (1) withdrew the mean and added the leave-one-out table. That table
remains the binding version; it is **not recomputed and not reinvented** here.
The study synthesis following the procedure of Figure S3B (study-wise z against
a joint donor-flip null) is introduced in **phase D** and applied there to the
donor-resolved cells; applying it retrospectively to four cohorts with two
different calibration states would suggest a finding status that point (2) has
just removed.

An arithmetic mean over cohorts does not appear in this work again — neither
here nor elsewhere.

### (4) What phase M-B now is

Phase M-B does **not** answer the question of the scissors. It answers the
question:

> Are patients with a skeletal diagnosis shifted, relative to controls, along
> the axis that the differentiation module describes?

That is a **legitimate secondary question with a measured limit** — and the
answer is: with 5-8 patients per cohort and an MDE80 of C 0.576-0.689, not
measurable. It is not the main question, because the contrast tested (patient
minus control median) is not the contrast on which the module is defined (§0.5 a
of the donor-resolved plan).

**The consequence for the paper:** phase M-B moves out of the main figures into
the supplement (Figure S4 of the new order) and is marked there as
**exploratory**. Its place among the main figures goes to phase D.

### (5) What this addendum expressly does not do

- It deletes nothing: `PREREG_M_B.md` stays unchanged, and §1 to §7 and
  addendum 1 stand word for word.
- It recomputes no number and adjusts no threshold.
- It does **not** remove GSE292600 and GSE77758 from the result files.
- It does not change `_module.py`, `_marker.py` or `51_patient_variability.py`.

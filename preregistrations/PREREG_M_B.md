> Translated from the German original of 2026-08-21. The content, the dates
> and every number are unchanged.

# Preregistration M-B — measuring individuality patient by patient

Written and dated **2026-08-21**, **before the first download of a patient
cohort** and before the first statistic of this phase.

`PREREG_whole_study.md` and `PREREG_S6.md` remain in force unchanged; this
document governs phase B only. Changes are added as a dated addendum at the
end, stating which numbers were known at the time of the change.

---

## 0. What was known before the date of this document

**0.1 The module is fixed.** The 173 convergent genes have stood unchanged in
`derived_data/reference_tables/S5_konvergente_gene.csv` (129 up, 44 down) since
S12. They are **not recomputed and not readjusted** here. The convergence axis
as such is closed by `PREREG_S6.md` §1; this phase looks for **no** new axis,
in any construction.

**0.2 The internal counter-test is known.** `iv` does not converge internally
(7 genes against 8.0 expected, p 0.45) and does not converge externally. That
is the starting point of the prediction tested here, not its result.

**0.3 The detection limit of the statistic is known** (0.35 z at 60 genes per
side). For the patient level it is **measured afresh**, because the cohort size
and the noise structure are different (§6).

**0.4 What is not known.** No patient cohort of this work package has been
downloaded, prepared or computed. No accession of the search axes in §2 has
been verified so far. No numerical value of a patient cohort exists.

---

## 1. The prediction being tested

From the existing finding it follows:

> In patient-derived transcriptome cohorts of skeletal disease the **173
> programme genes run in the same direction between patients**, while the
> **lesion-response genes run individually**.

The directional prediction, one-sided upwards in the statistic of §4:

- **H1 (programme):** `U_programme` lies above the background-drawn null.
- **H2 (lesion response):** `U_lesion` does **not** lie above the null.
- **H3 (the scissors):** `U_programme − U_lesion` lies above the null of the
  two-set contrast.

**Both outcomes are publishable.** A null result on H1 is not a failure but a
statement about the transferability of differentiation models — **provided the
positive control in §5 passes.**

---

## 2. Cohort selection — rules fixed before the search

**Search axes** (GEO, ArrayExpress, SRA; every search logged with its date,
search string and hit count):

osteogenesis imperfecta (COL1A1/COL1A2) · mucopolysaccharidoses ·
pseudoachondroplasia and MED (COMP, MATN3) · achondroplasia and the FGFR3
spectrum · SHOX deficiency, Leri-Weill · cleidocranial dysplasia (RUNX2) and
campomelic dysplasia (SOX9) · fibrodysplasia ossificans progressiva (ACVR1) ·
free text "skeletal dysplasia", "short stature", "growth plate",
"chondrodysplasia".

**Inclusion, all conditions simultaneously:**

- E1 human, patient-derived (primary cells, patient iPSC derivatives, tissue);
- E2 **n ≥ 5 patients** with an unambiguous sample-to-patient assignment;
- E3 a control group (healthy donors or isogenic correction), n ≥ 2;
- E4 a gene-level matrix obtainable (counts, RPKM/TPM or array intensities);
- E5 at least 60 of the 173 module genes measurable on that platform.

**Exclusion codes** (the same coding as
`derived_data/reference_tables/S1_sichtung_alle_datensaetze.csv`, extended):

| code | reason |
|---|---|
| `A2` | no control group, no contrast possible |
| `A3` | no gene-level matrix |
| `A4` | cell empty (the condition is not present in the data set) |
| `A7` | not human |
| `A9` | single-cell or multiome format without a per-patient pseudobulk |
| `M1` | **fewer than 5 patients** |
| `M2` | **immortalised line** |
| `M3` | **no patient assignment per sample** — n arises from cells, wells or replicates |
| `M4` | not a skeletal disease, wrong entity |
| `M5` | fewer than 60 module genes measurable (E5) |
| `M6` | positive control not passed — "the level does not carry" (§5) |

`M2` and `M3` are **hard**: cells are not biological units. The patient is the
unit, without exception.

**An undifferentiated arm is not required here.** That is the deliberate
difference in feasibility from phase C; the statistic in §4 needs no
differentiation contrast but a patient-against-control contrast.

**Every accession is verified when it is found**, not taken from memory. The
complete screen — including the exclusions — goes into
`derived_data/M_patienten/kohorten_sichtung.csv`.

---

## 3. Preparation — fixed before the computation

- Symbols through the platform annotation of the data set; multiple probes per
  symbol are combined by the **median**.
- Count matrices: CPM, `log2(x+1)`; arrays: the normalised value deposited by
  the authors, without renormalisation.
- Genes with zero variance across all samples drop out.
- **No batch-correction step** and **no matching to a covariate of the
  baseline** — the guard "matching also matches away the target quantity".
  Baseline expression appears **only** as a stratification variable of the null
  (§4).
- Several samples of the same patient are averaged to the patient level
  **before** the contrast.

---

## 4. Statistic and null model

The **gene key** is the HGNC symbol.

**Two fixed gene sets, both from the same existing computation:**

- **Programme** = the 173 genes from `S5_konvergente_gene.csv`, direction `ri`.
- **Lesion response** = the **173** genes with the highest lesion-response
  consistency from `derived_data/R_intern/R_interne_genkarte.csv`, direction
  `iv_vz`. The ordering is deterministic, by `(|iv_cons| descending, |iv_med|
  descending, ensembl ascending)`, then the first 173. The equal size is
  deliberate: only sets of equal size are comparable under the statistic below.
  The set is **not a discovery** — it is a deterministic re-sorting of a column
  that has already been reported.
  *(The hard threshold set `|iv_cons| ≥ 0.8` comprises only 6 to 7 genes and
  lies below the minimum size of `_module.konkordanz`; it is reported
  descriptively but carries no test.)*

**One contrast per patient.** For patient *p* and gene *g*:

    delta_pg  =  x_pg  −  median over the control group

**The primary statistic `U` — agreement in direction between patients.** For
every gene of the set:

    u_g  =  (share of patients with sign(delta_pg) == s_g)

with `s_g` = `ri` or `iv_vz`. The set statistic `U` is the mean of `u_g` over
the genes of the set. `U` is exactly the concordance from
`reference_implementations/_module.py`, averaged patient-wise; **no second
implementation is written** — `51_patient_variability.py` calls
`_module.konkordanz`.

**The null model.** A baseline-stratified background draw as in `_module.py`:
gene sets of the same size from the background measurable on that platform,
drawn **within the same deciles of mean control expression**, carrying **the
same set of signs** as the set under test. **20 000 draws. Seed 20260821**
(identical to `_module.SEED`). Empirical, one-sided upwards, with a boundary
correction and then doubled — as usual in this project.

**The secondary statistic** is the continuous variant from
`_module.konkordanz` (the mean signed rank), reported unchanged.

**The scissors (H3)** are computed as `_module.kontrast` between the two sets
on the series `g -> u_g`, against the same background draw. This statistic is
insensitive to any global shift of the cohort — the same argument that saves
the per-gene z scale in the core of the paper.

**The unit is the patient.** Cells, wells, replicates and technical repeats
never count as n and appear in no formula.

---

## 5. The mandatory positive control

Before every main computation, per cohort, with **the same** statistic and
**the same** null:

- The canonical lineage markers from `reference_implementations/_marker.py`,
  matched to the tissue type of the cohort. What is tested is whether the
  cohort shows the textbook separation of its own tissue type at all.
- The cohort **passes** if the marker set matching its tissue type reaches
  `p < 0.05` in the predicted direction.

**If a cohort fails, it is discarded — with a log entry (code `M6`, "the level
does not carry"), not silently.** **No** number from a failed cohort is
reported as a finding.

**The gate:** at least **two** cohorts with a passed positive control. With
only one, the single-cohort rule continues to apply: the finding is
**exploratory** and is marked as such in the manuscript.

---

## 6. Detection limit — to be measured before the analysis

Per cohort, permutation determines which excess concordance would be found with
**80 % power** at that patient count and that set size
(`MDE80 = null mean + 2.8 x null SD`, the project threshold from `_module.py`).
**Without this number a null result of this phase is not readable** and is not
reported — the same logic as MDE80 in the main part.

---

## 7. Stopping criterion

- If the screen finds **no** cohort satisfying E1 to E5, phase B ends with
  `PROTOCOL_M_B.md` and the screening table. **The screen alone is a result**
  and is reported.
- If **no** included cohort passes the positive control, the level is dead: an
  entry in the guard list of fallen hypotheses, and no rebuilding out of
  phase B.
- There is **no** searching on until the result fits. The search axes in §2 are
  final; every extension is a dated addendum.

---

## 8. What is explicitly not done

- No new convergence axis, no new gene sets, no readjustment of the module.
- No second implementation of the statistic.
- No matching to a covariate of the baseline.
- No mixing with the confirmatory 18-data-set cohort; the patient cohorts stand
  separately.
- The word **"specific"** is not used for the scissors.
- The existing null result `f4_krankheitsanreicherung.csv` is not overwritten.

---

## 9. Output

`derived_data/M_patienten/streuung.csv`, `streuung_null.csv`,
`kohorten_sichtung.csv`, `preregistrations/PROTOCOL_M_B.md`.

---

## Addendum 1 — 2026-08-21, making the calibration concrete (§5)

**Known at the time of this addendum:** the seven cohorts of the screen with
their metadata (entity, tissue type, patient count, platform) and the fact that
the matrices have been downloaded. **No expression value has been read and no
statistic computed.**

§5 fixes the positive control but does not say how it is to be carried out on a
patient-against-control contrast. The marker set describes a **differentiation
state**, not a disease contrast; the control is therefore not a contrast
control but a **tissue-identity control** — the same question as Figure S4A for
the orthogonal levels: *does this measurement level find its own textbook
markers again?*

**The prescription.** `_module.kontrast` is computed on the **mean expression
profile of the control group**:

    set A = the marker set belonging to the tissue type, from `_marker.py`
    set B = the remaining lineage marker sets (adipogenic + myogenic +
            chondrogenic, or + osteogenic), excluding set A

It passes at `p < 0.05` in the direction "set A higher". This establishes that
the platform, the annotation and the preparation represent the tissue identity
at all.

**The assignment, fixed before the computation:**

| cohort | tissue | set A |
|---|---|---|
| GSE186141 | primary osteoblasts (bone) | `OSTEOGEN` |
| GSE22855 | cartilage / enchondroma | `CHONDROGEN` |
| GSE292600 | dermal fibroblasts | `NAIV` (mesenchymal stromal) |
| GSE77758 | dermal fibroblasts | `NAIV` (mesenchymal stromal) |
| GSE160207 | whole blood | **none** |
| GSE228522 | CD14+ monocytes | **none** |
| GSE58435 | amniotic fluid, cell-free mRNA | **none** |

**For the last three, `_marker.py` holds no canonical set of their tissue
type.** They are therefore **not calibratable** and carry no finding — code
`M6`, with the note "not calibratable, no marker set of this tissue type".
Their numbers are computed and stored descriptively but are **not reported as
a finding**. This assignment stands **before** the first statistic and is not
readjusted.

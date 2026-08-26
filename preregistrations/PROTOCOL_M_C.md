> Translated from the German original of 2026-08-21. The content, the dates
> and every number are unchanged.

# Protocol M-C — the scissors against real diagnoses

Computed **2026-08-21**, following §4 of the medical-extension plan of
2026-08-21. Script: `reference_implementations/53_diagnosis_screening.py`,
building on the search from `50_cohort_search.py` (1424 series) and the
pre-screen from `50b_screening.py` (127 candidates for screening by hand).
The numbers are in `derived_data/M_diagnosen/`.

---

## 1. The question

The same computation as in the main part, but the lesion is a **diagnosis**
instead of an engineering intervention. That requires the same 2 x 2 design as
the eighteen points:

|  | undifferentiated / undifferentiated | differentiated |
|---|---|---|
| **control** | yes | yes |
| **diagnosis** | yes | yes |

Figure 1 predicts that most candidates fail at **`A1`** — no undifferentiated arm.
**That failure is itself the result** and belongs in the main text.

---

## 2. The audit

All **127** candidates for screening by hand were audited against their **GSM
metadata** (`data_raw/_meta/<GSE>_proben.csv`), with the same exclusion codes as
`derived_data/reference_tables/S1_sichtung_alle_datensaetze.csv`. Per series the
check asked:

- **the diagnosis axis** — is there any disease assignment per sample at all?
- **the undifferentiated arm** — is there an undifferentiated starting state (day 0,
  undifferentiated, expansion medium, an iPSC or MSC starting point)?
- **the differentiated arm** — is there a differentiation arm?

Seven series survived the mechanical audit and were **checked individually by
hand**; the verdicts are named in the script (`HANDURTEIL`) and are therefore
reproducible.

---

## 3. The balance

`sichtung.csv`, all 127 series with their verdict.

| code | reason | series |
|---|---|---|
| `A2` | no diagnosis as a lesion axis | 77 |
| **`A1`** | **no undifferentiated arm** | **46** |
| `A7+M2` | mouse or immortalised line in the differentiation arm | 1 |
| `M4` | not a skeletal entity | 1 |
| `-` | a complete 2 x 2 | **2** |

**Of the 50 candidates that have a diagnosis as a lesion axis at all, 46 fail at
`A1` — 92 %.**

**The two complete 2 x 2 designs of the whole screen are both already part of
the confirmatory cohort:**

| GSE | entity | design |
|---|---|---|
| GSE218101 | MPS VI (ARSB) | day 0/14, 4 patient lines against an isogenic correction — the "ARSB / MPS VI" point of the eighteen |
| GSE221128 | FOP (ACVR1) | iMSC day 0/6, FOP against resFOP — the "ACVR1 / FOP" point of the eighteen |

**New points for the metric: none.** `punkte.csv` is empty, and that is not an
omission but the finding.

**The seven hand checks in detail:**

| GSE | verdict | why |
|---|---|---|
| GSE218101 | ALREADY | a complete 2 x 2, already one of the eighteen |
| GSE221128 | ALREADY | a complete 2 x 2, already one of the eighteen |
| GSE206213 | OUT, `A2` | Werner syndrome as a **WRN knockdown in hESC** — an engineering intervention, not a diagnosis; a line instead of a patient |
| GSE206214 | OUT, `A2` | the same study, the ChIP arm |
| GSE241507 | OUT, `A7+M2` | a human OPLL arm without a differentiation axis; the day 0/14 axis lies in **MC3T3-E1** — mouse, immortalised |
| GSE188698 | OUT, `M4` | retinal organoids, Mueller glia — not a skeletal entity |
| GSE222109 | OUT, `A2` | iPSC against sclerotome, wild type only — no diagnosis |

---

## 4. What follows from it

> **The skeletal dysplasia literature almost never sequences the
> undifferentiated state.** Of 1424 screened series and 50 candidates with a
> real diagnosis axis, **four** have an undifferentiated arm at all, and **two** a complete
> 2 x 2 design — and both are already analysed in this work.

That is the practically most important sentence for the target community, and it
is a **design finding**, not an aside: the question whether a diagnosis shifts
the differentiation response **cannot** be answered with the available public
data, because the reference state is missing. Anyone who wants to answer it must
sequence the undifferentiated arm as well — which costs nothing extra and is
nonetheless almost never done.

It also explains why phase B had to take the detour through the patient against
control contrast (`PREREG_M_B.md` §2: "an undifferentiated arm is not required here"):
phase B is exactly the design that can be run with the available data, and
phase C shows why there is no other.

**The cohort of eighteen stays untouched.** No point was added and none
removed; the screen stands beside `S1_sichtung_alle_datensaetze.csv`, not inside
it.

---

## 5. What was expressly not done

- No candidate was included in order to raise the number of points.
- No `A1` exclusion was circumvented by substituting the undifferentiated arm from another
  study — points must come from **one** study.
- The two complete 2 x 2 designs were **not** computed a second time and issued
  as new evidence; they are the same points.
- The screen was not broken off once it was clear that nothing new would be
  added.

---

## 6. Output

`derived_data/M_diagnosen/sichtung.csv` (127 series with their verdict) and
`punkte.csv` (empty, with a header row). The metadata per series are in
`data_raw/_meta/`.

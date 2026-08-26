> Translated from the German original of 2026-08-21. The content, the dates
> and every number are unchanged.

# Protocol M-A — the human-genetics anchor, at full power

Computed **2026-08-21**, following `PREREG_M_A.md` including its addendum 1
(both dated before the first figure). Scripts:
`reference_implementations/52a_go_sets.R`, `52b_fetch_panels.py`,
`52_human_genetics_anchor.py`. Seed 20260821, 20 000 draws. The numbers are in
`derived_data/M_humangenetik/`.

---

## 1. The question

Are the 173 convergent programme genes enriched for skeletal dysplasia, short
stature or body-height genes — at **markedly better power** than in the earlier
version (357 and 50 panel genes respectively, with a detection limit around
OR 1.6)?

---

## 2. The panels — source, version, retrieval date

`panels.csv`. Retrieved 2026-08-21 except where noted otherwise.

| panel | source | genes in total | **in the background** |
|---|---|---|---|
| `NOSO` | HPO `HP:0002652` plus 10 subterms | 188 | 154 |
| `NOSO_BREIT` | PanelApp 309 **all levels** (v10.3) united with the HPO subtree | 689 | 539 |
| `KLEIN` | HPO `HP:0004322` plus 24 subterms | 1483 | 1130 |
| `KLEIN_BREIT` | PanelApp 1471 all levels (v2.9) united with the HPO subtree | 1187 | 1144 |
| `GWAS` | GWAS catalogue, `MAPPED_TRAIT = "body height"` (identifier `OBA_VT0001253`), 45 209 associations, `MAPPED_GENE` | 12 435 | **5654** |
| `PA309` | PanelApp 309, level 3 — **the earlier version** | 440 | 355 |
| `PA1471` | PanelApp 1471, level 3 — the earlier version | 73 | 62 |

**Two deviations from the plan, both documented before the first number
(addendum 1):**

1. The plan names the *Nosology of Genetic Skeletal Disorders* (about 770
   genes). A machine-retrievable, versioned gene list of the Nosology does not
   exist; `NOSO_BREIT` (689 genes) is the substitute and matches the order of
   magnitude.
2. The plan names MAGMA for the GWAS aggregation. MAGMA is not available in this
   working environment; what is used is the fallback rule foreseen in the plan —
   the catalogue's own `MAPPED_GENE` assignment (the overlapping or nearest
   gene). That is the weaker assignment and is reported as such.

**The gain in power is real:** short stature goes from **50** to **1130** panel
genes, and a height GWAS axis with **5654** genes is added that did not exist at
all before.

---

## 3. Background and matching

The background is the gene pool of the computation: **11 581 genes** with both a
baseline expression (`basis_med`) and a union exon length (GENCODE v46). For
each set gene the null draws a background gene from the **same cell** of a
10 x 10 grid (expression decile x length decile). Without that matching the odds
ratio is not interpretable — disease genes are longer and more highly expressed
than the average.

**Why this is not a forbidden adjustment.** The guard "matching also matches the
target quantity" concerns covariates of the **baseline on the z scale**: there
the baseline correlates structurally with `dWT` (median cor −0.566), and
matching removes part of the signal. Here the target quantity is **panel
membership** — an external annotation from HPO, PanelApp and the GWAS catalogue,
which is not a quantity of these data sets and stands in no structural relation
to `dWT`. As a control, the **raw, unmatched four-field odds ratio** is reported
throughout; it stands in the same table and changes no verdict.

---

## 4. Gate A — the positive control, both parts passed

`eichung_A.csv`.

**(a) Lineage markers in the dysplasia panel.** `OSTEOGEN` plus `CHONDROGEN`
from `_marker.py`, 17 of the 30 symbols measurable in the background:

| panel | in the panel | null expectation | OR | z | p |
|---|---|---|---|---|---|
| `NOSO` | 4 of 17 | 0.30 | 17.1 | +6.84 | 0.0004 |
| `NOSO_BREIT` | 11 of 17 | 0.85 | 34.8 | +11.33 | 0.0001 |
| `PA309` | 11 of 17 | 0.58 | 51.6 | +13.91 | 0.0001 |

**(b) The anchor.** `S_DISTAL` against `S_BIOSYN`, the GO sets taken
**verbatim** from `03_metric.R` and exported with `52a_go_sets.R` (545 and 2711
genes; 416 and 1616 in the background):

| panel | distal | biosynthetic | OR | p |
|---|---|---|---|---|
| `NOSO` | 15/416 | 16/1616 | **3.74** | 4.3 x 10⁻⁴ |
| `NOSO_BREIT` | 44/416 | 50/1616 | **3.70** | 4.0 x 10⁻⁹ |
| `PA309` | 37/416 | 27/1616 | **5.75** | 2.9 x 10⁻¹¹ |

The anchor of the earlier analysis (OR 4.97, p 6.7 x 10⁻¹¹ on the 357-gene
panel) **reproduces** — on the same panel, in the version computed here,
OR 5.75, p 2.9 x 10⁻¹¹. The design and the background are therefore set
correctly.

**GATE A PASSED.** The numbers of this phase carry.

---

## 5. The result

`anker.csv`, with the detection limits separately in `anker_power.csv`.
Bonferroni over 7 panels x 2 sets = 14 comparisons, threshold **p < 0.0036**.

**The programme (173 genes, 147 in the background)**

| panel | observed | null expectation | OR matched | z | p | **OR at 80 % power** | OR raw |
|---|---|---|---|---|---|---|---|
| `GWAS` | 80 | 79.99 | **1.00** | +0.00 | 0.92 | **1.59** | 1.26 |
| `KLEIN` | 11 | 13.72 | 0.79 | −0.78 | 0.36 | **1.85** | 0.75 |
| `KLEIN_BREIT` | 12 | 13.81 | 0.86 | −0.52 | 0.53 | 1.85 | 0.81 |
| `NOSO` | 2 | 2.42 | 0.82 | −0.27 | 0.61 | 2.86 | 1.02 |
| `NOSO_BREIT` | 8 | 7.11 | 1.13 | +0.35 | 0.84 | 2.13 | 1.18 |
| `PA309` | 4 | 4.66 | 0.85 | −0.31 | 0.62 | 2.37 | 0.88 |
| `PA1471` | 1 | 0.92 | 1.09 | +0.09 | 0.80 | 3.98 | 1.28 |

**The lesion response (173 genes, of equal size)**

| panel | observed | null expectation | OR matched | z | p | OR at 80 % power |
|---|---|---|---|---|---|---|
| `GWAS` | 94 | 84.13 | 1.26 | +1.53 | 0.15 | 1.52 |
| `KLEIN` | 18 | 16.90 | 1.07 | +0.28 | 0.86 | 1.76 |
| `KLEIN_BREIT` | 18 | 17.09 | 1.06 | +0.23 | 0.90 | 1.76 |
| `NOSO` | 4 | 2.43 | 1.66 | +1.02 | 0.45 | 2.84 |
| `NOSO_BREIT` | 5 | 8.54 | 0.57 | −1.26 | 0.13 | 2.02 |
| `PA309` | 2 | 5.46 | 0.36 | −1.52 | 0.051 | 2.25 |
| `PA1471` | 2 | 0.85 | 2.36 | +1.24 | 0.42 | 4.09 |

**None of the fourteen comparisons reaches the Bonferroni threshold.** The
largest single value is `lesion response x PA309` at p 0.051 in the
**opposite** direction (OR 0.36) — under 14 comparisons exactly what chance
produces.

---

## 6. What follows from it

**The existing negative finding is confirmed and considerably sharpened.** The
earlier version said: no enrichment at 357 and 50 panel genes and a coarse
detection limit around OR 1.6. This version says the same — but against a short
stature panel with **1130** genes (detection limit OR 1.85), against a height
GWAS axis with **5654** genes (detection limit **OR 1.59**) and against a
dysplasia panel of the Nosology order of magnitude, each matched on expression
and length.

> On the body-height GWAS axis, where the detection limit at OR 1.59 is
> sharpest, the enrichment of the programme lies at **OR 1.00** — 80 observed
> against 79.99 expected genes. A number cannot come closer to nothing.

**No contradiction with the earlier version.** All signs and orders of magnitude
are the old ones; `PA309` and `PA1471` were carried along here unchanged and
deliver the same non-findings. `f4_krankheitsanreicherung.csv` stays unchanged
in the paper and is **extended, not replaced**, by this table.

**The lesion response does not behave differently.** It too is enriched on no
panel. The human-genetics anchor therefore does not separate programme from
lesion response — unlike every metric at the transcriptome level.

**And the anchor itself still stands.** Distal against biosynthetic secretion
separates skeletal dysplasia genes at OR 3.7-5.8 (p down to 2.9 x 10⁻¹¹) — the
same computation, the same background, the same matching in which the convergent
genes show nothing. **Both stand side by side and are not harmonised.** The
reading offered earlier — that disease genes are defined through their
*consequences of loss*, not through their *regulatory dynamics* — remains a
hypothesis and is not declared a measurement. It now has a number at its side,
though: a GO axis of secretion logic finds the disease genes, a data-driven axis
of differentiation dynamics does not.

---

## 7. What was expressly not done

- No panel was exchanged, extended or trimmed once its result was known; all
  seven are reported, including the weak ones.
- No change of background once the numbers were known.
- The module was not readjusted and no new axis was sought.
- The GO sets of the positive control were taken **verbatim** from
  `03_metric.R`, not redefined.
- The raw Fisher odds ratio is reported throughout, so that the matching can be
  checked.
- The word "specific" was not used.

**One implementation decision deserves naming:** the minimum size for the
matched computation stands at 8 measurable genes — the same bound as in
`_module.konkordanz`. It was lowered from 20 to 8 **before** any marker
enrichment was visible; at 20, positive control (a) would have failed for the
sole reason that the canonical marker sets are small (17 of the 30 symbols
expressed in the computation pool).

---

## 8. Output

`derived_data/M_humangenetik/`: `panels.csv`, `go_saetze.csv`, `eichung_A.csv`,
`anker.csv`, `anker_power.csv`, `_exonlaengen.csv`, `52_log.txt`. The raw panel
data are in `data_raw/_panels/`.

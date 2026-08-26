> Translated from the German original of 2026-08-21. The content, the dates
> and every number are unchanged.

# Preregistration M-A — the human-genetics anchor, at full power

Written and dated **2026-08-21**, **before the reference panels were
downloaded** and before the first statistic of this phase. Phase M-B is
complete at this point (`PROTOCOL_M_B.md`); its result enters **no** decision
of this document.

Changes are added only as a dated addendum at the end.

---

## 0. What was known before the date of this document

**0.1 The existing null result.**
`derived_data/manuscript/f4_krankheitsanreicherung.csv`: the dysplasia panel of
357 genes (Genomics England PanelApp 309, confidence level 3), the short
stature panel of 50 genes (PanelApp 1471, level 3), with the overlap removed.
All odds ratios 0.71 to 1.65, all p ≥ 0.27; the coarse detection limit is
around OR 1.6. **This result stays in the paper and is not overwritten.** This
phase may extend it; a contradiction is reported as a contradiction.

**0.2 The anchor.** "Distal against biosynthetic secretion" separates skeletal
dysplasia genes: 39/523 against 35/2 192, **OR 4.97, p 6.7e-11**. The two gene
sets are the GO sets `S_DISTAL` and `S_BIOSYN` from
`reference_implementations/manuscript/methods/03_metric.R` and are taken over
**unchanged**.

**0.3 The module is fixed** (173 genes) and is not readjusted. The convergence
axis is closed by `PREREG_S6.md` §1.

**0.4 What is not known.** None of the new reference panels has been
downloaded. No enrichment number against a new panel has been computed.

---

## 1. The question

Are the 173 convergent programme genes enriched for skeletal dysplasia or short
stature genes — **at markedly better power** than in the earlier version? And
does the lesion response behave differently?

The directional prediction: **no**, in continuation of 0.1. A null result is
readable only if (a) the positive control in §5 passes and (b) the detection
limit in §6 is measured alongside.

---

## 2. Gene sets — the test and its counterpart

- **Programme** = the 173 genes from `S5_konvergente_gene.csv`.
- **Lesion response** = the 173 genes of highest `iv` consistency,
  deterministically by `(|iv_cons| down, |iv_med| down, ensembl up)`, as fixed
  in `PREREG_M_B.md` §4 and stored in
  `derived_data/M_patienten/laesionssatz_173.csv`.

Both sets are of equal size and come from the same computation.

---

## 3. Reference panels — each with its retrieval date and version

| key | source | retrieved |
|---|---|---|
| `NOSO` | HPO, term **HP:0002652 "Skeletal dysplasia"** including its subterms, genes from `phenotype_to_genes.txt` | 2026-08-21 |
| `KLEIN` | HPO, term **HP:0004322 "Short stature"** including its subterms, from the same file | 2026-08-21 |
| `GWAS` | GWAS catalogue (EBI), trait **body height**, aggregated per gene over the genes the catalogue assigns | 2026-08-21 |
| `PA309` | Genomics England PanelApp 309, confidence level 3 — **the earlier version**, unchanged | 2026-08-20 (existing) |
| `PA1471` | PanelApp 1471, level 3 — the earlier version | 2026-08-20 (existing) |

`PA309` and `PA1471` run alongside so that the new computation is calibrated
against the old one and any contradiction becomes visible.

**For `GWAS` the window rule applies:** the `MAPPED_GENE` entries deposited by
the GWAS catalogue itself are used (the nearest or overlapping gene, by the
catalogue's rule). MAGMA is not available in this working environment; the
deviation from the plan ("MAGMA if available, otherwise the nearest gene with
the window stated and the rule documented") is therefore the fallback the plan
provides for, and it is reported as such.

All panels are brought onto **Ensembl gene identifiers** and intersected with
the computation pool.

---

## 4. Background and statistic

The **background** is the gene pool of the computation
(`derived_data/R_intern/R_interne_genkarte.csv`), **matched on expression and
on length**:

- the expression stratum is the decile of `basis_med` (mean baseline
  expression);
- the length stratum is the decile of the union exon length from GENCODE v46;
- genes lacking either quantity drop out, and they do so equally in **all**
  sets.

**Two statistics per (set x panel):**

1. **The raw two-by-two odds ratio** with a Fisher p — the same computation as
   in `f4_krankheitsanreicherung.csv`, so that the numbers can be read side by
   side. It is **not** corrected for matching and is reported only for
   continuity.
2. **The matched null (the decisive one).** 20 000 draws: for every set gene a
   background gene is drawn from **the same expression-by-length cell**, and
   the overlap with the panel is counted. Reported are the observation, the
   null mean, the null standard deviation, z, the empirical two-sided p and
   `MDE80 = null mean + 2.8 x null SD`, expressed as an odds ratio.
   **Seed 20260821**, 20 000 draws — as everywhere in this project.

The stratification is **not a matching of the target quantity**: panel
membership is an external annotation, not a quantity of these data sets. The
guard "matching also matches away the target quantity" does not apply here, and
the protocol states the reason explicitly.

---

## 5. The mandatory positive control

Both parts must pass. The computation is **the same** as in §4.

**(a) Lineage markers within the dysplasia panel.** The canonical sets
`OSTEOGEN` and `CHONDROGEN` from `_marker.py` must be enriched in the panel
`NOSO` against the matched null (`p < 0.05`, upward). If they are not, either
the panel or the background is set up wrongly.

**(b) The anchor.** The contrast `S_DISTAL` against `S_BIOSYN` must reproduce
on `NOSO`: expected **OR about 5** in the distal direction. It passes at
OR > 2 and p < 0.001. The two GO sets are exported unchanged from
`org.Hs.eg.db` by `reference_implementations/52a_go_sets.R`.

**If either part fails, no number of this phase is reported** — the phase then
counts as "the level does not carry".

---

## 6. Detection limit

A permutation computation per (set x panel): which true odds ratio would be
found with **80 % power** (alpha = 0.05, two-sided) at this panel size and this
set size? Determined by simulation over a grid of odds ratios against the same
matched null. **Without this number the null result of this phase is not
reported.**

---

## 7. Stopping rules and rules of honesty

- No panel is exchanged, extended or trimmed after its result is known.
- No change of background after the numbers are known.
- The earlier version `PA309`/`PA1471` is reported alongside, and especially so
  if it contradicts the new version.
- The word "specific" is not used.
- A result of phase M-B changes nothing in this document.

---

## 8. Output

`derived_data/M_humangenetik/anker.csv`, `anker_power.csv`, `panels.csv`,
`preregistrations/PROTOCOL_M_A.md`.

---

## Addendum 1 — 2026-08-21, the nosology substitute and the GWAS trait key

**Known at the time of this addendum:** only the **sizes** of the panels fixed
in §3 and the column structure of the source files. **No enrichment number, no
contrast and no positive control has been computed.**

**(a) The HPO term is narrower than the nosology.** `HP:0002652` with its 10
subterms yields 188 symbols (154 in the computation pool) — fewer than the
existing 357-gene panel, and therefore **no** gain in power on the dysplasia
axis. The nosology holds about 770 genes. Two broad panels are therefore added,
**in addition to those fixed in §3, which run on unchanged**:

    NOSO_BREIT  = PanelApp 309 (ALL confidence levels, 641 genes)
                  united with the HP:0002652 subtree
    KLEIN_BREIT = PanelApp 1471 (all levels, 172 genes)
                  united with the HP:0004322 subtree

**All seven panels are reported** — `NOSO`, `KLEIN`, `GWAS`, `PA309`,
`PA1471`, `NOSO_BREIT`, `KLEIN_BREIT`. None is selected after its result is
known; the Bonferroni correction in §4 runs over all seven.

**(b) The GWAS trait key.** In the version retrieved, the EBI catalogue lists
"body height" under `OBA_VT0001253` rather than under `EFO_0004339`. The filter
is therefore on `MAPPED_TRAIT == "body height"`, and the key found is named in
the protocol. The window rule from §3 (the catalogue's `MAPPED_GENE`) is
unchanged.

**(c) Bonferroni.** The main statement is made over **7 panels x 2 gene sets =
14** comparisons; the threshold is accordingly 0.05/14 = 0.0036. It stands
here, before the first number.

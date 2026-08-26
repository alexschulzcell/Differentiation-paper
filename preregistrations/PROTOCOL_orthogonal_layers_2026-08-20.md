> Translated from the German original of 2026-08-20. The content, the dates
> and every number are unchanged.

# Protocol — orthogonal measurement levels (chromatin, methylome)

Date: **2026-08-20**
Status: **exploratory, not preregistered.** This protocol was written AFTER the
computation and documents what was computed — it is not a preregistration and
is expressly not presented as one in the manuscript. Figure 5 carries the same
status as Figure 3: a hypothesis for future data, not a confirmation.

Why it has to be this way: the draft S13 preregistration aimed at a single
locus (SERPINA3) and was not usable for the module question. Rather than making
it fit after the fact — which would have made the preregistration worthless —
the whole level is openly carried as exploratory.

---

## 1. The question

Can the direction of the fixed 173-gene module
(`derived_data/reference_tables/S5_konvergente_gene.csv`, column `ri`) be
recovered on measurement levels that measure no RNA — and does the same
asymmetry hold there as internally, namely differentiation (`dWT`) yes, lesion
response (`iv`) no?

## 2. The data sets

| level | accession | cell type | axis | biological units | use |
|---|---|---|---|---|---|
| ATAC-seq, peak count table | `GSE224251` | hMSC | osteogenic against growth medium | **3 biological replicates** per arm, two surfaces | **the main data set** |
| H3K27ac ChIP-seq, BigWig | `GSE129031` | BM-MSC | chondrogenic day 14 against day 0 | 2 donor lines (`8A`, `2454e`) | **the second carrying data set** |
| ATAC-seq, BigWig | `GSE332758` | an MSC line | osteogenic / adipogenic against 0 d | 1 line, 3 time points | descriptive |
| Illumina 450K | `GSE129266` | BM-MSC | chondrogenic day 14 against day 0 | 2 paired donors | supplementary |
| Illumina 27K | `GSE33896` | hASC | osteogenic / myogenic against undifferentiated | 3 paired donors | supplementary |

**On the follow-up search of 2026-08-20 (a second replicated chromatin
cohort).** A systematic follow-up search over GEO (three independent search
strategies, 395 unique series, 88 with a matching assay and differentiation
reference) found **no** second ATAC data set with an undifferentiated arm, primary cells and
n of at least 3. Checked and rejected: `GSE270602` (2 donors), `GSE151315` and
`GSE151311` (2 replicates), `GSE113253` (hTERT-immortalised lines), `GSE239277`
(no undifferentiated arm), `GSE331503` (2 samples, no differentiation arm), `GSE310733`
and `GSE310734` (no growth-medium arm). `GSE224265` is the SuperSeries of the
already used `GSE224251` and therefore not an independent data set. The ENCODE
collection could not be queried from the working environment (DNS blocked) and
remains an open gap.

Instead of a second ATAC cohort, `GSE129031` was taken up: H3K27ac of the same
in vitro chondrogenesis, from the same laboratory and with the same protocol as
`GSE129266`, the methylation data set of that axis. The comparison of levels
(chromatin against promoter methylome) thereby becomes possible **within one
differentiation axis** instead of across two different axes. Two donor lines are
below the n-of-3 mark; the level is descriptive at cohort level, and the
inferential level is the gene. **The donor labels of the two series do not
agree** (`8A` and `2454e` against `Donor 1..4`) — it is the same system, but not
demonstrably the same donors, and it is stated that way.

**Excluded, and why**
- `PXD073297` (proteome, ADSC osteogenic): the deposited `proteinGroups.txt`
  contains only **153 protein groups** and 43 LFQ columns without a resolvable
  assignment to arms and donors. Not a viable proteome level.
- `PXD003978` (secretome, 3 tissues x 5 donors): **no differentiation axis**.
  Not usable for the module question; not reported.
- Cell lines (MG-63, TE 32.T, RD in `GSE33896`) and unpaired arms (`GSE33896`
  donor 5; `GSE129266` donors 3 and 4) are excluded. Cells and wells are not
  biological units.

## 3. Preparation

**Chromatin `GSE224251`.** Peaks with a gene assignment from the ChIPseeker
annotation supplied with the series. The gene signal is the sum of the peak
counts. CPM against the full library, `log2(x + 1)`. The axis is
`mean(osteogenic) − mean(undifferentiated)`, formed separately per surface and then
averaged. Eight preparations in a cross: (all peaks / promoter peaks only) x
(unfiltered / peak sum at least 20) x (CPM / CPM plus quantile normalisation);
**all eight are reported**, and none is selected.

**H3K27ac `GSE129031`.** The same geometry as the accessibility level
(promoter, TSS ± 10 kb, TSS ± 50 kb), so that the levels are compared with the
same window. Median normalisation, the `log2` ratio `CHON/MSC` per donor line,
then averaged over the two lines. The expected sign is `+ri` — H3K27ac marks
active promoters and enhancers and runs in the same direction as transcription.
In addition the axis is repeated on the ratio of mark to input per sample — the
input measures mappability, copy number and the tendency to fragment, and an
H3K27ac difference that also appears in the input is not a difference of the
mark.

**Chromatin `GSE332758`.** Four windows (promoter −2000/+500, TSS ± 10 kb,
TSS ± 50 kb, gene body), the signal through
`rtracklayer::summary(..., defaultValue = 0)`, median normalisation, the `log2`
ratio against `MSC-0d`.

**Methylome.** β from the non-normalised signals, `β = M / (M + U + 100)`. The
probe-to-gene assignment only near the promoter (27K: |d(TSS)| at most 1500 bp;
450K: `UCSC_RefGene_Group` in TSS1500, TSS200, 5'UTR or 1stExon). The difference
is formed donor-paired, then averaged over donors, then averaged per gene over
its probes.

## 4. Sign rules (fixed before the computation)

- **Chromatin:** open chromatin and transcription run in the same direction →
  the expected sign is `+ri`.
- **Methylation:** promoter methylation and transcription run in opposite
  directions → the expected sign is `−ri`.

Neither rule was turned around afterwards.

## 5. The metrics

All levels run through **one** implementation
(`reference_implementations/_module.py`, functions `konkordanz` and `kontrast`).

1. **Directed agreement.** The mean signed rank of the level difference over the
   module genes, against a null that draws gene sets of the same size from the
   measurable background and assigns them **the same set of signs**. Global
   drift of the level and the sign imbalance of the module (129 up, 44 down) are
   thereby neutralised.
2. **Baseline stratification.** The same null, but drawn decile-wise by the
   internal starting level `basis_med`. That is the control against the
   confounder Figure 2 warns of: highly convergent genes start systematically
   low, and a level that sees only the starting level would produce the same
   finding without any convergence at all.
3. **Calibration.** The same test, applied to the canonical lineage marker set
   of the respective axis (`reference_implementations/_marker.py`, textbook
   knowledge, fixed before the computation, disjoint). A level that does not
   find its own markers can say nothing about the module.
4. **Dose-response.** Eighths of the internal convergence **strength** `|kons|`
   against the mean directed rank, with 2 000 gene bootstraps per eighth, plus a
   trend test (Spearman over the eight eighths).
5. **A size control.** 2 000 random subsets of the module at marker size; what
   is reported is the share that beats the marker set.

The seed is **20260821** throughout, with 20 000 draws per null.

## 6. What was expressly NOT done

- No new gene list, no post-selection within the 173, no relaxation of a
  threshold.
- No second implementation of the metric.
- No selection of a preparation after the result; all eight stand in
  `derived_data/manuscript/f7_robustheit.csv`.
- No combination of the levels into a joint p value. The levels are not
  independent enough of one another (`GSE332758` delivers two axes from one
  sample, `GSE33896` two from the same donors), and there is no joint null model
  for that.

## 7. The state of the results (the numbers are in `derived_data/reference_tables/S9_orthogonale_ebenen.csv`)

- **The module finding stands on the replicated chromatin level:** z = +3.08
  (p = 0.0031), across all eight preparations z = +3.13 to +3.72, and it
  survives the baseline stratification.
- **The dose-response carries it:** the trend over the convergence eighths is
  ρ = +0.88 (p = 0.0039) for `dWT` against ρ = +0.55 (p = 0.16) for `iv`.
- **The canonical osteogenic marker set finds nothing on the same level**
  (z = −1.11); 98.2 % of the size-matched module subsets beat it.
- **On the methylome the module is not resolvable**, not even where the level
  finds its own markers cleanly (450K chondrogenic: markers z = +2.94, module
  z = +1.84, p = 0.066).
- **`GSE332758` (a single line) behaves inconsistently** — there the `iv` dose
  lies above the `dWT` dose. That is the reason the level is carried as
  descriptive and not as evidence.

## 8. The stopping and honesty rules that were observed

- The null finding on the methylome is reported as a null finding, not computed
  away by changing the window. The variants examined stand in
  `derived_data/A_dnam/A_dnam_sensitivitaeten.csv`.
- The first version of the chromatin analysis was faulty (it averaged only over
  covered base pairs instead of over the whole window). It is archived and is
  not reported.
- The first version of the calibration tested the separation of lineages, while
  the test concerned the differentiation axis. It too is archived; the
  calibration reported tests the axis at issue.

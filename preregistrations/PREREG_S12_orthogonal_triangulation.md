> Translated from the German original of 2026-08-20. The content, the dates
> and every number are unchanged.

# Preregistration S12 — orthogonal data-set triangulation

Date: **2026-08-20**
Status: written before any S12 number and before any S12 download.

## 1. Purpose and scope

S12 is a new, self-contained computational evidence session. It is meant to
raise the stability of the positive differentiation finding through genuinely
orthogonal data sets: a different study, different biological units and, where
possible, a different measurement modality. S11 remains unchanged. S12 may not
relabel the S11 result retrospectively and may not change any S11 threshold,
null or gene list.

The negative result of the S11 routes A and B is not quietly overwritten. Given
sufficient new evidence, S12 may deliver a new, clearly delimited triangulation
finding; a change to the paper's narrative or a submission decision is not part
of this session.

## 2. A fixed order and layers of evidence

The routes are tested in this order. The decision which route is computed may
not depend on the result of a previous route; a qualified candidate is
analysed even if another route was negative.

### Route A — donor-resolvable transcriptome replication (priority)

The question: does the fixed 173-gene set reproduce its known direction in
independent human differentiation studies?

Admissible are human bulk RNA-seq, microarray or scRNA-seq data sets with a
clear undifferentiated and differentiated axis. For scRNA-seq the biological
unit is the donor or sample; cells are aggregated to pseudobulk within donor
and time point. At least three independent biological units per time point are
required for a strong replication finding. Two different studies count as two
carriers of evidence only if they do not use the same donors or the same source
matrix.

Data sets already analysed in S11 do not count as new confirmation:
`GSE200492`, `GSE166824`, `GSE324998`, `GSE255646`. `GSE200492` may be
referenced as an S12 anchor but not counted again as an independent success.

### Route B — an orthogonal regulatory level

The question: do independent ATAC-seq, H3K27ac ChIP-seq, H3K4me1 ChIP-seq or
comparable epigenomic differentiation data show a matching regulatory direction
for the fixed 173-gene set?

Admissible are human mesenchymal osteogenic or chondrogenic differentiation
data with at least two biological units per state. A study may count as an
orthogonal second level only if the measurement was not merely derived from the
same RNA matrix.

Peaks are assigned in advance, using hg38 and GENCODE, to the promoters from
−2 kb to +500 bp relative to the transcription start site. Where the data set
provides gene activity values from the authors, those are used unchanged. Per
gene the direction is differentiated minus undifferentiated; the same `ri`
direction from S5 is tested. Enhancer or peak assignments outside the
predefined rule are not added after the result is known.

### Route C — an independent phenotype link

The question: does the fixed module score predict a measured differentiation
phenotype in independent samples — mineralisation, ALP or matrix formation, for
instance?

Admissible are human data sets with a joint expression measurement and a
quantitative phenotype for at least five independent donors or samples. The
score is defined in advance as the mean of the 173 per-gene expression values,
standardised within the data set and each multiplied by `ri`. The Spearman
correlation between score and phenotype is computed; the unit remains the donor
or sample. A data set without independent samples or without a quantitative
phenotype is descriptive only.

### Route D — independent perturbation data

The question: does an independent CRISPR, RNAi or chemical perturbation study
change the fixed differentiation module in the expected direction?

Admissible are human mesenchymal differentiation studies with at least three
independent biological units per perturbation arm or time point. The
perturbation may not be selected from the S5 gene list after the result is
known. The direction of the perturbation is defined in advance from the logic
of the study metadata; no new gene list and no new module may be formed.

## 3. The fixed gene list and the roles

The only confirmatory gene list is unchanged:
`derived_data/reference_tables/S5_konvergente_gene.csv`. All 173 genes and `ri`
are read in. There is no new gene discovery, no new cut-off, no rank-based
post-selection and no selection of individual genes after seeing the S12
results.

The roles are fixed:

- confirmatory: the 173 genes and their S5 direction `ri`
- not confirmatory: new markers, modules formed afterwards, and freely chosen
  candidate genes

## 4. Route A: statistic, null and success

For every qualified study the same statistic as in S11 is computed as primary:
the share of the measurable 173 genes whose external direction reproduces the
known `ri`. The comparison is differentiated minus undifferentiated; with
several differentiation time points the earliest clearly differentiated one is
used. Missing genes stay in the denominator. The biological unit is the donor
or sample, not the cell.

The noise expectation uses 10 000 permutations with a joint sign flip of the
complete biological unit or sample columns. The sample structure is preserved.
Per study the share, the Wilson 95 % interval, the null mean, the null SD, the
z value, the two-sided permutation p value and MDE80 are reported. MDE80 is, as
in the project standard, `null mean + 2.80 x null SD`.
A leave-one-biological-unit-out computation is run for every study.

A study is a **strong route A carrier** if it has at least three independent
units per time point, the share lies above MDE80, the permutation p value is at
most 0.05 and no single leave-one-out step brings the finding to or below the
null mean.

The new triangulation finding "external transcriptome replication" may be
stated only if at least **two independent studies** are strong route A
carriers. A single study is reported as an indication, not as confirmation.
Routes B, C and D can deliver independent layers of evidence but are not
quietly redeclared as route A replicates.

## 5. Route B statistic

For every gene with measurable promoter or gene activity the direction of
differentiation is determined and tested against `ri`. The denominator is only
those genes with predefined measurability; there is no additional peak or
fold-change threshold. The null, the Wilson interval, MDE80 and the
leave-one-unit-out rule follow route A.

An orthogonal regulatory layer counts as strong if at least one independent
data set shows the fixed direction above MDE80 with p ≤ 0.05 and has at least
two biological units per state. A bare enrichment number without a biological
replicate structure is descriptive only.

## 6. Route C statistic

The fixed signed module score is computed per sample. The primary number is the
Spearman correlation with the pre-existing quantitative differentiation
phenotype. Significance is determined by 10 000 permutations of the complete
sample or donor assignment with the sample structure preserved; genes are not
permuted individually. MDE80 is reported for the magnitude of the correlation
as the predefined null-mean-plus-2.80-null-SD threshold. With fewer than five
biological units the layer is not confirmatory.

## 7. Candidate search

The search runs in GEO/NCBI, ArrayExpress, ENCODE and suitable proteomics and
phenotype repositories with fixed terms:

- `human MSC osteogenic differentiation donor RNA-seq`
- `human mesenchymal chondrogenic differentiation replicate`
- `human MSC osteogenic ATAC-seq H3K27ac donor`
- `human MSC mineralization transcriptome donor`
- `human MSC osteogenic CRISPR RNA-seq`

The metadata classification happens before any expression or peak download. The
GEO candidates to be tested with priority are `GSE202080` and further hits with
a clear donor or replicate structure. `GSE288316` is not counted as an
independent donor study because of its single immortalised hBMSC-TERT4 line.
`GSE241130` is not counted as confirmatory without a documented biological
replicate structure. Further candidates are taken up by the same rule, not by
their result.

## 8. Shared quality-control and replicate rules

- The authors' quality control is adopted and logged; no result-driven quality
  control.
- Cells are never treated as biological replicates.
- Several wells of the same preparation are not donors.
- For scRNA-seq, pseudobulk is aggregated within donor and time point.
- Ensembl versions are unified technically at base level only; multiple hits
  are not double counted.
- Unclear donor identity, missing time or arm metadata, or a missing biological
  unit lead to exclusion or to purely descriptive status.

## 9. Stopping criteria

- **S12-AB1:** no two independent route A studies reach strong carrier status.
  Route A is then reported as not confirmatory; B, C and D continue to be
  tested by the rules fixed in advance.
- **S12-AB2:** an analysis needs the gene list, the replicate rule, a
  threshold, the mapping window, a population or a reference to be relaxed
  after a result is known. Stop the session and require a decision by the
  authors.
- **S12-AB3:** a robust orthogonal finding contradicts a statement in the
  paper. Report the finding; no autonomous manuscript or submission decision.
- **S12-AB4:** no layer has a sufficient biological replicate structure. Report
  the limitation; do not upgrade any level of evidence.

## 10. Technical rules and outputs

This preregistration must exist before any download. Archive size and SHA-256
go into the download log. RDS files are loaded and released one at a time;
existing partial files are skipped and not overwritten. There is no second
implementation of `kern()`, `figure_style/publication_style.R` stays unchanged,
and the word "specific" is not used for the scissors.

Results are written with a seed, package versions and sources, together with
reproducible CSV and log files per data set, a report, and candidate, exclusion
and metadata lists. No main figure and no automatic change to the manuscript.

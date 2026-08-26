# Supplementary figure and table legends

Conventions are those of the main figure legends: every legend states n, the
null model, the detection limit (MDE80) and the status of the panel.

---

# Supplementary figures

## Figure S1 · Scale critique and self-test of the metric, Related to Figure 1

Methodological; confirmatory as a calibration of the metric itself. It documents
how the metric behaves, and why several intuitive-looking analyses are therefore
absent from the paper.

**(A)** Why no covariate of the baseline is adjusted on the z scale. Median z of
each of 105 external GO sets against the median shift of its baseline
expression; red marks sets flagged as baseline artefacts. The per-set
correlation is r = −0.112, but the rule rests on the gene-level correlation
between baseline expression and dWT, median −0.566 across all eleven data points
and negative at every one of them: adjusting on this scale destroys the positive
control (plasma-cell machinery z +9.10 to −2.43) while making 10 of 11 points
null-intact. The rule is general over baseline covariates, not specific to
expression.

**(B)** What offset the metric actually finds. An additive offset injected into
gene sets of 8 to 60 genes per side is recovered as 0.35 z at 60 genes per side,
0.50 at 30 and 0.75 at 15; every negative result in the paper therefore means no
effect of that magnitude. Colour: set size. Dashed line: z = 2.

**(C)** Self-test against the known null rate. 520 neutral contrasts give
|z| > 2 in 6.3 % of cases against a nominal 5 %: slightly anticonservative, and
stated wherever a z = 2 threshold is used.

**(D)** Why two convergence counts cannot be placed side by side. Convergent
genes against the concordance threshold for dWT (teal) and iv (orange); solid
observed, dashed the noise floor, band its 5th to 95th percentile. The iv curve
lies inside its own floor at every threshold, and compressing the shared
component of dWT by 35 % costs 91 % of the convergent genes, 120 down to 11.
Hence dWT serves only as a positive control and a tool, never as a result.

## Figure S2 · Every calibration in the study, in one place, Related to Figure 1

Confirmatory as bookkeeping. All 44 calibration attempts on one axis: the 18
perturbation datasets, the 14 donor cells, and the 12 combinations of ATAC axis
and window (8 of calibration type D, 4 of type L). Dashed line: the threshold
z = 2. Green passed, red failed, grey not calibratable because too few genes are
measurable. Each row is named on the y axis with its accession in parentheses
where one exists; the four rows without an accession are data first reported
with this study.

The 44 attempts are not 44 independent experiments. All six study units of the
donor level are also among the 18 datasets, and the four ATAC windows are four
window definitions of one dataset. The load-bearing head counts are therefore
the two study-level ones: 2 of 18 datasets and 7 of 14 donor cells.

## Figure S3 · External triangulation, all three versions, Related to Figure 2

Preregistered follow-up.

**(A)** The three synthesis versions computed for this level, with their
statistic, P value and n. Only the first has a detection limit, and only it is
reported in the main text: the fixed 173-gene set, pooled directional share
0.682 over 6 studies and 1 038 observations, z +2.84, P = 0.0015, MDE80 0.679.
The pooled observation-level version (z +2.81, P = 0.0027) and the study-level
synthesis (z +3.18, P = 0.0002) are shown for completeness and marked grey,
because neither carries a detection limit.

**(B)** Directional share per study, each against its own MDE80 (vertical bar).
No single study reaches its own limit. The external triangulation carries as a
synthesis and only as a synthesis; it is not a per-study replication and is not
described as one.

## Figure S4 · Patient against control: a level that carries a limit, not a result, Related to Figure 2

Exploratory, and downgraded from an earlier role in the analysis.

**(A)** Calibration of each of the 7 patient cohorts. Three cohorts carry no
marker set applicable to their tissue (whole blood, CD14+ monocytes and
cell-free amniotic-fluid mRNA) and are marked as not calibratable. Of the four
that can be calibrated, two are calibrated with the undifferentiated marker set,
which demonstrates tissue identity rather than position on the differentiation
axis; the decision gate that once rested on this level is therefore withdrawn.
Dashed line: the threshold z = 2.

**(B)** Programme concordance, patient against control, per cohort, against each
cohort's own MDE80 (vertical bar). None of the 7 cohorts reaches its own limit;
the limits run from 0.576 to 0.689. What this level supports is a statement
about the detection limit, valid under every leave-one-out, and not a statement
about an effect. Note also that this contrast, patient against control, is not
the contrast the module is defined on, which is undifferentiated against
differentiated.

## Figure S5 · Robustness, omission and the self-test, Related to Figure 2

Preregistered follow-up in panel A; methodological in panels B and C.

**(A)** Leave-one-out over the 7 calibrated donor cells. The programme runs
concordantly between donors: S1 = 0.349 against a null of 0.273 ± 0.025,
z +3.00, P = 0.0028, above an MDE80 of 0.344. Dropping any single cell leaves z
between +2.14 and +3.82, with P < 0.05 for every omission. Removing the cells
that contributed to defining the module raises the statistic to z +4.51, which
is the circularity control.

**(B)** Self-test of the metric in the donor set-up. Fraction of runs with
|z| > 2 against the injected effect, for each of the four statistics of that
set-up; dashed line marks the nominal 5 %. At zero injected effect the metric
returns 3.0 to 5.0 %.

**(C)** Sensitivity of the calibration balance to the marker set. The
preregistered marker sets against the subset reachable in vitro, with *SOST*,
*DMP1*, *PHEX*, *MEPE* and *PTH1R* omitted. Dashed lines mark z = 2 on both axes. The
balance is unchanged, 2 of 18 under both marker lists: the head count does not
depend on the choice of marker list.

## Figure S6 · The orthogonal levels in full, including the ones that carry nothing, Related to Figure 2

Exploratory.

**(A)** ATAC (GSE332758), complete: four windows against three axes under both
nulls, background (grey) and H1 stratified by baseline accessibility (blue);
vertical bar, the corresponding MDE80. H1 is the harder null and the one quoted
in Figure 2. Under H1 the adipogenic axis lies above its limit in 3 of 4 windows
(background null: 2 of 4); the lineage contrast lies below its limit under both
nulls and fails calibration L, so it is not measurable rather than zero.

**(B)** H3K27ac (GSE129031), module concordance per window and cell line. The
calibration passes; concordance runs from 0.590 to 0.636 across all nine
chondrogenic-axis tests, z +2.97 to +3.63, each above its own limit.

**(C)** Promoter methylome 27K (GSE33896), a null result with a limit. Target
osteogenic axis 0.349 against a limit of 0.478 (z −0.46); myogenic lineage
control 0.532 against 0.720 (z −1.76). Both lie below their limits. A paired
Wilcoxon test on the same donors and the same 126 genes shows more movement on
the target axis (W = 2 774, P = 0.0028), but that comparison has no detection
limit of its own and is used nowhere.

**(D)** Promoter methylome 450K (GSE129266). The calibration passes (z −3.91,
correctly signed); module concordance K = 0.298 sits below its limit of
K = 0.312 (z +1.85, P = 0.11). Again a limit, not a result.

## Figure S7 · The screens in detail, Related to Figure 1

Confirmatory.

**(A)** Exclusion codes of the diagnosis screen. Code A2, no diagnosis as the
lesion axis, accounts for 77 series; code A1, no undifferentiated arm, for 46;
two further single-series codes and the 2 complete 2 × 2 designs make up the
rest. The reason recorded for each excluded series is listed in full in
Supplementary Table 2 rather than in the panel, where it would not fit at
legible size.

**(B)** The 127 hand-checked series by the disease axis searched for.

**(C)** The independent screen by design rather than by entity: 22 series with
their exclusion codes.

## Figure S8 · Day zero falls as a predictor, and the publication matching as a methods control, Related to Figure 2

Exploratory.

**(A)** The undifferentiated state does not predict later differentiation.
Day-0 module value against later amplitude, n = 7 calibrated cells: z −1.68,
P = 0.084; direction P = 0.50; raw scale z +0.83, P = 0.67. Dashed lines mark
|z| = 2.

**(B)** The same question from the other side: the 7 cells that later fail
their calibration do not differ at day 0 from the 7 that pass (U = 23.5,
P = 0.95). Day zero carries no competence information; the undifferentiated arm
remains necessary as a contrast, not as a predictor.

**(C)** The confounder, measured directly. Median publications per gene
relative to a background median of 73 (NCBI gene2pubmed, retrieved 2026-08-22).
The cell cycle, the negative control, is the most-studied set at 2.03 times
background, which defeated two candidate axes; the programme sits at 1.00 times
background and is untouched.

**(D)** What the third matching axis does. Absolute expression per set without
publication matching (red) and with it (blue): every set collapses, the
negative control furthest, z +5.01 to −3.05. By the decision rule fixed before
the run, an axis that pulls the negative control to the null measures study
intensity rather than biology, and absolute expression appears nowhere in this
paper as a result.

## Figure S9 · The lineage contrast across three chromatin cohorts, and an independent observation of the decoupling, Related to Figure 2

Exploratory throughout. Panels A and B follow up the ATAC lineage contrast
that Figure 2 reports as not measurable on a single cohort; panels C and D
report an unplanned observation from the same run, whose preregistered version
is Figure 2F.

**(A)** Calibration L, osteogenic minus adipogenic markers on the difference
axis, per window and cohort, against each row's own detection limit (vertical
bar). None of the 12 windows reaches it, so by the decision rule fixed before
the first number the lineage contrast stays out of the main figure.

**(B)** Why each cohort fails. The osteogenic marker set never sits above
z = +2: in ATAC (GSE332758) it does not move (z −0.47 to +1.66), in ATAC (GSE151311)
underpowered (MDE80 0.25 to 0.68 for a contrast near zero;
replicate coverage differs threefold, median 223 against 77), and in H3K27ac (GSE151315) the
contrast carries the wrong sign (−0.19 to −0.27). The question is not
answerable with the chromatin data currently public, and the limits state how
large an effect would have had to be.

**(C)** GSE151315, the one cohort with an undifferentiated arm as well as two
lineages, calibration against module per window. The calibration passes 0 of 8
while the module lies above its own detection limit in 8 of 8, z +3.10 to
+7.02.

**(D)** Where each marker set sits on each axis. The undifferentiated markers
fall (z −1.46 to −2.80), the osteogenic markers do not rise (z −1.16 to
−1.30), and on the osteogenic axis it is the adipogenic markers that rise
(z +1.70 to +2.51): the cultures differentiate without reaching their label,
and the module runs regardless. Calibration D itself is underpowered
(MDE80 0.29 to 0.55), so its failure means partly not reached and partly not
testable.

---

# Supplementary tables

**Table S1. The 18 perturbation datasets.** Accession, arm (osteogenic or
chondrogenic), design, sample counts and perturbation, one row per dataset.

**Table S2. Every screened series with its verdict and exclusion code.** Verdicts of the preregistered search screening, one row per series; the source column records where the candidate surfaced (preregistered search; priority or secondary list plus search).

**Table S3. Every calibration, per dataset and per donor cell.** Contrast, null
mean and standard deviation, z, P, MDE80, status and pass or fail, for all 32
study-level and donor-level attempts, in both the unfiltered and the filtered
gene space. Study level 2 of 18; donor level 7 of 14.

**Table S4. The 173 module genes.** Ensembl identifier, symbol, direction (+1 in
129 genes, −1 in 44), number of datasets contributing and median effect. This
set was fixed before every analysis in this paper.

**Table S5. In vivo, per zone and per specimen.** Cell count, module contrast
with its own MDE80 and z, and the positive control with its own MDE80 and z, for
each of the 66 (zone, sample) pairs, with the specimen identifier.

**Table S6. Gene panels and mechanism classes.** Which panel genes fall into
which of the seven mechanism classes, with counts, and the sources and versions
of the disease and trait gene sets used throughout.

**Table S7. Every statistic of the paper with its detection limit.** Module per
dataset, ATAC, H3K27ac, both methylomes, the gene-set matching, the mechanism
classes and the complementarity tests, each with its statistic, its null and its
MDE80. This is the table to read against any number in the text.

**Table S8. The preregistrations.** Every preregistration and protocol document
with type, date and title. They are reproduced unchanged in the code archive,
including the ones documenting hypotheses that fell.

**Table S9. Figure 2D in full: narrow sets, broad sets and a second independent
source.** All six categories against every set computed, with source,
version-bearing set name, set size in the background, k, odds ratio, 95 %
confidence interval, raw and Bonferroni-corrected P, and the verdict of the
pre-set decision rule. Three things are visible here that the main figure
compresses: the cell-cycle exit number that the decision rule removed from the
running text (broad OR 1.80 [0.93–3.20] against narrow 2.11 [1.09–3.77]); that
no category reaches Bonferroni significance in either pass; and that the second
independent source, Reactome, reproduces the direction of both matrix categories
without reaching the criterion (OR 1.55 [0.18–5.84] and 1.45 [0.39–3.84]).
**Table S9b** lists which module genes carry which category under which set, so
that k = 5 can be read as five named genes.

**Table S10. The three-way decomposition, all 18 datasets, both gene spaces.**
One row per dataset and gene space, with the undifferentiated-marker,
own-lineage-marker and other-axis-marker contrasts, the module concordance and
its own detection limit, the dataset's calibration, and the verdict of the
preregistered rule. The unfiltered block is the primary run, in the same gene
space as the calibration of Figure 1D; the filtered block is the sensitivity run
in the space of the internal gene map. The verdict is the same under both, 8
against 7 confirmations and 0 refutations either way, so the decomposition does
not depend on the gene-space decision.

**Table S11. Is the in-vivo trend carried by a handful of genes?** The
donor-stratified trend test repeated after removing the module genes with the
largest absolute zone-to-zone difference, each run against the reduced module's
own detection limit. The 10 % row is the pre-set decision rule; the 5 %, 20 % and
30 % rows were declared descriptive before the run and decide nothing. The trend
clears its limit at 10 % (ρ 0.144 against 0.125, z +2.74) and falls below it at
20 % and 30 %. The series is not monotone in the number of genes removed, because
the statistic is a rank correlation across samples, and that is shown rather
than smoothed. **Table S11b** is the ranking itself: all 173 module genes by
absolute difference, with each gene's own concordance flag. The 17 largest carry
57 % of the total.

**Table S12. Primary publication of every GEO series used.** Author, year,
journal, volume, pages and DOI, covering the 18 perturbation datasets, the
patient cohorts and the orthogonal chromatin and methylome series. Series
without a linked publication in GEO are marked as such rather than attributed to
a guess.

**Table S13. Does the in-vivo anchor hang on the single hypertrophic sample?**
The donor-stratified trend test repeated with every hypertrophic point removed,
each run against the reduced selection's own detection limit; the rule was fixed
before the run. Under the rule the module trend holds at ρ 0.430, z +4.65
against its own new limit of 0.250, and the positive control holds at ρ 0.880
against 0.615, so the axis stays calibrated without the terminal zone. A second
variant, declared descriptive in advance, also removes the prehypertrophic zone
and holds as well.

**Table S14. Levels, detection limits and what each level carries.** One row per
level (the same levels as Figure 6, in the same order), with its calibration,
its own detection limit, the result and the verdict as text columns. It carries
the calibration column that the figure leaves to colour, and the per-level
wording of limits that live in units other than the headline statistic (dWT
units for the dynamics axis, LOEUF units for the constraint contrast).

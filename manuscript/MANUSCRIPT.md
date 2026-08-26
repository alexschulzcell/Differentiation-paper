# What skeletal differentiation models actually measure: a lineage-independent matrix programme, and why the disease genes lie orthogonal to it

Alexander Schulz^1^ (ORCID 0009-0009-2605-4350), Christian T. Thiel^1,\*^ (ORCID 0000-0003-3817-7277)

^1^ Institute of Human Genetics, Universitätsklinikum Erlangen, Friedrich-Alexander-Universität Erlangen-Nürnberg, 91054 Erlangen, Germany

\* Lead contact and corresponding author: Christian T. Thiel
(Christian.Thiel@uk-erlangen.de)

## Summary

Skeletal dysplasias are genetically heterogeneous, and the field's standing
hope is that they converge on a shared downstream signal that mesenchymal
differentiation models could report. We asked what those models measure,
calibrating every level against its own positive control and its own measured
detection limit. Across 18 published perturbation datasets only 2
pass a calibration on the lineage markers of the axis they model, and none
reaches its own detection limit on it. Yet a fixed 173-gene matrix
programme runs above its own limit in 18 of 18, as strongly where the
calibration fails as where it passes, and in chromatin also on the adipogenic
lineage. It tracks the human fetal growth plate to the prehypertrophic zone.
The dysplasia genes lie orthogonal: defined by distal secretory localisation
and gene dosage, ordinary on the differentiation axis. The transcriptional
downstream marker these models are searched for is structurally not there.

## Keywords

skeletal dysplasia; short stature; mesenchymal differentiation; detection
limit; calibration; extracellular matrix programme; growth plate; gene
constraint; reanalysis

---

## Highlights

- Only 2 of 18 published differentiation datasets pass their own lineage
  calibration
- A fixed 173-gene matrix programme runs in 18 of 18, lineage-independently
- Skeletal dysplasia genes are ordinary on the differentiation axis, 0 of 40
  tests
- Every null reports its own detection limit; limits, not absences

---

## Introduction

Skeletal dysplasias and disorders of growth are among the most genetically
heterogeneous groups in clinical genetics. Several hundred genes are established
as causal [Unger et al. 2023], they span structural collagens, glycosylation and linker enzymes,
transcription factors, cilium components, signalling molecules and vesicular
transport, and the phenotypes they produce overlap heavily. The clinical
consequence is familiar: a diagnosis often follows the variant rather than the
presentation, and there are few handles that generalise across lesions.

The standing hope is that this heterogeneity resolves downstream. If many
lesions eventually disturb one shared programme, then that programme is a
biomarker, a stratification tool and possibly a target, and one would not need
a separate approach per gene. This is not a straw man: it has produced both
a licensed drug and an active drug-development line. Vosoritide, a CNP
analogue acting on the MAPK pathway *downstream* of *FGFR3*, has been licensed
since 2021 on the strength of a phase 2 trial [Duggan 2021; Savarirayan
et al. 2024, phase 2 trial], is now used in thousands of children with the
therapy monitored under international consensus guidelines [international
consensus guidelines on implementation and monitoring, 2024], and is
described in its own literature as potentially useful in any growth
disorder with increased MAPK signalling. Independently, a large group
of chondrodysplasias and collagenopathies is treated as protein-folding
disease: misfolded matrix protein is retained in the ER, the unfolded protein
response is engaged through *IRE1*, *ATF6* and *PERK*, and chondrocyte
differentiation is disturbed, with the explicit conclusion that relief of ER
stress is a lesion-crossing therapeutic avenue [Cameron et al. 2015; Cameron et al. 2011; Farhan et al. 2026].

It is the *transcriptional* version of this hope that motivates a large body of
work in which mesenchymal stromal cells carrying a disease-relevant
perturbation are differentiated towards bone or cartilage and profiled. That is
the version we test here.

The expectation is testable, and to our reading it has not been tested against
a noise expectation. Two things are needed for the test and are usually
absent. The first is a positive control internal to each dataset: a
demonstration that the culture in question actually reached the lineage it is
supposed to model, measured on the same scale and against the same null as the
quantity of interest. The second is a detection limit: for every negative
result, the size of effect the analysis would have found. Without the first, a
null result can mean the biology is absent or that the experiment did not
happen. Without the second, a null result carries no information at all.

We therefore built the test around those two requirements and applied it at
every level we could reach: transcriptome across 18 published perturbation
datasets and, within one of them, across donors; chromatin accessibility and
H3K27ac; two promoter methylomes; a human fetal growth plate atlas; and human
genetics through curated disease panels, gene constraint and mode of
inheritance. Each level carries its own positive control and its own measured
detection limit, and the levels that carry nothing are reported alongside the
ones that carry something.

Two questions organise the results. What do these models reliably measure?
And are the disease genes part of it? The answers turn out to be
independent of each other, and the second is a negative result of the kind that
is only worth reporting because a positive result in the same set-up shows the
set-up works.

---

## Results

### The calibration that almost no published experiment passes

Before asking what the 18 perturbation datasets measure, we asked whether they
reach the lineage they claim to model. The rule was fixed in advance and is reproduced unchanged in the code
archive (preregistration M-D, §6): each dataset is tested for the textbook markers of its
own axis (osteogenic or chondrogenic) against a null that matches background
genes on expression decile and union exon length, and passes at z ≥ 2.

In this preregistered test, two of eighteen datasets pass (Fig. 1D): *LAMA5*-KO chondrogenic (z +2.67)
and *ERCC6L2*-KD (z +2.13). One further dataset is not calibratable at all (it
retains a single marker of its axis), so counted strictly the figure is 2 of 17
calibratable datasets. Two more come close without reaching the threshold
(*ACVR1*/FOP z +1.90, *RB1*-mut isogenic z +1.93) and are counted as failures.
No dataset reaches its own MDE80 on this calibration (that limit sits at
z ≈ 2.8), so the two that pass do so on the preregistered z ≥ 2 rule alone, and
we report it that way rather than as a strong positive.

The calibration also says *where* these cultures stop, and this is the
one place in the paper where the positive control is itself informative. The
markers that separate the two passing datasets from the rest are the
terminal ones. In *RB1*-mut isogenic the early osteogenic regulators rise
(*RUNX2* +1.47, *ALPL* +1.88, *POSTN* +2.02, *COL1A1*/*COL1A2* +1.35/+1.51) while *BGLAP*,
*SOST* and *MEPE* sit at exactly zero and *DMP1* falls (−2.22); in *TP53* LFS,
*RUNX2* (+2.67) and *DLX5* (+1.97) rise while *SP7*, *SOST*, *DMP1* and *MEPE* all fall.
The cultures engage the early osteogenic programme and stop before
mineralisation. This is visible only because the marker sets are not filtered
for cross-dataset measurability, a decision explained in Methods, and it is the
same conclusion the decomposition of Figure 2F reaches by a different route.

Resolved to
individual donor cells, 7 of 14 pass (Fig. 1E); of the cells carrying a
genuine patient lesion rather than an engineered perturbation, exactly one does.
The two counts differ in resolution, not in material, and the load-bearing one
is the study-level count.

The material situation behind this head count matters in itself. Of 1 424 GEO
series screened for a disease diagnosis as the lesion axis, 127 were checked by
hand and 50 carry such an axis; 46 of those 50 (92 %) have no undifferentiated arm
(Fig. 1C). An independent screen by design rather than by entity found 22 series
and 2 complete 2 × 2 designs, and both are already among the 18. A missing
data type is an observation, not a finding, and we treat it as one. But it
explains why the calibration is so rarely checked: without an undifferentiated arm there is
nothing to calibrate against.

This raises the question the rest of the paper answers. If most of these
cultures do not demonstrably reach their lineage, what is it that they measure?

### One matrix programme, running independently of lineage

A fixed 173-gene programme, defined before any analysis reported here
(Supplementary Table 4), runs above its own detection limit in 18 of 18
datasets, z +5.25 to +13.10 (Fig. 2B). Each point is tested against its own
null and its own MDE80; not one is a borderline case.

The programme is not lineage-specific. Cross-arm concordance of the module genes
between the osteogenic (12 datasets) and the chondrogenic (6 datasets) arm is
ρ = +0.622 against a null of 0.091 ± 0.076, z = +7.03, P = 0.0004
(Fig. 2A). A programme specific to bone or to cartilage would diverge between
the arms.

Nor is it coupled to lineage commitment. Pooled by calibration
status, the datasets that fail their calibration give concordance 1.000,
z +13.13 (n = 16, MDE80 0.617); those that pass give 0.994, z +12.79
(n = 2, MDE80 0.612) (Fig. 2C). There is no difference. The calibration
observation is therefore not a methodological complaint about the field: it is
the evidence that two separable layers exist, and that published experiments
routinely capture one without the other.

What the programme is made of we tested against external gene sets rather than
sets derived from the result. Because the first pass used single, narrow GO
terms whose backgrounds are 33 to 104 genes, we repeated it against broad,
independently curated sets (the matrisome for the two matrix categories,
Reactome for the signalling and senescence categories), fixed with source,
version and retrieval date before the run, under a decision rule written into
the script head before the first number (`code/24_gene_sets_v2.R`). Both passes
are reported side by side (Fig. 2D).

The direction survives and the intervals tighten: matrix remodelling OR 3.62
[1.14–8.91] against 265 matrisome ECM regulators (narrow GO: 7.92 [1.53–25.9]),
matrix components OR 3.38 [1.20–7.73] against the core matrisome (narrow:
4.03 [1.26–9.95]), and, in the opposite direction, secretory machinery
OR 0.32 [0.10–0.78], depleted; this is the one category for which no broader set
exists, so narrow and broad are the same test. TGFβ/BMP and hypoxia/stress are
null in both passes. Cell-cycle exit does not survive: it keeps its
direction but its interval now covers 1 (OR 1.80 [0.93–3.20]), and by the
pre-set rule we no longer carry its number in the text (Supplementary Table 9).

The programme is therefore matrix output, depleted for the machinery that
secretes matrix; the cell-cycle component is a direction we can state and not
a magnitude we can defend. That depletion turns out to matter (below). Three
things limit this result and are stated rather than footnoted: k remains 5 to
13 genes per category, because only 147 of the 173 module genes lie in the
background of 11 581 at all; no category reaches Bonferroni significance in
either pass (smallest p_Bonf 0.042), the criterion applied being the pre-set
one that the 95 % CI excludes 1; and a second independent source for the
two matrix categories (Reactome ECM organisation, degradation and collagen
formation) reproduces the direction but not the criterion (OR 1.55 [0.18–5.84]
and 1.45 [0.39–3.84], Supplementary Table 9). This is a claim about what the programme is
*made of* in direction, and it should not be read as a quantitative estimate.

The same question can be put to a different assay. GSE332758 carries an osteogenic
and an adipogenic differentiation axis in the same cells, measured by ATAC-seq,
in four window definitions (Fig. 2E, F). The adipogenic axis passes its own
calibration in all four windows (z +3.73 to +5.01) and the module runs there
above its own limit in 3 of 4 windows under the harder, baseline-stratified
null (2 of 4 under the background null), z +2.38 to +4.76. The module runs on a
non-skeletal lineage, in an assay that measures accessibility rather than
transcript abundance. That is the lineage independence, positive and
calibrated.

The osteogenic axis of the same dataset fails its calibration in all four
windows (z −0.54 to +1.66) while the module runs in all four (z +3.50 to +4.51).
This is the decoupling in its sharpest single-dataset form, but because the
calibration fails we report it as an observation about decoupling and not as a
calibrated module result.

We also tested the lineage contrast itself (osteogenic minus adipogenic), where
the module is absent in all four windows. We do not report that as a
finding. The lineage axis fails its own calibration (type L) in every window
(z +0.93 to +2.26), which means the assay does not demonstrably separate the two
lineages there. A zero on an axis whose positive control fails is *not
measurable*, not *no effect*, and that distinction is the one this analysis is
built around.

Because that verdict rested on one dataset, we obtained two further
chromatin cohorts carrying both an adipogenic and an osteogenic arm from the
same cells, ATAC-seq (GSE151311) and H3K27ac (GSE151315), and applied the
same test under a decision rule fixed before the first number
(`code/23_lineage_contrast.py`). Calibration L fails in 0 of 12 windows across
all three cohorts, for three different reasons: in GSE332758 the osteogenic
markers do not move; in GSE151311 the test is severely underpowered
(MDE80 0.25–0.68 for a contrast of ≈ 0, two replicates differing threefold in
coverage); in GSE151315 the contrast carries the *wrong sign*. The lineage
contrast is therefore not answerable with the chromatin data currently public,
and the limits state how large an effect would have to be (Fig. S9A, S9B). We
report that as a detection limit for the level, not as a null result.

The same run produced something we did not go looking for.
GSE151315 is the only chromatin cohort here that carries an undifferentiated arm *as
well as* two lineages. That makes the calibration decomposable: instead of
one pass/fail, one can ask separately whether the culture left the undifferentiated
state and whether it arrived at its lineage. Both questions have an
answer, and they differ. The undifferentiated markers go down (z −1.46 to −2.80 on the
osteogenic axis, −1.11 to −2.77 on the adipogenic): the cells do leave the
undifferentiated state. The osteogenic markers do not go up (z −1.16 to −1.30); on the
osteogenic axis it is the adipogenic markers that rise (z +1.70 to +2.51).
And the module runs regardless: calibration passes 0 of 8, the module is above
its own detection limit in 8 of 8, z +3.10 to +7.02 (Fig. S9C, S9D).

This observation was unplanned, and we did not want to build on it in that
state. We therefore wrote it down as a hypothesis, with its decision rule,
before computing it anywhere else (preregistration F2F,
`preregistrations/`). The hypothesis says that the undifferentiated markers
fall, at z ≤ −2 against their own null, that the own-lineage markers do not
rise above z = +2, and that the module lies above its own detection limit.
The decomposition is confirmed when the first and the third of these hold
while the second does not. When instead the second and the third hold, the
lineage is reached and the module runs, which the rule states in advance is
the *other* case of decoupling and not a contradiction. It is refuted when
the cells leave the undifferentiated state and the module fails to run. Preregistering it afterwards
does not make the chromatin observation confirmatory (it is not, and it keeps
its label); it makes the replication confirmatory.

The preregistered replication needed no new data. All 18 perturbation datasets
carry an undifferentiated arm, so the same decomposition is computable in every one of
them, with the calibration data the paper already reports (Fig. 2F). The result
is the pattern the model predicts, point for point. The undifferentiated
markers fall in 10 of the 18 datasets; the other 8 are reported as not
decomposable and carry no verdict in either direction. The own-lineage
markers rise above z = +2 in 2 of 18. This criterion is the same statistic
at the same threshold as the calibration, so its agreement with the
calibration status is definitional and we do not present it as a result;
it is carried because the decomposition needs the middle step, and because
it fixes what a failed calibration means: not that nothing happened, but
that the lineage was not reached. And the module lies above its own
detection limit in all 18 datasets.

Applying the pre-set rule to the 10 decomposable datasets gives 8
confirmations, 2 instances of the other case of decoupling, and 0 refutations
(Wilson 0.49–0.94 for the confirmed share). There is no dataset in which the
cells leave the undifferentiated state and the programme fails to run.

The cultures differentiate. They do not become what the label says. The
matrix programme runs anyway. That is the claim of this paper stated as a
three-step decomposition rather than as a correlation, now on 18 transcriptome
datasets under a rule fixed in advance, with the chromatin cohort where it was
first seen kept separate and kept exploratory.

Three qualifications travel with it. The 18 are not 18 independent
experiments (six of the study units reappear at donor resolution), so no
pooled statistic and no mean over the 18 z values is computed; the result is a
count with a Wilson interval and every point is listed individually
(Supplementary Table 10). The module test in each dataset is licensed by no
axis calibration; it is the *pattern* that is the claim, not any individual
z. And the gene space matters: the primary run uses the unfiltered marker sets
and background, the same gene space in which the calibration of Figure 1D is
computed. In the filtered space of the internal gene map the same rule gives 7
confirmations rather than 8 and, again, 0 refutations. Both runs are reported
per dataset in Supplementary Table 10.

### Not a culture artefact: the human fetal growth plate

If the programme were an artefact of plastic, serum and dexamethasone,
it should have no counterpart in a human embryo. It does. In a human fetal limb
single-cell atlas (Pcw 5.1–9.3), along the chondrogenic axis
MesCond → ChondroProg → Resting → Prolif → Prehypertrophic → Hypertrophic.

The positive control passes clearly. Chondrogenic minus undifferentiated markers per
sample, 64 of 66 samples above their own detection limit (65 of 66 at z > 2),
median z +10.92 (Fig. 3B). This is the most cleanly passed calibration of
any orthogonal level in the project.

The module then tracks the same axis: Spearman ρ = +0.456, z = +4.80,
P = 5 × 10⁻⁵, at a detection limit of ρ = 0.274 (Fig. 3C). The unit is the
specimen: 66 evaluable (zone, sample) pairs from 16 specimens across
9 developmental stages (Supplementary Table 5). The null permutes zone labels within the
specimen, so no difference in specimen means can produce the trend. The zone medians run 0.062 → 0.046 → 0.055 → 0.065 → 0.126 → 0.149.

Two limitations belong next to that number rather than in a discussion section.
First, the effect is quantitatively modest: about five times smaller than
the positive control on the same axis. Second, `HyperChon` contributes a single
evaluable sample of 8 cells (14 in the whole atlas). The trend carries to the
prehypertrophic zone and no further (Fig. 3D), and nothing in this paper rests
on the hypertrophic end point, which we verified rather than asserted: under a
rule fixed before the run, dropping every `HyperChon` point leaves the module
trend at ρ 0.430, z +4.65 against its own recomputed limit of ρ 0.250,
and leaves the positive control above its own new limit as well
(`code/34_hypertrophic_zone_sensitivity.py`, Supplementary Table 13). A gene-by-gene breakdown is also less clean than
the aggregate: only 76 of 164 module genes (46.3 %, 95 % CI 38.7–54.0 %) are
individually concordant, an interval that covers chance.

That leaves an obvious question (whether the aggregate trend is the work of a
minority of large, aligned excursions), and we asked it under a rule fixed
before the run (`code/27_in_vivo_gene_decomposition.py`): remove the 10 % of module
genes with the largest |Δ| and repeat the identical test against the reduced
module's own detection limit. The trend survives, and narrowly: without its
18 largest genes it is ρ 0.144, z +2.74, above its own limit of
ρ 0.125. By the pre-set rule the trend is therefore *not* the work of a few
genes. The rest of the picture qualifies that sentence and is reported with it:
ρ falls threefold, the 17 largest genes carry 57 % of the total |Δ|, and in
descriptive runs declared descriptive in advance, removing 20 % or 30 % of the
genes puts the trend below its own limit. What we claim is that a handful of
outlier genes does not produce this trend. We do not claim that it is spread
evenly across the module, and we make no statement about which parts of the
programme are in-vivo-confirmed (Supplementary Tables 11, 11b).

### Where the disease genes actually are

The set-up's preregistered positive controls pass first (Fig. 4A). Lineage markers
inside the disease panels give OR 17.0 to 51.6 (z up to +13.9), and the
secretion anchor separates distal from biosynthetic secretion at OR 3.70 to
5.75, P down to 2.9 × 10⁻¹¹. The anchor is a positive control precisely
because it is already known: COPII defects cause skeletal disease through
failed collagen export (*SEC23A* in cranio-lenticulo-sutural dysplasia, *SEC24D* in
zebrafish skeletal morphogenesis) [Sarmah et al. 2010; Tang and Ginsburg 2023; Claeys et al. 2021], Golgi tethering defects cause
Saul-Wilson syndrome (*COG4*), and congenital disorders of glycosylation carry
skeletal involvement [Lipiński et al. 2021]. We use it to show the machinery detects enrichment when enrichment
is there, not as a finding of our own.

**Axis 1 is localisation.** Skeletal dysplasia genes are enriched in the
distal secretory compartment (PanelApp 309: OR 2.84, z +6.30,
P = 1 × 10⁻⁴, limit OR 1.79) and depleted in the biosynthetic one (OR 0.54,
z −3.18), and the finding repeats in two independent nosology-derived panels
(Fig. 4B). The programme is in neither: distal OR 0.81, z −0.41, n.s. at a
limit of OR 2.34. The programme is depleted for secretory machinery (Fig. 2D);
the disease genes are enriched in it. The two sets occupy complementary
compartments of the same cell.

**Axis 2 is gene dosage.** Splitting PanelApp 309 by mode of inheritance
gives a textbook result: purely monoallelic genes are strongly
haploinsufficiency-constrained (LOEUF 0.283, n = 120), purely biallelic
genes are not (0.826, n = 247), P = 6 × 10⁻²⁴ (Fig. 4C). (LOEUF is
measurable in 247 of the 249 biallelic genes; |dWT| in all 249.) The result
holds within every publication tertile (P 0.017 / 2.2 × 10⁻⁵ /
1.7 × 10⁻¹¹), so it is not an artefact of study intensity, a control that
matters here, because monoallelic panel genes are far better studied than
biallelic ones (median 313.5 versus 84 publications, P = 2 × 10⁻²⁰).

**The same calculation yields a second contrast** (Fig. 4D). After matching on
expression, exon length and publication count, common height variation sits
in constrained genes (height GWAS z −4.67, limit 0.024 LOEUF units) while
monogenic dysplasia genes do not (PanelApp 309 z +1.83, nosology core +1.30,
both n.s. at limits of 0.095 and 0.139); if anything, in the opposite direction.
The architecture of common size variation and that of monogenic dysplasia differ
in gene constraint, and the difference survives the publication control. The
negative control collapses under that control (z −6.99 → −2.54), which is
precisely what shows the control is doing its job while the GWAS signal
survives.

**The remaining axis does not separate them.** On differentiation
dynamics (|dWT| at equal absolute expression, exon length and publication
count), the disease genes are ordinary: PanelApp 309 z −0.71 at a
detection limit of 0.073 |dWT| units; height GWAS z −0.81 at a limit of
0.016 (Fig. 4E). In the same run, the programme sits at z +18.10 and
the cell-cycle negative control at z +3.72. The same genes split by mode of
inheritance, the split that halves them on the dosage axis, differ by
P = 0.86 on the dynamics axis, and by P 0.90 / 0.77 / 0.99 within the three
publication tertiles.

**The closing negative result concerns mechanism classes.** Seven mechanism classes × seven panels
gives 49 tests, of which 40 are computable; none exceeds Bonferroni
(α = 0.00125) (Fig. 4F). Four tests exceed their own detection limit: all four
in the glycosylation/linker class, OR 6.45–8.59, nominal P 0.042–0.076, and
effectively two independent values because two panel pairs are identical. At 40
tests on a nominal 5 % level that is the chance expectation, and we report it
rather than omit it. Pooled across classes, the sharpest axis (height GWAS,
5 649 genes) gives OR 1.00 at a detection limit of OR 1.59, and on none of
the seven panels does the pooled odds ratio exceed its own limit.

This negative result is worth reporting only because the same calculation, on
the same genes, with the same null, finds the lineage markers at OR 17–52 and
the secretion anchor at OR 3.7–5.8. The set-up carries. The disease genes
are simply not on this layer.

An earlier version of this idea, that the disease genes are *constitutively
expressed infrastructure* rather than dynamic, was tested and fell
(z −0.35 at a limit of 0.069). The correct statement is not that they sit at the
low end of the dynamics axis. It is that they are not on it.

### Both layers meet at the prehypertrophic transition

The two layers are orthogonal in what defines them, but they are not in
different places. Contrasted against the module, per zone of the growth plate,
the disease genes rise along the same axis and peak in the same zone: 0.075 →
0.112 → 0.123 → 0.137 → 0.168 → 0.132, with median z +2.66 → +3.79 → +4.09
→ +4.42 → +4.98 → +4.17 (Fig. 5A). Both curves have their maximum in the
prehypertrophic zone (Fig. 5B).

This is a statement about spatial co-localisation, not about enrichment
(the enrichment test is the one that fell), and not about direction of
regulation. Human genetics and developmental biology point at the same place in
the growth plate, and it is a place that can be sampled.

### What each level carries, and what it does not

Figure 6 lists every level with its own calibration, its
own detection limit and its verdict, including the levels that carry nothing:
both promoter methylomes (27K null on both axes at z −0.46 / −1.76, both below
their own limits; 450K K = 0.298 below a limit of K = 0.312) and the undifferentiated state
as a predictor of later differentiation, which falls (z −1.68, P = 0.084;
the 7 cells that later fail their calibration do not differ at day 0 from the 7
that pass, P = 0.95).

The distinction the figure enforces is between a good negative result (a
null with a positive control in the same set-up and a stated limit) and a
not measurable (a null on an axis whose own calibration fails). The ATAC
lineage contrast is in the second category and is labelled as such.

---

## Discussion

### 1 · What these models measure, and what they do not

The single most consequential number in this paper is the pair 18/18 against
2/18. A fixed matrix programme runs above its own detection limit in every one
of eighteen published perturbation datasets, while only two of them
demonstrably reach the lineage identity they were built to model; neither
of those two reaches its own detection limit on that calibration. These are two
separable things, and in the datasets we could check they are routinely read as
one: a change in the shared programme is taken as evidence that the lineage
axis was perturbed. We put it that way, as a statement about the eighteen
datasets in front of us, rather than as a claim about the practice of the
field, which we have not measured.

The existence of a shared component is not itself new. Liu *et al.* compared
bone-marrow and adipose MSCs across osteogenic, adipogenic and chondrogenic
differentiation and concluded that a set of common genes is needed for early
differentiation into all three lineages [Liu et al. 2007]. What has been
missing is the second half: a null expectation and a detection limit for
that shared set, and a per-experiment positive control that says whether the
lineage was reached at all. Without those, a shared programme and a
lineage-specific one are not distinguishable in the data, and that is the gap
this paper addresses.

A second neighbour lies in the same direction and is cited and delimited here:
lineage assignment in these cultures is weakly determined by the induction
protocol and strongly by the matrix [He et al. 2017]. That
weakens the protocol as a guarantee of lineage identity, which is consistent
with what we measure, but it does not test a fixed shared programme against a
noise expectation, and it reports no detection limit.

We want to be careful about what this does and does not say about the field. It
does not say that these experiments are wrong or that their conclusions are
void; many of them ask questions for which lineage identity is not the point. It
says that when such a dataset is used to ask whether a lesion perturbs
*osteogenesis* or *chondrogenesis*, the reader currently has no way of knowing
whether the culture got there; in the datasets we could check, it usually
did not. The programme that *is* reliably measured is shared between
lineages and is not, by itself, evidence of lineage-specific biology.

### 2 · The separation of the two axes, and why the calibration proves it

The strongest single piece of evidence for the separation is not a large effect
size but a non-difference in the right place: pooled by calibration status,
failed (z +13.13) and passed (z +12.79) are indistinguishable. If the module
were a downstream consequence of lineage commitment, it should be weaker where
commitment demonstrably did not occur. It is not weaker at all.

The chromatin data make the same point on a foreign assay and in a single
dataset: the module runs on the adipogenic axis (a lineage that is not
skeletal) above its own limit, on an axis that passes its own calibration in
all four window definitions.

The lineage contrast was given a second chance, and it did not take it.
Two further chromatin cohorts were obtained specifically to test the elegant
version of the argument below. Calibration L fails in 0 of 12 windows across
three cohorts, each for a different reason (Fig. S9A, S9B). We take that
seriously as an answer: with the chromatin data now public this question cannot
be asked, and the useful output is the limit rather than a number. A cohort
with more than two replicates per state and a demonstrably separable pair of
lineages would settle it.

**The attempt paid for itself in a different currency.** The cohort obtained to
test the lineage contrast turned out to be the only chromatin dataset in this
study with an undifferentiated arm *and* two lineages, and it converts the decoupling from
a correlational statement into a decomposition: the cells leave the undifferentiated
state, they do not arrive at the lineage, and the module runs anyway
(Fig. S9C, S9D). That observation was unplanned, and rather than build on it in
that state we wrote it down as a hypothesis with its decision rule and then
tested it where it could be tested confirmatorily: on all 18 transcriptome
datasets, every one of which carries an undifferentiated arm (Fig. 2F). It
replicates: 8 confirmations, 2 instances of the other case of decoupling,
0 refutations, and the two datasets that reach their lineage are precisely the
two that pass their calibration.

We consider the decomposition the strongest single argument in the paper for
the central claim, and we are explicit about how it earned that status: the
chromatin cohort where it was first seen remains exploratory and stays in the
supplement, and what is preregistered is the replication, not the original
observation. The replication is also not an independent assay: it is the
same transcriptome material the rest of the paper uses, re-read along a
different question. What it adds is that the pattern is not an artefact of one
cohort or one measurement technique; what it does not add is a second
laboratory.

We deliberately do not use the elegant version of this argument. Subtracting
the adipogenic from the osteogenic axis leaves nothing of the module in any
window, which reads as a clean positive-logic demonstration that the module *is*
the shared component of both lineages. But the lineage-contrast axis fails its
own calibration in all four windows: in that dataset the two lineages are not
demonstrably separated at the level of accessibility, so the subtraction has
nothing to subtract. The result is *not measurable*, and presenting it as a
finding would be precisely the error the calibration rule exists to prevent.
A properly calibrated lineage contrast (an ATAC dataset in which osteogenic and
adipogenic markers are cleanly separable) would test the claim directly, and
that experiment is worth doing.

### 3 · The developmental anchor and its limit

The programme tracks the maturation axis of the human fetal growth plate to the
prehypertrophic zone, with a positive control that is the cleanest in the
project. That rules out the most obvious deflationary reading: that it is a
culture stress response.

The limits are real and stated rather than deplored. The effect is about five
times smaller than the positive control on the same axis. Only 46.3 % of module
genes are individually concordant, an interval covering chance, so we make no
claim about which parts of the programme are in-vivo-confirmed. We did test
whether the trend is the work of a minority of large excursions, under a rule
fixed before the run, and it is not: it survives the removal of the 10 % of
genes with the largest |Δ| (ρ 0.144 against its own limit of 0.125). But it
survives narrowly, ρ falls by a factor of three, and removing 20 % of the genes
takes it below its limit. The trend is broader than a handful of genes and
narrower than the whole module, and that is as far as this dataset goes. `HyperChon`
contributes 8 cells in one sample, so nothing rests on the hypertrophic end
point. And this is one atlas: a second independent in-vivo source, ideally a
postnatal growth plate rather than a prenatal limb bud, would move this from
"one dataset" to "reproduced", and its absence is why this section is exploratory
rather than confirmatory.

### 4 · Why there is no transcriptional downstream convergence to find

Putting the layers together gives a structural rather than an empirical answer
to the field's question. The disease genes are defined by where their products
work (distal secretion, OR 2.84) and by how much of them is needed
(LOEUF 0.283 versus 0.826), and by differentiation dynamics not at all, at a
detection limit of 0.073 |dWT| units in a run where the programme sits at
+18.10. The programme, meanwhile, is *depleted* for secretory machinery. These
are complementary compartments, not two depths of the same one.

**What the broader gene sets did and did not change here.** The depletion that
carries this argument is the one category with no broader set available: the GO
secretory background is already 1 127 genes, and OR 0.32 [0.10–0.78] is the
same number in both passes. The matrix side, where broader and independently
curated sets do exist, keeps its direction with a tighter interval (OR 3.38 and
3.62 against the matrisome). What did not survive the broader sets is the
cell-cycle component: its interval now covers 1, and we have withdrawn its
number accordingly. That withdrawal costs this section nothing, because the
complementarity argument never ran through the cell cycle: it runs through
*secretion against matrix*. But it does narrow the description of the programme
from "matrix output plus cell-cycle exit" to matrix output, with cell-cycle
exit as a direction we can name and not a quantity we can defend. A reader who
wants the strong form of the claim (that the disease genes sit on precisely
the layer the programme is depleted for) should note that it rests on a single
category with k = 5 and on a confidence interval, not on a Bonferroni-surviving
test; no category in either pass reaches Bonferroni significance.

If that picture is right, then a transcriptional marker shared across lesions is
not merely hard to find in MSC differentiation models: it is not the kind of
thing these models measure. A gene whose product mislocalises, misfolds or
moves too slowly through the secretory pathway need not change its own
transcript level, nor the transcript levels of the matrix programme, in any
consistent direction. What changes is throughput.

And here the literature has already arrived, from the other side. The
ER-stress/UPR account of the collagenopathies and chondrodysplasias is a
lesion-crossing convergence claim, and it sits on exactly the layer our data
point to: protein folding and secretion, not transcript abundance of the
matrix programme. Our two axes say the same thing in genetic coordinates: the
disease genes are enriched in distal secretion (OR 2.84) and the programme is
*depleted* for secretory machinery (OR 0.32). We therefore do not claim the
field is looking in the wrong place. We claim something narrower and, we think,
more useful: the part of the field that looks at folding and secretion is on
the right layer, and the transcriptional readout of MSC differentiation models
cannot reach the same answer, not because the effect is small but because
that readout does not measure that layer. Our detection limits say how small an
effect we would have seen if it were there (0.016 to 0.113 |dWT| units).

The practical redirection follows from positive results rather than from a row
of nulls: measure the prehypertrophic transition, a place both layers point
at and one that can be sampled in human tissue, and measure secretory
throughput rather than transcript abundance.

### 5 · The methodological contribution

Three practices carried this analysis and we think they generalise beyond it.

Calibrate every differentiation experiment against its own lineage markers,
on the same scale and against the same null as the quantity of interest. In
our hands this reclassified 15 of 18 published datasets and, more importantly,
turned what looked like a quality problem into the observation that carries
the paper.

**Control for study intensity when comparing gene sets.** Two candidate axes in
this work, absolute expression and pooled constraint, were defeated by a
negative control that turned out to be the most-studied gene set of all
(2.03 × background). Adding publication count per gene as a third matching axis
decided both cases: absolute expression is an artefact and is reported nowhere;
constraint became interpretable and produced a result. The programme itself sits
at 1.00 × background, which is what secures its z +18.10 against the same
objection.

**State a detection limit for every null result.** Our metric recovers an
additive offset of 0.35 z at 60 genes per side, 0.50 at 30 and 0.75 at 15, at a
measured null rate of 3.0–5.0 % (donor set-up) and 6.3 % across 520 neutral
contrasts, against a nominal 5 %. Every negative result in this paper therefore
means "no effect of that magnitude", and where the limit is unreachable the
finding is called *not measurable*, never *no effect*.

### 6 · The wet-lab prediction

The picture makes a directed, falsifiable prediction, and it follows from the positive results.

> If the disease genes sit on the secretory layer and the programme on the
> matrix-output layer, then patient cells should show reduced secretory
> throughput while the matrix programme runs normally, measurable as
> retention in the ER and not as a transcriptome shift.

The experiment that would test it: an isogenic series of several patient lesions
in one laboratory, with a co-sequenced undifferentiated arm, a pulse-chase secretion
measurement, and the module score computed per donor. If patient cells show a
normal module score and normal secretory throughput, the prediction fails. If
they show a reduced module score, the orthogonality claim fails. Either outcome
is informative, which is the point.

The undifferentiated arm is not a formality. Its absence in 92 % of the diagnosis cohorts
we screened is the single practical obstacle to answering this question with
public data, and it costs almost nothing to include.

---

## Limitations of the study

The main result is exploratory. These analyses were carried out after the
preregistered questions had been closed, and we make no preregistered claims for
them.

The categories that define what the programme *is* rest on k = 5–13 genes each
(Fig. 2D), and that is still the weakest load-bearing point in the paper. It is
now tested twice, against narrow GO terms and against broad, independently
curated sets fixed before the run, which halves the width of the matrix
intervals but cannot raise k, because only 147 of the 173 module genes lie in
the measurable background at all. No category reaches Bonferroni significance
in either pass, and a second independent source (Reactome) reproduces the
direction of the two matrix categories without reaching the criterion. One
category, cell-cycle exit, did not survive the broader sets and its number has
been withdrawn from the text.

The in-vivo anchor is single-source, modest in size, and not carried by a
majority of individually concordant genes. It does not depend on the single
hypertrophic sample (Supplementary Table 13). We also attempted a second,
postnatal source and could not use it: in a pubertal human growth plate the
textbook positive control fails, because the hypertrophic zone switches the
cartilage-matrix programme off rather than intensifying it, so that level is
reported as not calibratable and changes nothing here
(Supplementary Table 14). It survives the removal of its
largest-|Δ| genes at the preregistered 10 % level and fails at 20 %; we
therefore describe it as broader than a few genes and not as uniformly
distributed.

The external triangulation, the preregistered pooling of the programme effect
across the published cohorts (Fig. S3), carries only as a synthesis: no single one of the
eleven contributing studies reaches its own detection limit (Fig. S3B). We
report z +2.84, P = 0.0015 at MDE80 0.679 and describe it as what it is.

Several levels carry nothing and are reported anyway: both promoter methylomes,
the patient-versus-control comparison (downgraded, gate B withdrawn), and the
undifferentiated state as a predictor.

The versions of this question that fell are reported in full. Three earlier
framings failed against their own controls, each preregistered before the data
were seen: that differentiation converges while lesion responses do not
(the convergence count is effectively a step function of the signal-to-noise
ratio: compressing the shared component by 35 % costs 91 % of convergent genes,
120 → 11); that patients diverge further still (a different contrast was tested,
and two of four cohorts were calibrated with the wrong marker set); and that
day zero carries competence information (P = 0.95 between cells that later pass
and fail). The preregistrations are reproduced unchanged in
`preregistrations/`, including the ones documenting what fell, and the current
version of the question is stronger because the weaker ones are on the record.

---

## Acknowledgments

We thank the authors of the eighteen perturbation studies, the two single-cell
atlases and the chromatin and methylome series reanalysed here for depositing
their data publicly. This paper exists because they did.

This work was supported by the Deutsche Forschungsgemeinschaft (DFG, German
Research Foundation) grant TH896/7-1 (C.T.T.).

## Author Contributions

Conceptualization, A.S.; Methodology, A.S.; Software, A.S.; Validation, A.S.;
Formal Analysis, A.S.; Investigation, A.S.; Data Curation, A.S.; Writing –
Original Draft, A.S.; Visualization, A.S.; Writing – Review & Editing, A.S. and
C.T.T.; Supervision, C.T.T.; Project Administration, C.T.T.; Funding
Acquisition, C.T.T.

## Declaration of Interests

The authors declare no competing interests.

## Declaration of Generative AI and AI-assisted technologies in the writing process

During the preparation of this work, the authors used Claude Opus 5
(Anthropic) as the principal tool for language editing, for drafting sections
of the manuscript, and for implementing, documenting and checking the analysis
and figure code. Output from other large language models was used
occasionally during earlier stages and was revised before use. After using
these tools, the authors reviewed and edited all content as needed and take
full responsibility for the content of the publication. All scientific
content, study design, data analysis, interpretation, and conclusions are the
sole responsibility of the authors.

---

## STAR Methods

### Resource Availability

#### Lead contact

Requests for further information should be directed to the lead contact,
Christian T. Thiel (Christian.Thiel@uk-erlangen.de).

#### Materials availability

This study did not generate new unique reagents. It is a reanalysis of
publicly available data together with data first reported here; all materials
underlying it are the datasets listed in Supplementary Table 1 and in the Key
Resources Table.

#### Ethics

This study generated no new human material and no new human data. It analyses
exclusively datasets that are already publicly deposited; ethical oversight
and participant consent for each of them are documented in the primary
publication of that dataset, which is cited in the reference list and mapped
to its accession in Supplementary Table 12. No separate institutional review
board approval was therefore required for the work reported here.

One consequence should be stated explicitly: the
urine-derived stem cell series (E-MTAB-16566) derives from a single donor, so
the influence of sex or gender on the quantities computed from it cannot be
assessed. The remaining datasets carry the donor structure their original
authors reported, which Supplementary Table 1 lists per series.

#### Data and code availability

All data reported in this paper are public. The 18 perturbation series, the
patient cohorts, the chromatin and methylome series and the two single-cell
atlases are listed with their accession numbers in Supplementary Table 1 and
in the Key Resources Table; no dataset used here is unpublished. Sixteen of
the eighteen perturbation datasets come from fourteen Gene Expression Omnibus
series, and the remaining two are the chondrogenic and osteogenic arms of the
authors' own urine-derived stem cell series, deposited at ArrayExpress under
E-MTAB-16566.
All original code is publicly available as of the date of publication in the
GitHub repository (https://github.com/alexschulzcell/Differentiation-paper). The archive holds `code/`, the analysis
and figure pipeline, with a per-script header stating purpose, inputs, outputs
and runtime; `reference_implementations/`, one implementation of each metric;
`figure_style/`, the publication style of the figures; `figures/data/`, one
CSV per figure panel; `results/`, logs, session information and self-tests;
and `preregistrations/`, every preregistration unchanged, including the ones
that document questions that fell. Every number quoted in this manuscript can
be reproduced from the delivered files by `python code/70_check_numbers.py`,
which checks each value against the panel data and stops on the first
mismatch, and the reference apparatus is cross-checked in both directions by
`python code/71_check_references.py`. Any additional information required to
reanalyse the data reported here is available from the lead contact on
request.

### Additional resources

This study is preregistered. The preregistration and protocol documents,
dated before the computations they license, are reproduced unchanged in the
code archive under `preregistrations/` and listed with type, date and title in
Supplementary Table 8. They include the documents recording the three earlier
framings of this question that fell against their own controls. No clinical
trial is associated with this work.


### Reference implementations

**One implementation per metric, and no second one.** The module concordance,
the contrast statistic and the MDE80 are implemented once each:
`reference_implementations/_module.py` for everything computed here, and
`reference_implementations/manuscript/methods/03_metric.R` for Figures 1, 2
and S1-S5. Panel data files are produced by `code/50_panel_data.py` and
`code/51_supplement_data.py`, which read those outputs and reshape them without
recomputing anything. The single exception is `code/20_in_vivo_donor_test.py`
(see *In vivo*), which implements the donor-stratified trend test.

### Data provenance

All primary data are public. Accession numbers, arm, design and sample counts
for the 18 perturbation datasets are in Supplementary Table 1; the verdicts of
the preregistered search screening, with the exclusion code and reason for
each series, are in Supplementary Table 2.

Two levels of screening are kept apart by name throughout. The machine screen
is an automated keyword query over GEO metadata; it returned 1 424 candidate
series and decided nothing by itself. Checked by hand means that a candidate
was opened individually on its own GEO page, abstract and design read before
any computation, and a verdict with, where applicable, one of the
preregistered exclusion codes A1-A9 recorded at that moment; entry before
analysis was binding. Candidates reached the hand check from three sources,
recorded per series in Supplementary Table 2: the automated query itself, a
priority candidate list carried over from the predecessor project, and a
secondary candidate list assembled the same way.

<!-- PRIMARY_CITATIONS:START -->

Every reanalysed GEO series carries its accession in Supplementary Table 1
and in Supplementary Table 12. The twenty-four primary publications belong to
datasets rather than to claims, so they are cited once, collectively, at this
single point [Pansuriya et al. 2011; Berdasco et al. 2012; Massingham et al.
2014; Chiarelli et al. 2016; Kim et al. 2018; Barter et al. 2020; Cheung et
al. 2020; Zhytnik et al. 2020; Hernández et al. 2021; Reichenbach et al.
2021; Armes et al. 2022; Broeders et al. 2022; Chen et al. 2022; Novak et al.
2022; Tye et al. 2022; Gordon et al. 2023; La Manna et al. 2023; Schoenmaker
et al. 2023; Barter et al. 2024; Dinesh et al. 2024; Sun et al. 2024;
Vincent et al. 2024; Morales et al. 2025; Chu et al. 2026] and mapped
accession by accession in Supplementary Table 12. The 2 series without a
publication (GSE145235, GSE332758) are cited by accession.

<!-- PRIMARY_CITATIONS:END --> Raw data
(96 GB) are not redistributed; `code/00_setup.md` documents the download by
accession, and the processed matrices are included in the GitHub repository.

Gene annotation: GENCODE v46, hg38. Gene constraint: gnomAD LOEUF. Disease
panels: Genomics England PanelApp (panel "Skeletal dysplasia", 309 green genes;
a broader 1 471-gene version used as a sensitivity) and two panels derived from
the Nosology of Genetic Skeletal Disorders (core and broad)
[Unger et al. 2023]. Height GWAS genes [Yengo et al. 2022] and short-stature
gene sets as listed in Supplementary Table 6. Publication counts
per gene: NCBI `gene2pubmed`, tax_id 9606, retrieved 2026-08-22, mapped to
Ensembl through `gene2ensembl`; where the Entrez-to-Ensembl mapping is
ambiguous the maximum count was taken, which is the conservative direction
for a matching variable. Human fetal limb atlas: Nature 2023, Pcw 5.1–9.3,
136 311 cells × 26 522 genes, `X` log-normalised.

Second chromatin cohort (obtained 2026-08-23 for the lineage-contrast
test): ATAC-seq (GSE151311; AC and OB, two replicates each) and H3K27ac
(GSE151315; hMSC, AC and OB, two replicates each) from one study of
lineage-specific chromatin rearrangement in adipocyte and osteoblast commitment
[Hao et al. 2022], processed BigWig tracks, both in hg19. Window definitions are built from GENCODE v46 exactly as for the
first cohort and lifted point by point to hg19 with the UCSC hg38-to-hg19
chain (20 012 of 20 036 TSS uniquely liftable; gene bodies retained only when
both ends lift to the same contig and the body does not grow more than
threefold). Signal extraction, normalisation and every statistic follow
`reference_implementations/22_atac_window_calibration.R` and
`reference_implementations/_module.py` without modification
(`code/22_second_cohort_windows.R`, `code/23_lineage_contrast.py`).

### The module

The 173-gene programme was fixed before any analysis in this paper
(Supplementary Table 4), with direction `ri` (+1 in 129 genes, −1 in 44). It is
applied unchanged everywhere; it is never re-derived, re-fitted or re-selected
per level. Where a level cannot measure all 173 genes the measurable subset is
used and its size is reported per point.

Note that the enrichment of Figure 2D is computed on the 147 module genes
that appear in the background of 11 581, not on all 173. Both numbers are
correct for their context and must not be interchanged. The background is the
same for the narrow and the broad sets, so n_module = 147 holds for both passes;
only the gene sets differ. The sets themselves are deposited with source,
version and retrieval date in Supplementary Table 9, were fixed before the
run, and were not altered after any result was seen. Narrow and broad sets
are translated to Ensembl through the same bridge (Entrez →
`org.Hs.egENSEMBL`), so the two passes differ in set content and in nothing
else. The enrichment test has a single implementation
(`reference_implementations/_enrichment.R`), and both passes call it.

### Nulls

Every test has its own permutation null, and the null is always constructed on
the same data as the observation.

Gene-set tests draw background genes matched on expression decile and union
exon length, adding publication count as a third decile axis (6 × 6 × 6 = 216
cells) wherever it is a matching variable, and run 10 000 draws unless stated.
Module concordance per dataset permutes the gene-to-direction assignment within
the dataset's own measurable background. The in vivo trend permutes zone labels
within the specimen (20 000 draws, seed 20260823), so between-specimen
differences cannot produce a trend. ATAC carries two nulls: a background draw,
and H1, stratified by baseline accessibility in the undifferentiated state,
which controls for module genes simply being those closed at day 0 and able
only to open; H1 is the harder null and the one quoted in the main text, and a
second hardening, H2, additionally stratifies by window width for the
gene-body window. Seeds are fixed and recorded in each output file.

### Calibration, per level

**The rule is the same at every level and was fixed before the data were
seen** (preregistration M-D, §6): the textbook markers of the axis being measured are
tested against the same null as the module, and the level passes at z ≥ 2.
A level that fails its calibration can produce neither a positive nor a negative
finding; its result is *not measurable*. This is applied without exception,
including where it costs us the most attractive panel in the paper (the ATAC
lineage contrast; see Discussion).

The marker sets are given in `reference_implementations/_marker.py` and in
Supplementary Table 3.

The five sets are canonical throughout [Pittenger et al. 1999]. The osteogenic
set pairs transcriptional regulators (*RUNX2*, *SP7*, *DLX5*, *MSX2*,
*SATB2*, *ATF4*) with early matrix genes (*ALPL*, *BGLAP*, *IBSP*, *SPP1*,
*COL1A1*, *COL1A2*, *POSTN*) and mineralisation markers (*MEPE*, *PHEX*,
*DMP1*, *SOST*, *PTH1R*); the chondrogenic set covers the SOX trio and
cartilage matrix (*SOX9*, *SOX5*, *SOX6*, *ACAN*, *COL2A1*, *COL9A1*,
*COL11A1*, *COMP*, *HAPLN1*, *MATN3*, *PRG4*, *WWP2*) [Goldring et al. 2006];
the adipogenic set pairs core regulators (*PPARG*, *CEBPA*, *CEBPB*, *CEBPD*,
*SREBF1*, *NR1H3*) with lipid-handling genes (*FABP4*, *ADIPOQ*, *LEP*,
*LPL*, *PLIN1*, *PLIN4*, *CIDEC*, *CD36*, *GPD1*, *AQP7*, *LIPE*, *PNPLA2*)
[Gregoire et al. 1998]; the myogenic lineage control lists determinants and
structural genes (*MYOD1*, *MYOG*, *MYF5*, *MYF6*, *PAX7*, *DES*, *MYH3*,
*MYH8*, *TNNT1*, *TNNT2*, *ACTA1*, *CKM*, *MYL4*, *TTN*, *DMD*, *CAV3*); and
the ten-gene undifferentiated set combines mesenchymal stromal surface
markers (*THY1*, *ENG*, *NT5E*, *NGFR*) with stemness and proliferation genes
(*KITLG*, *LIF*, *CXCL12*, *MKI67*, *CCNB1*, *TOP2A*) [Dominici et al. 2006].

A sensitivity analysis omitting the
markers not reachable in vitro (*SOST*, *DMP1*, *PHEX*, *MEPE*, *PTH1R*) is in Fig. S5C
and gives the same head count.

The gene space of the calibration, and why it is not the gene space of the
module. The module is defined on the project's internal gene map: genes with
an evaluable differentiation response in at least 15 of the 18 datasets. That
threshold is necessary there: a gene measurable in three datasets cannot be
called convergent across eighteen. It is wrong for a per-experiment positive
control, and not neutrally so. Terminal differentiation markers are off in
the undifferentiated arm and become measurable only where differentiation actually runs,
so a "measurable in ≥ 15 of 18" filter removes them preferentially: it drops
*SP7*, *BGLAP*, *IBSP*, *DMP1*, *SOST* and *MEPE* from the osteogenic set, *ACAN*, *COL2A1*,
*COL9A1*, *COMP*, *HAPLN1* and *PRG4* from the chondrogenic set, and most of the mature
adipocyte set, while leaving the undifferentiated marker set complete at
10 of 10. The filter therefore
acts asymmetrically against the half of the calibration that measures
*arrival*, and it does so in the datasets where arrival occurred. The datasets
themselves measure these genes: one measures 18 of 18 osteogenic markers where
the filter passed 12, and four measure 12 of 12 chondrogenic markers where the
filter passed 5.

The calibration therefore runs on unfiltered marker sets and an unfiltered
background, which is also what exchangeability of the permutation null
requires: marker sets and background must come from the same population. A
decomposition of the two changes separately
(`code/29_calibration_gene_space.py`) shows that the marker set alone moves the
statistic (median |Δz| 0.378 across the 18 datasets, maximum 2.23) and the
background does not (median 0.017, maximum 0.22). On the unfiltered marker sets
the head count is 2 of 18; the per-dataset values for both gene spaces are in
Supplementary Table 3.

**What the restored markers show is not bookkeeping.** In the two datasets that
the correction reclassifies, the early osteogenic regulators move and the
mineralisation markers do not. In *RB1*-mut isogenic, *RUNX2* (+1.47), *ALPL*
(+1.88), *COL1A1*/*COL1A2* (+1.35/+1.51), *POSTN* (+2.02) and *ATF4* (+2.01) all rise
while *BGLAP*, *SOST* and *MEPE* sit at exactly zero and *DMP1* falls (−2.22); in
*TP53* LFS, *RUNX2* (+2.67) and *DLX5* (+1.97) rise while *SP7* (−0.74), *SOST* (−0.83),
*DMP1* (−0.73) and *MEPE* (−2.00) fall. These cultures engage the early
osteogenic programme and stop before terminal differentiation, which is the
same statement the decomposition of Figure 2F makes, arrived at independently,
and it is visible only when the terminal markers are left in.

### Detection limits

For every test we report MDE80: the smallest true effect the analysis would
detect in 80 % of repetitions, computed from the same permutation null as the
test (the 95th percentile of the null plus 0.8416 null standard deviations).
Supplementary Table 7 lists every statistic in the paper with its own limit;
Figure 6 is the level-by-level summary. No number in this paper appears
without one, and where the limit is unreachable the result is reported as
*not measurable* rather than as an absence of effect.

### Self-test against the known null rate

The metric was tested where no effect exists by construction. In the donor
set-up, |z| > 2 occurs in 3.0–5.0 % of runs across four statistics at zero
injected effect, against a nominal 5 %. Across 520 neutral contrasts of three
kinds it occurs in 6.3 % (6.9 / 8.3 / 3.9 % by kind), i.e. slightly
anticonservative; this is stated wherever a z = 2 threshold is used. Sensitivity:
the metric recovers an additive offset of 0.35 z at 60 genes per side, 0.50 at
30 and 0.75 at 15 (Fig. S1B, S1C, S5B).

### Rules that constrain the analysis, and the counter-example that earned each

Five rules constrain the analysis, each earned against a counter-example.

No adjustment is made for any covariate of the baseline on the z scale. The
per-gene z standardisation forces cor(baseline, dWT) = −0.566 (median, negative
at all eleven points), and the machinery genes are the highly expressed ones
(SMD +0.527); adjusting therefore adjusts away part of the target, destroying
the positive control (+9.10 → −2.43) while making 10 of 11 points null-intact.

No two convergence counts are compared at unequal signal-to-noise.
   Compressing the shared component of `dWT` by 35 % costs 91 % of convergent
   genes (120 → 11), and a lesion-free pseudo-contrast built from replicates of
   the control arm produces 4.3 convergent genes against the real 6.0, with a
   97.5th percentile of the noise floor at 15.0.
The donor is the unit, without exception. Clones and repeated cultures of
one line are not donors. This is the rule that costs us the most cells
   (GSE221128 becomes one donor rather than three "experiments"; the *LAMA5*
   series becomes one donor rather than six clones), and it is the rule the
   literature supports: in population-scale iPSC panels the dominant variance
   component is the donor genetic background, above culture condition,
   passage and sex, and lines from one individual resemble each other more than
   lines from different individuals [Kilpinen et al. 2017]. The same
literature also documents how strongly protocol and batch drive variability in
iPSC-derived models, and recommends reporting it explicitly
[Volpato and Webber 2020]; the analogous problem in another lineage (models
that do not reach the mature state they are read as) is well described for
iPSC-derived cardiomyocytes [Wu et al. 2021].
`dWT` is a positive control and a tool, never a result: a quantity that
licenses an analysis cannot also be the result of it. And no arithmetic mean
is taken across cohorts or studies.

### In vivo

Unit of analysis: the specimen (`adj_sample` prefix `Pcw<stage>[_s<i>]`), not
the cell: 16 specimens across 9 developmental stages. A (zone, sample) pair
enters only with ≥ 5 cells, giving 66 evaluable pairs. The contrast is
module-up genes (ri = +1, 129) against module-down genes (ri = −1, 44) against a
background-drawn null of 4 000 random genes expressed in the atlas
(`reference_implementations/_module.py`). The trend test is Spearman ρ between
zone rank
(1…6) and the per-sample contrast, against the within-specimen permutation null
described above. The mandatory positive control is the same trend logic applied
to the chondrogenic marker set against the undifferentiated one.

Two versions of this trend test exist and both are reported. Under a free
permutation of samples the trend is ρ 0.456, z +3.65; under the
specimen-stratified permutation used throughout this paper
(`code/20_in_vivo_donor_test.py`) the same ρ 0.456 gives z +4.80 at a limit of
ρ 0.274. The stratified version is the stricter test, and z rises rather than
falls because permuting within a specimen raises the null mean but shrinks its
spread more. The stratified value is the one quoted in the Results.

### Multiple testing

Bonferroni within each family, with the family stated: 6 tests for the GO
composition (α = 0.0083); 8 for the complementarity tests (α = 0.00625); 40
computable tests for the mechanism classes (α = 0.00125, from 49 attempted);
21 for the matched gene-set tests (α = 0.00238).

### Preregistrations

All preregistration and protocol documents are in `preregistrations/`, dated,
reproduced unchanged, and listed with type and title in Supplementary Table 8. They include the
registrations of hypotheses that subsequently fell, and the protocols recording
how they fell. Preregistered claims are made only for analyses covered by them:
the screens, the human-genetics anchor and its positive controls, the
donor-level programme statistic, and the external triangulation follow-up.
Everything in Figures 2, 3, 5 and in panels C, D, E, F of Figure 4 is
exploratory and is labelled as such in every legend.

### Scale critique

The full version, including why several intuitive analyses are absent, is in
Fig. S1 and in the rules above. In short: the per-gene z scale on which the
module is defined cannot support a question about absolute expression level, it
cannot support adjustment for baseline expression, and a convergence count on it
is a function of the signal-to-noise ratio rather than of biology. Where a
question needs the absolute scale, the absolute expression matrices are loaded
and the question is asked there, which is what made the publication-matched
analyses of Figure 4 possible.

### Software

R 4.4.3 with ggplot2, patchwork, ragg, systemfonts, DESeq2 and matrixStats;
Python 3.12 with numpy, pandas and scipy. Exact versions are written by each
script into `results/` as a session-information file alongside its log.
Figures are drawn by `code/60_figures_main.R` and
`code/61_figures_supplement.R` against `figure_style/publication_style.R`.

---

## References

1. Unger S. *et al.* Nosology of genetic skeletal disorders: 2023 revision.
   *Am J Med Genet A* **191**, 1164–1209 (2023). doi:10.1002/ajmg.a.63132
2. Yengo L. *et al.* A saturated map of common genetic variants associated with
   human height. *Nature* **610**, 704–712 (2022).
   doi:10.1038/s41586-022-05275-y

3. Duggan S. Vosoritide: first approval. *Drugs* **81**, 2057–2062 (2021).
   PMID 34694597.
4. Savarirayan R. *et al.* Vosoritide treatment for children with
   hypochondroplasia: a phase 2 trial. *eClinicalMedicine* (2024).
5. Savarirayan R. *et al.* International consensus guidelines on the
   implementation and monitoring of vosoritide therapy. *Nat Rev Endocrinol*
   (2024). doi:10.1038/s41574-024-01074-9
6. Cameron T. L. *et al.* XBP1-independent UPR pathways suppress C/EBP-β
   mediated chondrocyte differentiation in ER-stress related skeletal disease.
   *PLoS Genet* **11**, e1005505 (2015).
7. Cameron T. L. *et al.* Transcriptional profiling of chondrodysplasia
   growth plate cartilage reveals adaptive ER-stress networks that allow
   survival but disrupt hypertrophy. *PLoS One* **6**, e24600 (2011).
   doi:10.1371/journal.pone.0024600
8. Farhan H. *et al.* The endoplasmic reticulum proteostasis network and
   bone disease. *Trends Mol Med* **32**, 119–132 (2026).
   doi:10.1016/j.molmed.2025.06.005

9. Tang V. T. & Ginsburg D. Cargo selection in endoplasmic
   reticulum-to-Golgi transport and relevant diseases. *J Clin Invest*
   **133**, e163838 (2023). doi:10.1172/JCI163838
10. Sarmah S. *et al.* Sec24D-dependent transport of extracellular matrix
    proteins is required for zebrafish skeletal morphogenesis. *PLoS One*
    **5**, e10367 (2010). doi:10.1371/journal.pone.0010367
11. Claeys L. *et al.* Collagen transport and related pathways in
    osteogenesis imperfecta. *Hum Genet* **140**, 1121–1141 (2021).
    doi:10.1007/s00439-021-02302-2
12. Lipiński P. *et al.* Skeletal and bone mineral density features,
    genetic profile in congenital disorders of glycosylation: a review.
    *Diagnostics* **11**, 1438 (2021). doi:10.3390/diagnostics11081438

13. Kilpinen H. *et al.* Common genetic variation drives molecular
    heterogeneity in human iPSCs. *Nature* **546**, 370–375 (2017).
14. Volpato V. & Webber C. Addressing variability in iPSC-derived models of
    human disease: guidelines to promote reproducibility. *Dis Model Mech*
    **13**, dmm042317 (2020).
15. Wu P. *et al.* Maturation strategies and limitations of induced
    pluripotent stem cell-derived cardiomyocytes. *Biosci Rep* **41**,
    BSR20200833 (2021). doi:10.1042/BSR20200833

16. Hao R.-H. *et al.* Lineage-specific rearrangement of chromatin loops and
    epigenomic features during adipocytes and osteoblasts commitment.
    *Cell Death Differ* **29**, 2503–2518 (2022).
    doi:10.1038/s41418-022-01035-7

17. Liu T. M. *et al.* Identification of common pathways mediating
    differentiation of bone marrow- and adipose tissue-derived human
    mesenchymal stem cells into three mesenchymal lineages.
    *Stem Cells* **25**, 750–760 (2007). doi:10.1634/stemcells.2006-0394
18. He J. *et al.* Directing the osteoblastic and chondrocytic
    differentiations of mesenchymal stem cells: matrix vs. induction media.
    *Regen Biomater* **4**, 269–279 (2017). doi:10.1093/rb/rbx008

19. Pittenger MF *et al.* Multilineage potential of adult human mesenchymal
    stem cells. *Science* **284**, 143-147 (1999). doi:10.1126/science.284.5411.143

20. Dominici M *et al.* Minimal criteria for defining multipotent mesenchymal
    stromal cells. The International Society for Cellular Therapy position
    statement. *Cytotherapy* **8**, 315-317 (2006). doi:10.1080/14653240600855905

21. Goldring MB, Tsuchimochi K & Ijiri K. The control of chondrogenesis.
    *J Cell Biochem* **97**, 33-44 (2006). doi:10.1002/jcb.20652

22. Gregoire FM, Smas CM & Sul HS. Understanding adipocyte differentiation.
    *Physiol Rev* **78**, 783-809 (1998). doi:10.1152/physrev.1998.78.3.783

<!-- PRIMARY_SOURCES:START -->

1. Pansuriya TC *et al.* Genome-wide analysis of Ollier disease: Is it all in the genes? *Orphanet J Rare Dis* **6**, 2 (2011). doi:10.1186/1750-1172-6-2

2. Berdasco M *et al.* DNA methylation plasticity of human adipose-derived stem cells in lineage commitment. *Am J Pathol* **181**, 2079-93 (2012). doi:10.1016/j.ajpath.2012.08.016

3. Massingham LJ *et al.* Amniotic fluid RNA gene expression profiling provides insights into the phenotype of Turner syndrome. *Hum Genet* **133**, 1075-82 (2014). doi:10.1007/s00439-014-1448-y

4. Chiarelli N *et al.* Transcriptome-Wide Expression Profiling in Skin Fibroblasts of Patients with Joint Hypermobility Syndrome/Ehlers-Danlos Syndrome Hypermobility Type. *PLoS One* **11**, e0161347 (2016). doi:10.1371/journal.pone.0161347

5. Kim H *et al.* Oncogenic role of SFRP2 in p53-mutant osteosarcoma development via autocrine and paracrine mechanism. *Proc Natl Acad Sci U S A* **115**, E11128-E11137 (2018). doi:10.1073/pnas.1814044115

6. Barter MJ *et al.* DNA hypomethylation during MSC chondrogenesis occurs predominantly at enhancer regions. *Sci Rep* **10**, 1169 (2020). doi:10.1038/s41598-020-58093-5

7. Cheung K *et al.* Histone ChIP-Seq identifies differential enhancer usage during chondrogenesis as critical for defining cell-type specificity. *FASEB J* **34**, 5317-5331 (2020). doi:10.1096/fj.201902061RR

8. Zhytnik L *et al.* RNA sequencing analysis reveals increased expression of interferon signaling genes and dysregulation of bone metabolism affecting pathways in the whole blood of patients with osteogenesis imperfecta. *BMC Med Genomics* **13**, 177 (2020). doi:10.1186/s12920-020-00825-7

9. Hernández R *et al.* Impact of the Epigenetically Regulated Hoxa-5 Gene in Neural Differentiation from Human Adipose-Derived Stem Cells. *Biology (Basel)* **10**, 802 (2021). doi:10.3390/biology10080802

10. Reichenbach M *et al.* Differential impact of fluid shear stress and YAP/TAZ on BMP/TGF-beta induced osteogenic target genes. *Adv Biol (Weinh)* **5**, 2000051 (2021). doi:10.1002/adbi.202000051

11. Armes H *et al.* Germline ERCC excision repair 6 like 2 (*ERCC6L2*) mutations lead to impaired erythropoiesis and reshaping of the bone marrow microenvironment. *Br J Haematol* **199**, 754-764 (2022). doi:10.1111/bjh.18466

12. Broeders M *et al.* Modeling cartilage pathology in mucopolysaccharidosis VI using iPSCs reveals early dysregulation of chondrogenic and metabolic gene expression. *Front Bioeng Biotechnol* **10**, 949063 (2022). doi:10.3389/fbioe.2022.949063

13. Chen P *et al.* Phenotypic Spectrum and Molecular Basis in a Chinese Cohort of Osteogenesis Imperfecta With Mutations in Type I Collagen. *Front Genet* **13**, 816078 (2022). doi:10.3389/fgene.2022.816078

14. Novak R *et al.* RNF4~RGMb~BMP6 axis required for osteogenic differentiation and cancer cell survival. *Cell Death Dis* **13**, 820 (2022). doi:10.1038/s41419-022-05262-1

15. Tye CE *et al.* LncMIR181A1HG is a novel chromatin-bound epigenetic suppressor of early stage osteogenic lineage commitment. *Sci Rep* **12**, 7770 (2022). doi:10.1038/s41598-022-11814-4

16. Gordon JAR *et al.* LINC01638 sustains human mesenchymal stem cell self-renewal and competency for osteogenic cell fate. *Sci Rep* **13**, 20314 (2023). doi:10.1038/s41598-023-46202-z

17. La Manna F *et al.* Molecular profiling of osteoprogenitor cells reveals FOS as a master regulator of bone non-union. *Gene* **874**, 147481 (2023). doi:10.1016/j.gene.2023.147481

18. Schoenmaker T *et al.* Transcriptomic Differences Underlying the Activin-A Induced Large Osteoclast Formation in Both Healthy Control and Fibrodysplasia Ossificans Progressiva Osteoclasts. *Int J Mol Sci* **24**, 6822 (2023). doi:10.3390/ijms24076822

19. Barter MJ *et al.* SERPINA3 is a marker of cartilage differentiation and is essential for the expression of extracellular matrix genes during early chondrogenesis. *Matrix Biol* **133**, 33-42 (2024). doi:10.1016/j.matbio.2024.07.004

20. Dinesh NEH *et al.* Mutations in fibronectin dysregulate chondrogenesis in skeletal dysplasia. *Cell Mol Life Sci* **81**, 419 (2024). doi:10.1007/s00018-024-05444-4

21. Sun L *et al.* Oxidative phosphorylation is a pivotal therapeutic target of fibrodysplasia ossificans progressiva. *Life Sci Alliance* **7**, e202302219 (2024). doi:10.26508/lsa.202302219

22. Vincent A *et al.* Monoallelic loss of *RB1* enhances osteogenic differentiation and delays DNA repair without inducing tumorigenicity. *Differentiation* **140**, 100815 (2024). doi:10.1016/j.diff.2024.100815

23. Morales AA *et al.* Dysregulation of cell migration by matrix metalloproteinases in geleophysic dysplasia. *Sci Rep* **15**, 19970 (2025). doi:10.1038/s41598-025-04666-1

24. Chu TLN *et al.* A transcriptional atlas of the pubertal human growth plate reveals two populations of stem cells and direct effect of growth hormone. *Sci Transl Med* **18**, eadw3590 (2026). doi:10.1126/scitranslmed.adw3590

25. Time series of retinoblastoma patient-derived mesenchymal stem cells differentiated to the osteogenic lineage. Gene Expression Omnibus, GSE145235.

26. Alterations in chromatin accessibility during osteoblast and adipocyte differentiation in human mesenchymal stem cells. Gene Expression Omnibus, GSE332758.


<!-- PRIMARY_SOURCES:END -->

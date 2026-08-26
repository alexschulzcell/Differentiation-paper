# Main figure legends

Conventions used in every legend. MDE80 is the smallest true effect the
analysis would detect in 80 % of repetitions, estimated from the same
permutation null as the test itself. Calibration means the built-in positive
control of a level: the textbook lineage markers of the axis being measured,
tested against the same null as the module, passing at z ≥ 2. The unit of
analysis is the donor throughout; clones and replicate cultures are not donors.
Each panel is labelled confirmatory (preregistered), exploratory, or
preregistered follow-up.

---

## Figure 1 · What the datasets are, and the calibration that almost none of them passes

**(A)** The two layers this paper separates: a lineage-independent matrix
programme of 173 genes, fixed before any analysis reported here, and the
skeletal dysplasia genes of the PanelApp "Skeletal dysplasia" panel (309 green
genes). Schematic; no data.

**(B)** Screen of perturbation datasets. Of 89 candidate series, 36 were
excluded for having no undifferentiated arm and 35 for other reasons, leaving
the 18 datasets used throughout. Confirmatory; the exclusion codes were fixed
before the screen.

**(C)** Screen for series carrying a diagnosis as the lesion axis, log scale,
by two independent strategies. By entity: 1 424 GEO series machine-screened,
127 checked by hand, 50 carry a diagnosis axis, of which 46 (92 %) have no
undifferentiated arm. By design: 22 series, of which 2 complete 2 × 2 designs,
both already among the 18. Confirmatory.

**(D)** The calibration, per dataset: textbook markers of each dataset's own
axis against its own permutation null (10 000 draws; background genes matched
on expression decile and union exon length). Two of eighteen pass: *LAMA5*-KO
chondrogenic (z +2.67) and *ERCC6L2*-KD (z +2.13). *SERPINA3*-KD chondrogenic
is not calibratable at all (a single marker of its axis survives; open circle
at zero), so counted strictly 2 of 17 calibratable datasets pass. Two more come
close without reaching the threshold (*ACVR1*/FOP z +1.90, *RB1*-mut isogenic
z +1.93) and count as failures. No dataset reaches its own MDE80 here (that
limit sits at z ≈ 2.8), so the two passes rest on the preregistered z ≥ 2 rule
alone. Dashed line: the calibration threshold z = 2. The marker sets are
deliberately unfiltered; reason and effect in STAR Methods.

**(E)** The same calibration resolved by donor cell: 7 of 14 cells from 6
studies pass. Diamonds: genuine patient lesions; circles: engineered
perturbations; exactly one patient-lesion cell passes, which is why no
donor-resolved lesion-response number is reported anywhere in this paper. The
two head counts differ in resolution, not material (all six donor-level study
units are also among the 18), so the 32 calibration attempts are not 32
independent experiments. Panels D and E apply a rule fixed before the main
analysis and are not themselves results.


---

## Figure 2 · The programme is lineage-independent and decoupled from lineage commitment

Panels A to E are exploratory: they ran after the preregistered questions had
been closed, and no preregistered claims are made for them. Panel F is a
preregistered follow-up.

**(A)** Cross-arm concordance of the 173 module genes between the osteogenic
arm (12 datasets) and the chondrogenic arm (6 datasets): Spearman ρ +0.622
against a null of 0.091 ± 0.076 from permutation of the gene-to-arm assignment
over 35 572 background genes, z +7.03, P = 0.0004. Grey area: null density.

**(B)** Module concordance per dataset, 18 points, each against its own
detection limit (vertical bar); all 18 lie above their own limit, z +5.25 to
+13.10. Colour: each dataset's own calibration status (green passed, 2 of 18;
red failed; grey not calibratable). The n per point is module genes measurable
in that dataset, 29 to 173.

**(C)** Pooled by calibration status. Datasets that fail their calibration
(n = 16) give concordance 1.000, z +13.13; those that pass (n = 2) give 0.994,
z +12.79, both far above their MDE80 of 0.617 and 0.612. Grey bracket: null
mean ± SD; black bar: MDE80.

**(D)** What the programme consists of, against external gene sets (Fisher,
two-sided, 95 % CI; 6 primary tests, Bonferroni α = 0.0083). Filled circles:
broad, independently curated sets (matrisome via MSigDB 2026.1.Hs for the two
matrix categories; Reactome for TGFβ/BMP and hypoxia and, added to GO, for
cell-cycle exit); open circles: the narrow single GO terms of the first pass;
secretory stays on GO, having no broader set. Sets and decision rule were
fixed, with source, version and retrieval date, before the run. Broad sets, in
panel order: matrix remodelling OR 3.62 [1.14–8.91] (narrow 7.92 [1.53–25.9]);
matrix components OR 3.38 [1.20–7.73] (narrow 4.03 [1.26–9.95]); secretory
machinery OR 0.32 [0.10–0.78], depleted (the same set in both passes).
Cell-cycle exit keeps its direction but its interval now includes 1
(OR 1.80 [0.93–3.20]; narrow 2.11 [1.09–3.77]), so by the pre-set rule its
number lives in Supplementary Table 9. TGFβ/BMP and hypoxia/stress are null in
both passes. Three qualifications: k stays small, 5 to 13 genes per category,
because only 147 of the 173 module genes lie in the background of 11 581 at
all, so the panel supports a direction, not a magnitude; no category reaches
Bonferroni significance in either pass (smallest corrected P 0.042), the
pre-set criterion being a 95 % CI excluding 1; and Reactome ECM organisation,
degradation, collagen formation and proteoglycans, an independent second source
for the two matrix categories, reproduces the direction without reaching the
criterion (OR 1.55 [0.18–5.84], OR 1.45 [0.39–3.84]). k is printed for both
sets.

**(E)** The same question in chromatin. ATAC-seq (GSE332758) carries
osteogenic and adipogenic axes in the same cells, in four window definitions
(promoter; TSS ± 10 kb; TSS ± 50 kb; gene body). Module concordance per axis
against its own MDE80 (vertical bar); each row's calibration z printed in the
field, points faded where the axis fails. The null shown is H1, stratified by
baseline accessibility in the undifferentiated state, the harder of the two
nulls: it controls for module genes simply being closed at day 0 and able only
to open. The adipogenic axis passes its calibration in all four windows
(z +3.73 to +5.01) and the module lies above its limit in 3 of 4 windows under
H1 (background null: 2 of 4), z +2.38 to +4.76: a calibrated positive result on
a non-skeletal lineage. The osteogenic axis fails its calibration in all four
windows (z −0.54 to +1.66) while the module runs in all four (z +3.50 to
+4.51); reported as an observation about decoupling, not as a calibrated module
result. The lineage contrast, osteogenic minus adipogenic, fails its own
calibration in every window (z +0.93 to +2.26), so its 0 of 4 means not
measurable rather than nothing there. Background null: Figure S6.

**(F)** The decoupling taken apart across all 18 datasets; preregistered
follow-up. All 18 carry an undifferentiated arm, so the calibration decomposes
into whether the culture left the undifferentiated state, whether it reached
its lineage, and whether the module ran; hypothesis and decision rule were
written down before computing. One point per dataset, coloured by its own
calibration status (green = passed, 2; red failed, 15;
grey = not calibratable, 1, drawn without verdict); dashed lines mark z = ± 2, the preregistered
thresholds. Left field: the three marker-set contrasts against the usual
background null. Right field: module concordance z on its own axis, since it
runs to +13.1 while the marker contrasts lie between −4 and +4. The
undifferentiated markers fall (z ≤ −2) in 10 of 18 datasets; the other 8 are
not decomposable and carry no verdict. Own-lineage markers rise above z = +2
in 2 of 18; the module lies above its own limit in 18 of 18. Among the 10
decomposable datasets the rule gives 8 confirmations (Wilson 0.49–0.94),
2 instances of the other case of decoupling, and 0 refutations. Two things are
stated rather than implied: the middle criterion is the same statistic at the
same threshold as in panel D, so its agreement with calibration status is
definitional and not a finding; and the 18 are not independent experiments (six
study units also appear at donor resolution, Figure 1E), so no pooled statistic
and no mean over the 18 z values is computed, the panel reporting a count with
Wilson interval, every point listed individually (Supplementary Table 10).


---

## Figure 3 · Not a culture artefact: the human fetal growth plate

Exploratory. Human fetal limb single-cell atlas, post-conception weeks 5.1 to
9.3, 136 311 cells. The unit of analysis is the specimen, not the cell: 66
evaluable (zone, sample) pairs from 16 specimens across 9 developmental stages,
entering with at least 5 cells.

**(A)** The chondrogenic axis of the atlas (mesenchymal condensation,
chondroprogenitor, resting, proliferative, prehypertrophic, hypertrophic),
cells per zone on a log scale.

**(B)** Positive control: textbook chondrogenic minus undifferentiated markers,
per sample along the same axis. 64 of 66 samples lie above their own detection
limit (65 of 66 at z > 2), median z +10.92. Grey dashes: per-sample MDE80.

**(C)** Module along the same axis, one line per specimen. Spearman ρ between
zone rank and module contrast is +0.456, z +4.80, P = 5 × 10⁻⁵, MDE80 ρ 0.274;
the null permutes zone labels within the specimen (20 000 draws), so no
between-specimen difference in means can produce the trend. Maximum:
prehypertrophic zone. Two sensitivity analyses under rules fixed before the run
qualify the trend. Removing the 10 % of module genes with the largest absolute
zone-to-zone difference leaves ρ 0.144, z +2.74, above the reduced module's own
limit of ρ 0.125; it survives narrowly (ρ falls threefold, the 17 largest genes
carry 57 % of the absolute difference, and at 20 % or 30 % removal the trend
falls below its limit in runs declared descriptive in advance). Dropping every
hypertrophic point leaves the trend at ρ +0.430, z +4.65, above its recomputed
limit of ρ 0.250, with the positive control still calibrated (ρ +0.880 against
0.615). Only 76 of 164 module genes are individually concordant (46.3 %, 95 %
CI 38.7–54.0 %), an interval covering chance, so no claim is made about which
parts of the programme are confirmed in vivo. Supplementary Tables 11, 11b and
13.

**(D)** Evaluable samples per zone. The hypertrophic zone contributes a single
sample of 8 cells (14 in the whole atlas); the trend therefore carries up to
the prehypertrophic zone and no further, and nothing rests on the hypertrophic
end point.


---

## Figure 4 · The other layer: what defines the disease genes

Panels A and B are confirmatory: the human-genetics anchor and its positive
controls were preregistered. Panels C, D, E and F are exploratory. The figure
answers one question: what defines the skeletal dysplasia genes, if not
differentiation dynamics?

**(A)** Positive controls of the set-up. Top: lineage markers inside the disease
panels, OR 17.0 to 51.6 after matching on expression decile and exon length,
z up to +13.9. Bottom: the secretion anchor, distal against biosynthetic
secretion, OR 3.70 to 5.75, P down to 2.9 × 10⁻¹¹ (Fisher). Vertical bar: the
test's own MDE80.

**(B)** Axis 1, localisation. The same matched enrichment applied to the
programme and to three disease gene sets, in the distal and the biosynthetic
secretion compartment; 8 tests, Bonferroni α = 0.00625. Disease genes are
enriched in distal secretion (PanelApp 309: OR 2.84, z +6.30, P = 1 × 10⁻⁴,
limit OR 1.79) and depleted in biosynthetic secretion (OR 0.54, z −3.18). The
programme is in neither: distal OR 0.81, z −0.41, not significant at a limit of
OR 2.34. Vertical bar: own MDE80.

**(C)** Axis 2, gene dosage. PanelApp 309 split by mode of inheritance, purely
monoallelic against purely biallelic: LOEUF 0.283 (n = 120) against 0.826
(n = 247), P = 6 × 10⁻²⁴ (Mann–Whitney). On the very same genes, absolute dWT
0.406 against 0.361 (n = 120 and 249), P = 0.86; the two n differ because LOEUF
is missing for two biallelic genes. The dosage result holds within every
publication tertile (P 0.017, 2.2 × 10⁻⁵, 1.7 × 10⁻¹¹), so it is not an
artefact of study intensity.

**(D)** The constraint contrast. LOEUF after matching on expression decile,
exon length and publication count per gene. Common height variation sits in
constrained genes (height GWAS z −4.67, limit 0.024 LOEUF units); monogenic
dysplasia genes do not (PanelApp 309 z +1.83, Nosology core +1.30, both not
significant at limits of 0.095 and 0.139). The negative control, the cell cycle,
collapses from z −6.99 to −2.54 once publications are matched, which is what
the matching is for. Limits are printed at the right.

**(E)** The axis that does not separate them, which is the point of the figure.
Absolute dWT (differentiation dynamics) at equal expression, exon length and
publication count. In the same run the programme sits at z +18.10 and the
cell-cycle negative control at z +3.72, while the disease panels sit on the
null: PanelApp 309 z −0.71 at a limit of 0.073 dWT units, height GWAS z −0.81 at
a limit of 0.016. An effect of that size would have been found.

**(F)** The closing negative result. Seven mechanism classes against seven
panels give 49 tests, of which 40 are computable (nine have too few genes); at
Bonferroni α = 0.05/40 = 0.00125 none is significant. Grey bars: each test's own
MDE80. Four tests lie above their own limit (red), all four in the
glycosylation and linker class (OR 6.45 to 8.59, nominal P 0.042 to 0.076,
effectively two independent values because two panel pairs are identical): the
chance expectation at 40 tests on a nominal 5 % level, reported rather than
omitted. Pooled over all classes, the sharpest axis (height GWAS, 5 649 genes)
gives OR 1.00 at a limit of OR 1.59, and on none of the seven panels does the
pooled odds ratio exceed its own limit.


---

## Figure 5 · Both layers meet at the prehypertrophic transition

Exploratory. Same atlas, same unit of analysis and same null as Figure 3.

**(A)** Contrast of the disease genes (PanelApp 309, 355 measurable) against the
module, per zone. Grey points: individual samples (n = 66); blue point and
range: median and interquartile range; median z per zone printed above. Values
run 0.075, 0.112, 0.123, 0.137, 0.168, 0.132 with median z +2.66, +3.79, +4.09,
+4.42, +4.98, +4.17. The last zone rests on a single sample (Figure 3D).

**(B)** Both curves on the same axis, each as the median and interquartile range
of the per-sample contrast against a background-drawn null of 4 000 random genes
expressed in the atlas. Both have their maximum in the prehypertrophic zone.
The disease genes lie above the programme in every zone up to the prehypertrophic
transition. This is a statement about spatial co-localisation, not about
enrichment (Figure 4E is the test that fell) and not about direction of
regulation.

**(C)** The schematic that summarises the paper: the programme is what the cell
does, the disease genes are the machinery that executes it. Dysplasia is a
failure of the machinery, not of the programme, which is why there is no
transcriptional downstream convergence to be found. Schematic; no data.


---

## Figure 6 · Levels, detection limits, and what each level carries

Confirmatory as bookkeeping: this figure asserts limits, not effects, and
deliberately includes the levels that carry nothing.

One row per level: its estimate (a point, or a line with end points where
several windows or tests feed one level) and its own detection limit (grey bar)
in that level's own unit; for the donor level the null mean ± SD is drawn as a
thinner grey bar. Dot colour codes the verdict: carries (green); good negative
result (blue: a negative finding with a positive control in the same set-up);
not measurable (red: the level's own calibration fails, licensing no statement
of either sign); observation only (orange); carries nothing (grey). Three
stacked fields hold statistics that cannot share an axis: **(A)** unit-free
concordance, correlation and share; **(B)** odds ratio on a log scale, carrying
the pooled height-GWAS enrichment against its limit of 1.59; **(C)** marker and
module contrast z values. Rows whose detection limit lives in another unit than
their headline z (dynamics axis, 0.073 dWT units; gene constraint, 0.024 to
0.139 LOEUF units) carry no grey bar; their limits are quoted in Supplementary
Table 14 instead.

Red versus grey carries the point of the figure: a level whose calibration
fails cannot deliver a negative result, and calling its zero no effect would be
an error; the ATAC lineage contrast is the clearest case, below its own
per-window limits with its axis failing calibration in all four windows. Three
levels are listed although they carry nothing (both promoter methylomes; the
undifferentiated state as a predictor of later differentiation) because
checked-and-failed levels make the limits interpretable. The module row of the
18 datasets is drawn as its range, z +5.25 to +13.10, each dataset having its
own limit (Supplementary Table 7); the level-by-level text table is
Supplementary Table 14.


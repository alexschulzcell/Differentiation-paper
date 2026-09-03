# What skeletal differentiation models measure: a shared matrix programme, separate from skeletal dysplasia genes

*Companion repository for Schulz & Thiel (manuscript under submission). The manuscript itself is not part of this repository.*

A lineage-independent matrix programme runs in every published mesenchymal
differentiation model we could reach — including where the model fails its own
positive control — while the skeletal dysplasia genes are defined by
localisation and gene dosage and lie orthogonal to the differentiation-dynamics
axis the field searches on. Every level carries its own internal positive
control ("calibration") and its own measured detection limit (MDE80), and every
negative result is reported together with the effect size it would have found.

This repository is the submission companion: it holds the analysis code, the
derived data behind every figure and table, the figures themselves and the
preregistrations — everything needed to check every number in the manuscript,
which is under submission and
therefore not part of this repository.

---

## The result in one paragraph

Across 18 published perturbation data sets, only **2 pass a calibration on the
textbook lineage markers of the axis they claim to model** (7 of 14 when
resolved to individual donor cells) — yet an exploratory, data-derived
173-gene matrix programme, subsequently frozen for downstream use, runs above
its own detection limit in **18 of 18** of them (z +5.25 to +13.10),
as strongly where the calibration fails (z +13.13, n = 16) as where it passes
(z +12.79, n = 2). Under a rule fixed in advance, 10 of the 18 data sets are
decomposable into their separate steps — undifferentiated state left, lineage
not reached, module running — giving **8 confirmations, 2 instances of the
other case and 0 refutations**; no data set leaves the undifferentiated state
without the module running. Re-derived from the studies it never saw
(leave-one-study-out), the programme still clears its own detection limit in
**14 of 18** held-out data sets (16 of 18 leave-one-dataset-out), stands at
z +6.19 above 10,000 expression-, length- and constraint-matched random
173-gene sets, and no single gene carries it (leave-one-gene-out rho
0.616–0.633). Locked and scored on independent data it played no part in
selecting, it runs above its limit in **3 of 4** external differentiation
cohorts (osteogenic z +10.1, adipogenic z +8.2, vascular z +9.0; the
iPSC-derived chondrogenic model does not reach the limit). In chromatin accessibility the module runs on the
adipogenic lineage, on an axis that passes its own calibration. In the human
fetal growth plate it tracks maturation to the prehypertrophic zone
(rho +0.456, z +4.80, limit rho 0.274), and it does not hinge on the single
hypertrophic point (without it rho +0.430 against its own recomputed limit of
0.250). The dysplasia genes sit in the distal secretory machinery (OR 2.84,
z +6.30), for which the programme is itself depleted; they split cleanly by
gene dosage (LOEUF 0.283 against 0.826, P = 6e-24); and they are ordinary on
the dynamics axis, z −0.71 at a detection limit of 0.073 dWT units, in the same
run in which the programme sits at z +18.10.

---

## Reproduce it

Every script derives the repository root from its own location; set
`PAPER_V2_ROOT` only if you want to run against a copy somewhere else.

```bash
python 03_lineage_calibration/10_calibration_18_datasets.py  # the per-data-set calibration (2 of 18)
python 02_matrix_programme_derivation/31_derive_matrix_programme.py  # reproduce the exploratory 173-gene module
python 03_lineage_calibration/12_calibration_gene_space.py                     # which gene space belongs under it
python 07_in_vivo_growth_plate/12_fetal_donor_trend_test.py                         # donor-stratified in vivo test
Rscript 06_orthogonal_layers/60_gene_sets_build.R                        # fix the broad gene sets (writes once)
Rscript 06_orthogonal_layers/61_gene_set_enrichment.R                               # Figure 2D against broad external sets
python 04_programme_decomposition/10_decomposition_18_datasets.py                     # the decomposition on all 18 data sets
python 07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py                         # atlas pseudobulk, cached (writes once)
python 07_in_vivo_growth_plate/13_fetal_gene_decomposition.py                 # is the in vivo trend broadly carried?
python 07_in_vivo_growth_plate/20_postnatal_growth_plate_annotation.py          # needs the raw archive and scanpy
python 07_in_vivo_growth_plate/21_postnatal_growth_plate_test.py                # postnatal anchor: not calibratable
python 07_in_vivo_growth_plate/14_hypertrophic_zone_sensitivity.py              # terminal-zone sensitivity
python 05_programme_validation/10_heldout_and_robustness.py                          # held-out validation and robustness (Figure 2G,H)
python 05_programme_validation/11_external_differentiation_systems.py                        # external validation on independent data (Figure 2I)
python 09_figures/10_panel_data_main.py                                 # one CSV per main-figure panel
python 09_figures/11_panel_data_supplement.py                            # supplement panels and Tables S1-S14
python 09_figures/12_panel_data_second_cohort.py                                    # the panels of supplementary Figure S9
Rscript 09_figures/20_figures_main.R                               # F1 to F6, PDF and PNG
Rscript 09_figures/21_figures_supplement.R                         # S1 to S9, PDF and PNG
python 09_figures/30_graphical_abstract.py                         # graphical abstract, PDF and PNG
python 10_manuscript_checks/10_check_numbers.py                              # every number against its panel file
```

`10_check_numbers.py` is the self-test: every load-bearing number of the text
and of the legends stands in it as a required value against
`figures/data/*.csv`. If it exits 0, the repository is consistent with the
numbers quoted in the manuscript. (The manuscript-side checks — references
cited in both directions, language rules — live with the manuscript sources
and are not part of this repository.)

Everything from stage `04` downward reads only stored outputs under
`derived_data/` and `results/`; it recomputes nothing, runs in seconds to about
three minutes, and reproduces from this repository alone. Two steps need
external input the first time and are marked in [`00_setup.md`](00_setup.md):
`06_orthogonal_layers/60_gene_sets_build.R` (MSigDB and `org.Hs.eg.db`, writes
the frozen gene sets once) and `07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py`
(the 7.6 GB limb atlas, writes once). The bibliographic tables and the
growth-plate zone-marker sets are frozen in the repository under
`derived_data/reference_tables/`, so no literature-fetching script runs at build
time. Re-running the upstream analyses (stages `01`–`03`, `06`–`08`, and
`data_acquisition/`) additionally needs the roughly 98 GB of raw data under
`data_raw/`, which is not part of this repository; every accession is listed in
`figures/data/TS1_eighteen_datasets.csv`.

Environments: Python 3.12 (`requirements.txt`) and R 4.4 (`r_packages.txt`).
The exact versions used are recorded in the session-information files under
`results/`.

---

## Provenance of the 173-gene module

The module is an internal, data-derived result; it was not an external gene
set and was not predefined before this paper. It was identified in the
explicitly exploratory `20_Exploration` analysis of the same 18 perturbation
data sets and then frozen for every downstream analysis.

`30_gene_level_convergence_build.R` produces the per-gene, per-data-set `dWT` values. The archived
selection step is preserved in section (3) of
`02_matrix_programme_derivation/32_dexamethasone_confounder.R`:

1. retain genes with a valid `dWT` in at least 14 of 18 data sets (12,563 genes);
2. centre each data-set column on its median within that universe;
3. count the positive and negative centred signs per gene; and
4. retain genes with `v / n >= 0.90`, where `v` is the larger sign count.

This gives 173 genes, with 129 expected to increase and 44 to decrease. The
public snapshot of the `dWT` matrix is
`derived_data/reference_tables/20d_dWT_matrix.csv.gz`. The standalone
`20f_convergence_dwt.py` script applies the same rule and verifies gene
identity, `n`, `v` and `ri` against the frozen
`S5_konvergente_gene.csv` table. The selection is exploratory because the
same 18 data sets define it; the subsequent uses treat the resulting table as
fixed and do not reselect genes.

This formation step is distinct from `11_internal_gene_map.py`, whose
15-of-18 filter belongs to its separate continuous gene-map analysis and does
not define the 173-gene module.

---

## Layout

The pipeline is one numbered sequence of stages, each a folder at the
repository root, ordered the way the biology is built up. Every metric is
implemented exactly once, in `00_shared/`; no statistic is computed twice.

```
00_shared/                  the metric, the marker sets and the enrichment
                            test — the single implementation everything calls
01_expression_landscape/    the internal per-gene RNA layer (continuous)
02_matrix_programme_derivation/  derive and freeze the 173-gene programme
03_lineage_calibration/     can the 18 models even reach their own lineage?
04_programme_decomposition/  the three-way decomposition over all 18 data sets
05_programme_validation/    held-out re-derivation, robustness, external data
06_orthogonal_layers/       methylation, chromatin accessibility, H3K27ac
07_in_vivo_growth_plate/    fetal limb atlas and postnatal growth plate
08_disease_gene_orthogonality/  human-genetics anchor, patient and donor phases
09_figures/                 panel data, main and supplementary figures, GA
10_manuscript_checks/       number/reference/language checks and packaging
data_acquisition/           GEO cohort search and screening — metadata only,
                            no analysis, nothing downloaded here
figure_style/               the publication style of the figures and the
                            rules they follow
figures/                    F1-F6 and S1-S9 as PDF and PNG at 600 dpi;
                            figures/data/ holds one CSV per panel and per
                            supplementary table
derived_data/               stored analysis outputs that the pipeline reads
results/                    logs, session information and the self-tests
preregistrations/           every preregistration and protocol, dated and
                            unchanged, including the ones that fell
data_raw/                   raw data, about 98 GB, not in this repository
```

Raw data are not part of the repository; they are public at GEO, ArrayExpress
and EMBL-EBI under the accessions in Supplementary Table 1. The repository
itself is the archive for the code and the derived tables: it is made public
on GitHub at submission, and the URL is carried in the Data and Code
Availability statement of the manuscript.

A note on names: the analysis code, the figure scripts and the delivered
tables are English throughout. The stored intermediate files under
`derived_data/` and `results/` keep the short internal column names that the
analysis scripts read; `column_glossary.csv` in the repository root gives the
English meaning of each of them.

---

## How to read the numbers

Two conventions carry the whole paper.

**Every level has its own positive control.** Before any level may produce a
result it must find the textbook lineage markers of the axis it measures, on
the same scale and against the same null, at z >= 2 — a rule fixed before the
data were seen (preregistration M-D, §6). A level that fails produces neither a
positive nor a negative finding; its result is *not measurable*. The postnatal
growth plate is the live example: its calibration fails because hypertrophic
cells switch the cartilage-matrix programme off, so it carries no verdict in
either direction.

**Every number has a detection limit (MDE80).** That is the smallest true
effect the analysis would detect in 80 % of repetitions, computed from the same
permutation null as the test itself. Supplementary Table 7
(`figures/data/TS7_all_statistics.csv`) lists every statistic in the paper with
its own limit. A negative result without a limit is not reported.

Each analysis is labelled confirmatory (preregistered), exploratory, or
preregistered follow-up, in its figure legend and in Supplementary Table 8. The
main result is exploratory and says so in every legend.

---

## Figures

| | content | status |
|---|---|---|
| **F1** | material, the screens, and the calibration almost nothing passes (2 of 18; 7 of 14) | confirmatory (screens) |
| **F2** | the main result: lineage independence and decoupling (A carries the matched-null control); E chromatin; F the decomposition over all 18 data sets; G held-out re-derivation, H gene-dropout robustness, I external validation on independent data | exploratory (F a preregistered follow-up; G-I exploratory validation) |
| **F3** | the human fetal growth plate, including the terminal-zone sensitivity | exploratory |
| **F4** | where the disease genes are: localisation, dosage, constraint, dynamics | mixed, labelled per panel |
| **F5** | both layers meet at the prehypertrophic transition | exploratory |
| **F6** | levels, estimates and each level's own detection limit, including the levels that carry nothing | confirmatory as bookkeeping |
| **GA** | graphical abstract: the growth-plate curves with both layers peaking at the prehypertrophic transition, the 2/18 against 18/18 decoupling, and the orthogonality of the disease-gene definitions | — |
| **S1-S9** | scale critique · every calibration · external triangulation · patient cohorts · robustness · orthogonal levels · screen detail · day zero and publication matching · the lineage contrast in three cohorts | legends in the manuscript's supplemental information |

---

## What is deliberately not in here

Reporting what fell is part of the argument. The preregistered analyses and
protocols are reproduced unchanged in `preregistrations/`. The 173-gene module
formation described above is explicitly exploratory and is included as
provenance, not as a preregistered claim. No downstream analysis searches for a
new axis or changes the frozen module: convergence counts track
signal-to-noise, donor-resolved lesion-response numbers, day-zero competence,
disease genes as constitutive rather than dynamic, and absolute expression are
reported under their stated analysis status.

---

## Citation and licence

This repository accompanies a manuscript under submission. Citation information will
be added on acceptance;
[`CITATION.cff`](CITATION.cff) already carries the machine-readable form for
the repository itself.

Text, figures and derived data are under CC BY 4.0 ([`LICENSE`](LICENSE)); the
code is under the MIT licence ([`LICENSE-CODE`](LICENSE-CODE)). The primary
data sets remain under the terms of their public archives.

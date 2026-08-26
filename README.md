# What skeletal differentiation models actually measure

A lineage-independent matrix programme runs in every published mesenchymal
differentiation model we could reach — including where the model fails its own
positive control — while the skeletal dysplasia genes are defined by
localisation and gene dosage and lie orthogonal to the differentiation-dynamics
axis the field searches on. Every level carries its own internal positive
control ("calibration") and its own measured detection limit (MDE80), and every
negative result is reported together with the effect size it would have found.

The manuscript source is [`manuscript/MANUSCRIPT.md`](manuscript/MANUSCRIPT.md);
the target journal is *iScience* (Cell Press). The submission package is built
from the sources by `code/64_build_submission.py`.

---

## The result in one paragraph

Across 18 published perturbation data sets, only **2 pass a calibration on the
textbook lineage markers of the axis they claim to model** (7 of 14 when
resolved to individual donor cells) — yet a fixed 173-gene matrix programme
runs above its own detection limit in **18 of 18** of them (z +5.25 to +13.10),
as strongly where the calibration fails (z +13.13, n = 16) as where it passes
(z +12.79, n = 2). Under a rule fixed in advance, 10 of the 18 data sets are
decomposable into their separate steps — undifferentiated state left, lineage
not reached, module running — giving **8 confirmations, 2 instances of the
other case and 0 refutations**; no data set leaves the undifferentiated state
without the module running. In chromatin accessibility the module runs on the
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
python reference_implementations/56_calibration_eighteen.py  # the per-data-set calibration (2 of 18)
python code/29_calibration_gene_space.py                     # which gene space belongs under it
python code/20_in_vivo_donor_test.py                         # donor-stratified in vivo test
Rscript code/24a_gene_sets_v2_build.R                        # fix the broad gene sets (writes once)
Rscript code/24_gene_sets_v2.R                               # Figure 2D against broad external sets
python code/25_decomposition_eighteen.py                     # the decomposition on all 18 data sets
python code/26_in_vivo_pseudobulk.py                         # atlas pseudobulk, cached (writes once)
python code/27_in_vivo_gene_decomposition.py                 # is the in vivo trend broadly carried?
python code/32_postnatal_growth_plate_annotation.py          # needs the raw archive and scanpy
python code/33_postnatal_growth_plate_test.py                # postnatal anchor: not calibratable
python code/34_hypertrophic_zone_sensitivity.py              # terminal-zone sensitivity
python code/50_panel_data.py                                 # one CSV per main-figure panel
python code/51_supplement_data.py                            # supplement panels and Tables S1-S14
python code/52_s9_data.py                                    # the panels of supplementary Figure S9
Rscript code/60_figures_main.R                               # F1 to F6, PDF and PNG
Rscript code/61_figures_supplement.R                         # S1 to S9, PDF and PNG
python code/62_graphical_abstract.py                         # graphical abstract, PDF and PNG
python code/63_key_resources_table.py                        # the Key Resources Table
python code/64_build_submission.py                           # the whole submission package
python code/70_check_numbers.py                              # every number against its panel file
python code/71_check_references.py                           # references, both directions
python code/72_check_language.py                             # language rules of the material
```

`70_check_numbers.py` is the self-test: every load-bearing number of the text
and of the legends stands in it as a required value against
`figures/data/*.csv`. If it exits 0, the repository is consistent with itself.
Two further checks guard the manuscript rather than the numbers:
`71_check_references.py` verifies in both directions that every reference is
cited and every citation has an entry, and `72_check_language.py` enforces the
language rules of the material — English throughout, `undifferentiated` rather
than `naive`, the accession always in parentheses after the information it
belongs to, and no abbreviated identifiers.

The panel and figure steps (`50_` upward) read only stored outputs under
`derived_data/` and `results/`; they recompute nothing and take seconds to
about three minutes in total. Three steps need external input and are marked in
[`code/00_setup.md`](code/00_setup.md): `24a_gene_sets_v2_build.R` (MSigDB and
`org.Hs.eg.db`, writes once), `26_in_vivo_pseudobulk.py` (the 7.6 GB limb atlas,
writes once) and `28_geo_primary_publications.py` (NCBI E-utilities,
bibliography only). Re-running the upstream analyses additionally needs the
roughly 98 GB of raw data under `data_raw/`, which is not part of this
repository; every accession is listed in
`figures/data/TS1_eighteen_datasets.csv`.

Environments: Python 3.12 (`requirements.txt`) and R 4.4 (`r_packages.txt`).
The exact versions used are recorded in the session-information files under
`results/`.

---

## Layout

```
code/                       the pipeline: 00_setup, 20-34 analyses,
                            50-52 panel data, 60-62 figures, 63-64 packaging,
                            70-73 checks
reference_implementations/  one implementation per metric; nothing is
                            computed twice. manuscript/methods holds the
                            reference implementation of the metric itself
figure_style/               the publication style of the figures and the
                            rules they follow
figures/                    F1-F6 and S1-S9 as PDF and PNG at 600 dpi, the
                            graphical abstract GA at 300 dpi (1200 x 1200 px,
                            the iScience requirement);
                            figures/data/ holds one CSV per panel and per
                            supplementary table
derived_data/               stored analysis outputs that the pipeline reads
results/                    logs, session information and the self-tests
preregistrations/           every preregistration and protocol, dated and
                            unchanged, including the ones that fell
manuscript/                 MANUSCRIPT.md, the figure legends, the cover
                            letter, the bibliography and the citation styles
data_raw/                   raw data, about 98 GB, not in this repository
submission/                 build output, not in this repository
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
| **F2** | the main result: lineage independence and decoupling; E chromatin; F the decomposition over all 18 data sets | exploratory (F a preregistered follow-up) |
| **F3** | the human fetal growth plate, including the terminal-zone sensitivity | exploratory |
| **F4** | where the disease genes are: localisation, dosage, constraint, dynamics | mixed, labelled per panel |
| **F5** | both layers meet at the prehypertrophic transition | exploratory |
| **F6** | levels, estimates and each level's own detection limit, including the levels that carry nothing | confirmatory as bookkeeping |
| **GA** | graphical abstract: the growth-plate curves with both layers peaking at the prehypertrophic transition, the 2/18 against 18/18 decoupling, and the orthogonality of the disease-gene definitions | — |
| **S1-S9** | scale critique · every calibration · external triangulation · patient cohorts · robustness · orthogonal levels · screen detail · day zero and publication matching · the lineage contrast in three cohorts | see `manuscript/CAPTIONS_SUPPLEMENT.md` |

---

## What is deliberately not in here

Reporting what fell is part of the argument. All of it is preregistered and
reproduced unchanged in `preregistrations/`: convergence counts as evidence
(they track signal-to-noise, not biology), donor-resolved lesion-response
numbers, day-zero competence, disease genes as constitutive rather than
dynamic, absolute expression as a finding, and any new search for a convergence
axis.

---

## Citation and licence

The manuscript is under preparation for *iScience* (Cell Press). The reference
apparatus is format-agnostic: `manuscript/references.bib` plus a CSL switch
(`manuscript/csl/cell.csl` for iScience, `manuscript/csl/springer-vancouver.csl`
for BMC Genomics), selected with
`python code/64_build_submission.py --style {cell,vancouver}`. Citation
information will be added on acceptance; `CITATION.cff` carries the
machine-readable form.

Text, figures and derived data are under CC BY 4.0 ([`LICENSE`](LICENSE)); the
code is under the MIT licence ([`LICENSE-CODE`](LICENSE-CODE)). The primary
data sets remain under the terms of their public archives.

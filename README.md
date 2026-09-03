# A conserved matrix programme of skeletal differentiation, and the orthogonal disease genes

Code and derived data for the analyses of the manuscript *A conserved,
lineage-independent extracellular-matrix programme organizes skeletal cell
differentiation, orthogonal to the secretion- and dosage-defined
skeletal-dysplasia genes* (under submission at *BMC Genomics*).

A fixed 173-gene extracellular-matrix programme runs in every published
mesenchymal differentiation model we could obtain — including in models that
fail their own lineage positive control — while the skeletal dysplasia genes
are defined by secretory-pathway localisation and gene dosage, and lie
orthogonal to the differentiation-dynamics axis the field searches on.

**Scope.** This repository is the analysis companion: the code, the derived
data, the panel data, the figures and the self-tests. The manuscript text is
not part of it and is not needed to reproduce anything below. Where a step
does need the manuscript sources, this README says so.

---

## Quickstart

```bash
git clone https://github.com/alexschulzcell/Differentiation-paper.git
cd Differentiation-paper

python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt                     # Python 3.12
Rscript install_r_packages.R                        # R 4.4, CRAN + Bioconductor

python reproduce.py
```

`reproduce.py` re-derives every number of the paper, redraws every figure, and
then verifies the result against what it just computed. It reads **only files
committed to this repository** — no raw data, no network — and takes about
fifteen minutes on a normal laptop.

A successful run ends like this:

```
==============================================================================
All 19 steps passed.
The manuscript is consistent with the numbers just computed.
==============================================================================
```

If a step fails, `reproduce.py` prints the failing script and its traceback and
returns a non-zero exit code; nothing after it runs.

Useful variants:

| command | what it does |
|---|---|
| `python reproduce.py --list` | print the steps and what each produces, run nothing |
| `python reproduce.py --only checks` | just re-verify the stored panel data, in seconds |
| `python reproduce.py --from figures` | redraw the figures and re-verify, skip the analysis |
| `python reproduce.py --package` | rebuild the submission package — authoring only, needs the manuscript sources, pandoc and Word |

Every script is standalone and can also be run on its own, from anywhere:

```bash
python 03_lineage_calibration/10_calibration_18_datasets.py
```

Each finds the repository root from its own location. To run against a copy of
the repository somewhere else, set `PAPER_V2_ROOT` to that directory. (The
variable is named for an earlier version of the project and is kept so that
older logs remain valid; it points at this repository.)

---

## What `reproduce.py` runs, and in what order

Three groups, in dependency order.

### 1 · `analysis` — the results

| script | question it answers | headline number |
|---|---|---|
| `03_lineage_calibration/10_calibration_18_datasets.py` | can the 18 published models even reach the lineage they claim to model? | **2 of 18** pass; 7 of 14 donor cells |
| `03_lineage_calibration/11_calibration_sensitivity.py` | does that verdict depend on the thresholds? | no |
| `03_lineage_calibration/12_calibration_gene_space.py` | which gene space belongs under the calibration? | background swap moves z by 0.02 |
| `04_programme_decomposition/10_decomposition_18_datasets.py` | undifferentiated state left · lineage reached · programme running, separately | 10 decomposable → **8 / 2 / 0** |
| `05_programme_validation/10_heldout_and_robustness.py` | does the programme survive being re-derived without the studies it came from? | **14 of 18** held out; matched-null z +6.19 |
| `05_programme_validation/11_external_differentiation_systems.py` | does the locked programme run on data it never saw? | **3 of 4** cohorts |
| `07_in_vivo_growth_plate/12_fetal_donor_trend_test.py` | does it track the human fetal growth plate? | rho +0.456, z +4.80 |
| `07_in_vivo_growth_plate/13_fetal_gene_decomposition.py` | is that trend broadly carried, or by a few genes? | broadly |
| `07_in_vivo_growth_plate/14_hypertrophic_zone_sensitivity.py` | does it hinge on the terminal hypertrophic point? | no; rho +0.430 without it |
| `07_in_vivo_growth_plate/21_postnatal_growth_plate_test.py` | the postnatal anchor | fails its own calibration → carries no verdict |

### 2 · `figures` — panel data, then figures

`09_figures/10`, `11` and `12` reshape the stored outputs into **one CSV per
panel** under `figures/data/`, with English column names. The two R scripts
draw from those CSVs and compute nothing. So to redraw a figure you need only
the panel data, and to check a number you need only open its CSV.

| output | script | panel data |
|---|---|---|
| F1–F6, PDF + PNG at 600 dpi | `09_figures/20_figures_main.R` | `figures/data/F1*.csv` … `F6*.csv` |
| S1–S9, PDF + PNG at 600 dpi | `09_figures/21_figures_supplement.R` | `figures/data/S1*.csv` … `S9*.csv` |
| Tables S1–S14 | `09_figures/11_panel_data_supplement.py` | `figures/data/TS1*.csv` … `TS14*.csv` |
| graphical abstract, 300 dpi | `09_figures/30_graphical_abstract.py` | — |

### 3 · `checks` — the self-test

| script | checks |
|---|---|
| `10_manuscript_checks/10_check_numbers.py` | **177 of 177** load-bearing numbers, each against its panel CSV |
| `10_manuscript_checks/11_check_references.py` | every reference cited and every citation resolved, both directions |
| `10_manuscript_checks/12_check_language.py` | the naming and terminology rules the project follows |

`10_check_numbers.py` is the one that matters: every number the paper leans on
stands in it as a required value against `figures/data/*.csv`. If it exits 0,
the repository is consistent with itself.

Two of these three read the manuscript text where it is available, so they
behave slightly differently here than in the authoring repository:

* `10_check_numbers.py` runs **177** checks here. A further 10 compare figure
  legends against the panel data; those need the manuscript sources and are
  skipped, giving 187 in the authoring repository. The 177 cover every number
  that has a data source in this repository.
* `11_check_references.py` cross-checks the reference list against the running
  text. With no manuscript here it reports that it was skipped and exits 0. The
  result of the last run in the authoring repository stands in
  [`results/reference_check.txt`](results/reference_check.txt): every reference
  cited, every citation resolved, nothing open in either direction.

---

## What is *not* in the quickstart, and why

Three kinds of step are excluded, each for a stated reason. Every stage README
marks its own scripts.

| excluded | why | what stands in for it |
|---|---|---|
| everything marked _needs raw data_ | needs the ~98 GB of public raw data under `data_raw/` | its outputs are committed under `derived_data/` and `results/` |
| `06_orthogonal_layers/60_gene_sets_build.R`, `61_gene_set_enrichment.R` | need MSigDB and `org.Hs.eg.db` | outputs frozen under `results/`; the figure steps read those |
| `07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py` | needs the 7.6 GB limb atlas | its pseudobulk is frozen under `derived_data/` |

Raw data are not part of this repository. They are public at GEO, ArrayExpress
and EMBL-EBI; every accession is listed in
[`figures/data/TS1_eighteen_datasets.csv`](figures/data/TS1_eighteen_datasets.csv).
Obtaining and arranging them is described in [`00_setup.md`](00_setup.md).

The bibliographic tables and the growth-plate zone-marker sets that some steps
once fetched from NCBI are frozen under `derived_data/reference_tables/`
(`geo_primary_publications*.csv`, `growth_plate_zone_markers.csv`), so no
literature-fetching script runs at build time.

Environments: Python 3.12 (`requirements.txt`), R 4.4
(`install_r_packages.R`). `r_packages.txt` is a record of the R versions
actually used, in the `name_version` form `sessionInfo()` prints; the installer
reads it and routes each package to CRAN or Bioconductor. The exact versions
behind the published figures are in the session-information files under
`results/`.

---

## Layout

One numbered sequence of stages, each a folder at the repository root, ordered
the way the biology is built up. Every metric is implemented exactly once, in
`00_shared/`; no statistic is computed twice.

| | |
|---|---|
| `00_shared/` | the concordance metric, the marker sets, the gene maps and the enrichment test — the single implementation everything calls |
| `01_expression_landscape/` | the internal per-gene RNA layer |
| `02_matrix_programme_derivation/` | derive and freeze the 173-gene programme |
| `03_lineage_calibration/` | can the 18 models reach their own lineage? |
| `04_programme_decomposition/` | the three-way decomposition over all 18 data sets |
| `05_programme_validation/` | held-out re-derivation, robustness, external data |
| `06_orthogonal_layers/` | methylation, chromatin accessibility, H3K27ac |
| `07_in_vivo_growth_plate/` | fetal limb atlas and postnatal growth plate |
| `08_disease_gene_orthogonality/` | the human-genetics anchor, patient and donor phases |
| `09_figures/` | panel data, main and supplementary figures |
| `10_manuscript_checks/` | the number, reference and language checks, and packaging |
| `data_acquisition/` | the GEO cohort search and screening — metadata only, nothing downloaded here |
| `figure_style/` | the publication style of the figures and the rules they follow |
| `figures/` | F1–F6 and S1–S9 as PDF and PNG; `figures/data/` holds one CSV per panel and per supplementary table |
| `derived_data/` | stored analysis outputs that the pipeline reads |
| `results/` | logs, session information and the self-tests |
| `preregistrations/` | every preregistration and protocol, dated and unchanged, including those the data did not support |

Each numbered stage folder, and `data_acquisition/`, carries its own README
stating what its scripts do and which of them are runnable from this repository
alone.

A note on names: the analysis code, the figure scripts and the delivered tables
are English throughout. The stored intermediates under `derived_data/` and
`results/` keep the short internal column names the analysis scripts read;
[`column_glossary.csv`](column_glossary.csv) gives the English meaning of each.

---

## How to read the numbers

Two conventions carry the whole paper, and both are enforced in code.

**Every level has its own positive control.** Before any level may produce a
result it must find the textbook lineage markers of the axis it measures, on
the same scale and against the same null, at z ≥ 2 — a rule fixed before the
data were seen (preregistration M-D, §6). A level that fails this control
produces neither a positive nor a negative finding: its result is reported as
*not measurable* rather than as evidence in either direction. The postnatal
growth plate is the worked example. Its calibration fails because hypertrophic
cells switch the cartilage-matrix programme off, so it carries no verdict.

**Every number has a detection limit (MDE80)** — the smallest true effect the
analysis would detect in 80 % of repetitions, from the same permutation null as
the test itself.
[`figures/data/TS7_all_statistics.csv`](figures/data/TS7_all_statistics.csv)
lists every statistic in the paper with its own limit. A negative result
without a limit is not reported.

Each analysis is labelled confirmatory (preregistered), exploratory, or
preregistered follow-up, in its figure legend and in Supplementary Table 8. The
main result is exploratory and says so in every legend.

---

## The figures

| | content | status |
|---|---|---|
| **F1** | the data sets, the screens, and the calibration almost nothing passes (2 of 18; 7 of 14) | confirmatory (screens) |
| **F2** | the main result: lineage independence and decoupling. A carries the matched-null control; E chromatin; F the decomposition over all 18 data sets; G held-out re-derivation; H gene-dropout robustness; I external validation | exploratory (F a preregistered follow-up) |
| **F3** | the human fetal growth plate, including the terminal-zone sensitivity | exploratory |
| **F4** | where the disease genes are: localisation, dosage, constraint, dynamics | mixed, labelled per panel |
| **F5** | both layers meet at the prehypertrophic transition | exploratory |
| **F6** | every level with its estimate and its own detection limit, including the levels that carry nothing | confirmatory as bookkeeping |
| **S1–S9** | scale critique · every calibration · external triangulation · patient cohorts · robustness · orthogonal levels · screen detail · day-zero competence and publication matching · the lineage contrast in three cohorts | full legends in the published supplement |

---

## Analyses that did not survive

Reporting what did not hold is part of the argument. Each of these was
preregistered and is reproduced unchanged in `preregistrations/`, together with
the result that retired it: convergence counts as evidence (they track
signal-to-noise, not biology), donor-resolved lesion-response numbers, day-zero
competence, disease genes as constitutive rather than dynamic, absolute
expression as a finding, and any new search for a convergence axis.

---

## Citation and licence

The manuscript is under submission at *BMC Genomics*. Citation information will
be added on acceptance; [`CITATION.cff`](CITATION.cff) carries the
machine-readable form.

Text, figures and derived data are under CC BY 4.0 ([`LICENSE`](LICENSE)); the
code is under the MIT licence ([`LICENSE-CODE`](LICENSE-CODE)). The primary
data sets remain under the terms of their public archives.

Questions about the analyses, and anything that does not reproduce, are best
raised as an issue on this repository.

# 00_setup — environment, paths and data acquisition

## 1 · Paths

No path is hard-coded. Every script derives the repository root from its own
location, so a clone works wherever it is put. Two environment variables
override that, for the case where the code and the data live apart:

| variable | meaning | default |
|---|---|---|
| `PAPER_V2_ROOT` | the repository root | the parent directory of the script |
| `LIMB_ATLAS` | path to the human fetal limb atlas `.h5ad` | inside `data_raw/` |

```bash
export PAPER_V2_ROOT="/path/to/repo"
```

## 2 · Software

R 4.4.3 — `ggplot2`, `patchwork`, `ragg`, `systemfonts`, `DESeq2`,
`matrixStats`, `rtracklayer`, `GenomicRanges`.
Python 3.12 — `numpy`, `pandas`, `scipy`.

Every script writes `sessionInfo()` or the package versions into its log; the
figure scripts write a session-information file into `results/`.

## 3 · Raw data

**The 96 GB of raw data are not in this repository.** They are public. Every
series used is listed with its accession in Supplementary Table 1
(`figures/data/TS1_eighteen_datasets.csv`); every screened series, including
the excluded ones with their exclusion code, is in Supplementary Table 2.

Reference data with retrieval dates:

| source | what | retrieved |
|---|---|---|
| GENCODE v46 (hg38) | gene annotation, union exon lengths | — |
| gnomAD | LOEUF constraint | — |
| Genomics England PanelApp | skeletal dysplasia panels (309 green; 1 471 broad) | — |
| Nosology of Genetic Skeletal Disorders | core and broad panels | — |
| NCBI `gene2pubmed` and `gene2ensembl` | publications per gene, tax_id 9606 | **2026-08-22** |
| Human fetal limb atlas | post-conception weeks 5.1–9.3, 136 311 cells | — |

Processed matrices are deposited at Zenodo; the DOI goes here on submission.
**DOI placeholder: `10.5281/zenodo.XXXXXXX`** (a snapshot of this
repository).

## 4 · Order of execution

```bash
python reference_implementations/56_calibration_eighteen.py  # the per-data-set calibration, 2 of 18
python code/29_calibration_gene_space.py                     # which gene space belongs under it
python code/20_in_vivo_donor_test.py                         # the donor-stratified in vivo test
Rscript code/24a_gene_sets_v2_build.R                        # fix the broad gene sets (writes once)
Rscript code/24_gene_sets_v2.R                               # Figure 2D against broad external sets
python code/25_decomposition_eighteen.py                     # the decomposition on all 18 data sets
python code/26_in_vivo_pseudobulk.py                         # atlas pseudobulk, cached (writes once)
python code/27_in_vivo_gene_decomposition.py                 # is the in vivo trend broadly carried?
python code/32_postnatal_growth_plate_annotation.py          # postnatal zones (needs the raw archive)
python code/33_postnatal_growth_plate_test.py                # postnatal anchor: not calibratable
python code/34_hypertrophic_zone_sensitivity.py              # does the fetal anchor hang on the terminal zone?
python code/50_panel_data.py                                 # one CSV per main-figure panel
python code/51_supplement_data.py                            # supplement panels and Tables S1-S14
python code/52_s9_data.py                                    # the panels of supplementary Figure S9
Rscript code/60_figures_main.R                               # F1 to F6, PDF and PNG
Rscript code/61_figures_supplement.R                         # S1 to S9, PDF and PNG
python code/62_graphical_abstract.py                         # graphical abstract, PDF and PNG
python code/70_check_numbers.py                              # every number against its panel file
```

Submission packaging and the manuscript-side checks (`71_check_references.py`,
`72_check_language.py`, `73_primary_references.py`) are not part of this
public repository.

The panel and figure steps read only files that already exist under
`derived_data/` and `results/`; they recompute nothing and take seconds. Steps
`20_`, `24_`, `25_` and `27_` recompute one analysis each from stored
per-sample or per-gene values and take seconds to a few minutes.

Two steps are **exceptions to "reads only stored files"** and are marked as
such:

| step | needs | note |
|---|---|---|
| `24a_gene_sets_v2_build.R` | `msigdbr`, `org.Hs.eg.db`, `GO.db` | writes the frozen gene sets **once** and refuses to overwrite them; set `GENSAETZE_V2_NEU=1` to rebuild |
| `26_in_vivo_pseudobulk.py` | the 7.6 GB limb atlas `.h5ad` and about 3 GB of memory | writes `results/invivo_pseudobulk.csv.gz` **once**; set `INVIVO_PSEUDOBULK_NEU=1` to rebuild. Path override: `LIMB_ATLAS` |

`code/32_postnatal_growth_plate_annotation.py` additionally needs the unpacked
GEO archive of the postnatal growth-plate series and the frozen zone-marker
sets under `data_raw/`, plus scanpy with leidenalg. Its decision rule stands in
its own header, before any number. The outcome is fixed by guards in
`70_check_numbers.py`: that level is **not calibratable**, and should a re-run
pass its positive control, those guards fire and the manuscript text has to be
revisited.

Both cached outputs are checked against their own source on every use: `25_`
asserts that its own-axis contrast reproduces the stored calibration
character for character, and `27_` asserts the same against the stored
per-sample module values. If either assertion fails, the script stops.

Re-running the upstream analyses — the ones that touch raw data — needs the raw
data; those scripts are in `reference_implementations/`, each with its own
header.

## 5 · Deposit

The repository is archived on Zenodo on submission. After deposition the DOI
is entered here and in the *Data and code availability* section of the
manuscript:

- **DOI 10.5281/zenodo.XXXXXXX (assigned on deposition)**

The deposit contains no raw data; the accessions are in Supplementary
Table 1.

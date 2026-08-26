# Follow-up search extension: independent donor cohort

Date: 2026-08-20
Status: locked before the extension search and before opening any new
expression matrix.

## Purpose

The original follow-up search yielded two eligible accessions,
`GSE247491` and `GSE247528`, but both share the same three donors. This
extension attempts to find a second independent donor cohort for the fixed
173-gene differentiation programme and, if available, an independent isogenic
2 x 2 lesion study.

This is a search extension, not a license to add a result-selected dataset. A
candidate is eligible only by metadata and biological-unit rules below. Its
expression result is not inspected before classification.

## Fixed extension queries

Use GEO/NCBI GDS Entrez metadata search, with at least 200 hits per query when
available, and deduplicate accessions:

1. `("mesenchymal stem cell" OR hMSC OR BMSC) AND (osteogenic OR osteoblast*) AND (RNA-seq OR transcriptome OR expression) AND (donor OR primary OR patient) AND Homo sapiens[Organism]`
2. `("mesenchymal stem cell" OR hMSC OR BMSC) AND (chondrogenic OR chondrogenesis OR cartilage) AND (RNA-seq OR transcriptome OR expression) AND (donor OR primary OR patient) AND Homo sapiens[Organism]`
3. `human BMSC osteogenic day0 day7 donor RNA-seq`
4. `human BMSC chondrogenic day0 donor RNA-seq`
5. `human MSC undifferentiated differentiated control perturbation siRNA RNA-seq`
6. `human mesenchymal stem cell 2x2 control knockdown differentiation RNA-seq`

Record query, date, hit count, accession, title, summary, platform, sample
count, and URL. Search metadata only in this phase.

## Eligibility before expression analysis

- Human primary MSC/BMSC or equivalent primary skeletal stromal cells.
- Clear undifferentiated versus differentiated axis in the same lineage.
- At least three independent donors or biological sources per relevant state.
- Donor/sample identity recoverable from GEO metadata; cells and technical
  wells are never biological units.
- No accession already analysed in S11, S12, or the SERPINA3 cohort.
- No shared donors, cell source, or starting matrix with the SERPINA3 cohort.
- For a lesion candidate, a complete same-unit undifferentiated/differentiated x
  control/perturbed 2 x 2 is required; perturbation direction is classified
  from metadata before expression analysis.

Exclude cell lines, immortalized/TERT/iPSC-only systems, pooled donors without
unit resolution, missing undifferentiated arms, unclear donor identity, non-human data,
and studies without the required expression design. Keep excluded candidates in
the ledger with an exclusion code and evidence URL.

## Analysis lock

If an eligible independent cohort exists, use the unchanged S5 173-gene module
and `ri`. Use the corrected transcript-to-gene aggregation, not the previous
`drop_duplicates` implementation. Primary atlas inference is one study-level
endpoint per independent cohort with whole-donor sign flips. Do not pool
gene-by-dataset rows as independent replications.

If no second independent cohort exists, report the search as negative. Do not
promote the two SERPINA3 lineages to two independent studies.

## Outputs

- `26_Orthogonal_S12/data_raw/followup_extension_search_log.csv`
- `26_Orthogonal_S12/data_raw/followup_extension_search_raw.csv`
- `26_Orthogonal_S12/data_raw/followup_extension_kandidaten.csv`
- `26_Orthogonal_S12/data_raw/_download_protokoll_followup_extension.md` only if a
  candidate passes metadata eligibility
- an update to `BERICHT_FOLLOWUP.md` documenting whether a second independent
  cohort was found

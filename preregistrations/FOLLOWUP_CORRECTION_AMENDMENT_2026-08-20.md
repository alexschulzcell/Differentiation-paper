# Follow-up correction and extension amendment

Date: 2026-08-20
Status: locked before corrected expression parsing and before any new candidate
expression result is inspected.

## 1. Corrections to frozen follow-up implementation

The audit found two implementation defects that are independent of the result:

1. Salmon transcript rows were deduplicated instead of summed by mapped gene.
   The corrected operation is `groupby(gene).tpm.sum()` for every sample.
2. The preregistered donor module-score sensitivity endpoint standardized each
   gene separately within each state and then tested the same donors. Its state
   mean is therefore identically zero up to floating-point error. It is marked
   **not estimable under the registered normalization** and will not be used as
   evidence or plotted as a biological null result.

The fixed module, `ri`, primary atlas endpoint, 2x2 `iv` endpoint, time points,
donor definitions, and existing seeds remain unchanged. The corrected existing
SERPINA3 primary analyses rerun with the registered seeds 20260825 and 20260827.

For the three newly eligible microarray cohorts, platform annotations are
locked before expression parsing: retain probes mapping unambiguously to one
fixed-module gene symbol; if multiple retained probes map to the same gene,
average their processed expression values within sample. Probes with multiple
gene symbols or no fixed-module symbol are excluded. The GEO series matrix is
used as deposited, with the mapping table and hash recorded.

The literal S12 rule that missing module genes remain in the denominator is
also enforced for all corrected external atlas outputs. The denominator is
therefore the fixed 173 genes; unavailable or incomplete genes contribute a
non-match and are reported separately as coverage. The earlier measurable-only
outputs are discarded as audit outputs, not paper results.

## 2. Extension candidates locked after metadata-only classification

The search extension was locked in
`FOLLOWUP_SEARCH_EXTENSION_2026-08-20.md`. Metadata eligibility identified:

- `GSE18043`: 3 independent bone-marrow MSC donors, Day0/Day7 osteogenic axis;
- `GSE63754`: 3 independent adipose-derived MSC donors, paired control/osteo;
- `GSE12266`: 4 independent bone-marrow hMSC patients, confluent/Day7
  mineralization axis.

`GSE12265` is the same four-patient historical series as `GSE12266` and is not
counted separately. No expression result was used for this classification.

New atlas permutation seed: **20260828**. A new study-level joint-null
synthesis, if run, uses seed **20260829** and is labelled exploratory because
it was not the original S12 primary endpoint.

## 3. New inferential rules

- One study/cohort contributes one primary atlas effect.
- No gene-by-dataset row is treated as an independent study.
- The new candidates are analysed as dWT-only atlas datasets; none is silently
  treated as an external `iv` study.
- The chondrogenic and the osteogenic SERPINA3 arm remain two lineage
  endpoints from one shared donor cohort.
- A second independent donor cohort can support a stronger atlas statement only
  after the corrected primary endpoint and threshold rule are met; a directional
  result alone is not confirmation.
- The study-level exploratory synthesis uses equal study weights and a joint
  whole-unit sign-flip null. GSE37521 remains a sensitivity dataset because it
  has only two biological units.

## 4. Required provenance

All corrected matrices, platform mappings, source URLs, file sizes, SHA256
hashes, seeds, candidate verdicts, and rerun logs must be written under
`26_Orthogonal_S12/` before the results are integrated into the paper.

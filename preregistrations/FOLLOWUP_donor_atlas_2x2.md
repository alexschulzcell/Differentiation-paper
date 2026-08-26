> Translated from the German original of 2026-08-20. The content, the dates
> and every number are unchanged.

# Preregistration FOLLOW-UP — a donor-resolved differentiation atlas and an independent isogenic 2 x 2 lesion study

Date: **2026-08-20**
Status: written **before** any result was looked at and before any candidate
was downloaded.
Related: `PREREG_S12_orthogonal_triangulation.md`,
`derived_data/reference_tables/S5_konvergente_gene.csv`.

## 1. Purpose and scope

This session attempts a **stronger external test** of the two core statements
of the paper:

1. **The differentiation programme:** a donor-resolved, independent human
   transcriptome atlas is to test the fixed 173-gene programme.
2. **The lesion response (`iv`):** an independent isogenic 2 x 2 perturbation
   study is to test the lesion interaction separately.

Forbidden: any change to the fixed module
`derived_data/reference_tables/S5_konvergente_gene.csv` (173 genes, `ri`), to
`kern()`, to `figure_style/publication_style.R`, to the S11 or S12 results, or
to any existing null. There is no new gene list, no post-selection and no
threshold-dependent relabelling.

**The inferential unit:** the donor or biological unit, and the study.
Gene-by-data-set rows may **never** be treated as independent study
replication. A pooled or sensitivity finding that does not use the
preregistered unit is explicitly **exploratory**.

## 2. Candidate search and metadata classification (before any download or result)

Fixed search terms (GEO/NCBI gds, plus ArrayExpress and ENCODE where
available):

- `human MSC donor undifferentiated differentiated RNA-seq`
- `human mesenchymal stem cell osteogenic differentiation donor transcriptome`
- `human MSC chondrogenic differentiation biological replicates`
- `human MSC isogenic perturbation undifferentiated differentiated control RNA-seq`

The procedure per candidate, **before** any expression or peak analysis:
1. record the title, organism, assay, cell source, study type and donor count
   from the metadata;
2. record the undifferentiated against differentiated status and the exact time
   point;
3. record the control and perturbation arms, if present;
4. document biological against technical replicates (the authors' statement,
   the GSM structure);
5. assess independence from S11 and S12 (cell source, donors, source matrix);
6. assign an inclusion or exclusion code **before** any computation; record the
   source URL and, where applicable, the download name, size and SHA-256 in the
   candidate ledger.

Neither the title nor the track record acts as a result selection: a candidate
is **not** taken up because of a probable positive result, and a negative
appearance does not lead to exclusion. Cell lines, technical replicates, pooled
donor samples or unclear donor identity are **not** strong donor-resolved
carriers.

## 3. Inclusion and exclusion rules (the atlas route; independent of the result)

Inclusion for the donor-resolved atlas analysis:

- a human bulk RNA-seq, microarray or scRNA-seq data set with a clear
  **undifferentiated** and **differentiated** axis of the same differentiation
  lineage (osteogenic or chondrogenic);
- at least **3 independent biological units** (donors or donor cell sources;
  for scRNA-seq, donor or sample pseudobulk) per relevant time point;
- the biological unit is the donor, sample or cell source; cells and wells of
  the same preparation are not units.

Exclusion:

- cell lines, immortalised or clonal lines as the only source (as in S12: hFOB,
  TERT lines, iPSC clones without a donor structure);
- only technical replicates, or unclear donor identity;
- pooled donor matrices without resolution into units;
- mouse or other species (not Homo sapiens);
- data sets already analysed in S11 (`GSE200492`, `GSE166824`, `GSE324998`,
  `GSE255646`) as well as all analysed in S12 (`GSE185951`, `GSE161176`,
  `GSE113253`, `GSE210984`, `GSE202080`, `GSE37521`, `GSE125166`, `GSE220162`,
  `GSE286540`) do not count anew;
- if two studies use the same donors or the same source matrix, only the one
  classified first counts as carrying evidence (noted in the independence
  record).

## 4. The definition of the donor or biological unit

Per data set, the biological unit is fixed from the metadata **before** the
module computation:

- bulk with n labels per state that are declared as n independent donors and
  differ in the GSM metadata → n biological units;
- scRNA-seq → pseudobulk: the sum or aggregate of all cells of a sample or
  donor at one time point; cells are never units;
- several wells of the same preparation → averaged within the unit;
- an unclear assignment → exclusion or purely descriptive status.

## 5. The primary endpoint — the donor-resolved transcriptome atlas

For every eligible independent study **one** primary endpoint is fixed and is
**not** exchanged for another after the result:

- **P = `share_match`** (the unchanged route A endpoint from S12 and S11): the
  share of the measurable 173 genes whose direction, differentiated minus
  undifferentiated on the unit-averaged gene difference, reproduces the S5
  `ri`. Measurable means module genes that occur in the data set and have a
  value in all units (complete cases). Missing genes stay in the denominator.
  The earliest clearly differentiated time point is fixed and logged **before**
  the result is looked at.

- **The assessment of the atlas:** the differentiation atlas carries the
  programme only if **at least two independent studies** meet the strong route
  A carrier status (the rule below). A single study, or a tendency, is reported
  as an indication, not as confirmation.

### The secondary (preregistered) endpoint — donor module-score sensitivity

Only as a **secondary, explicitly named sensitivity analysis**:

- per biological unit and state the signed module score is formed:
  `score = mean( (x_g - mean_g)/sd_g * ri_g )` over the module genes measurable
  in the data set (a per-gene z standardisation within the state is forbidden;
  instead, log2 expression values standardised per gene over all units and
  multiplied by `ri`);
- the test: `diff_score = score(differentiated) - score(undifferentiated)` per
  unit; a Wilcoxon signed-rank or paired mean test over the units, with p and
  the effect. The unit is the donor. (With several time points: the earliest
  differentiated one, identical to P.)
- this endpoint **cannot** rescue a missed primary endpoint and never runs as a
  primary study replication (the inferential unit in P remains the study; here
  it is the donor level within one study).

## 6. The primary endpoint — an independent isogenic 2 x 2 lesion study

Preconditions for an eligible 2 x 2 candidate:

- an undifferentiated and a differentiated state **and** a control and a
  perturbation arm;
- the same biological unit (donor or sample) in as many of the four cells as
  possible; at least **3 independent biological units per cell**, or a
  defensible paired design (the same unit in all four cells);
- the perturbation is chosen and named from the study metadata **before** the
  fixed-module result (the direction follows logically from the title and
  background, not from the result);
- isogenic: the perturbation disturbs the differentiation of the same cellular
  background; not a comparison of different cell lines as the perturbation.

The primary endpoint:

- `dWT` = the differentiation programme in the control arm: per gene and per
  unit, differentiated minus undifferentiated in the control arm; the study
  endpoint as for the atlas P = `share_match` against the fixed `ri`.
- `iv` = the perturbation-by-differentiation interaction: per gene and per unit,
  `(perturbed_diff - perturbed_undiff) - (control_diff - control_undiff)`. The
  endpoint is a convergence test against chance: the number of genes with a
  consistent `iv` sign in **one direction** over the units, against the
  permutation expectation (as in the S11 `iv` count: the same 90 % rule, n
  against expectation; the expectation from whole-unit sign flips).
- no new genes, no post hoc markers, no result-dependent choice of arm; the
  fixed module and the fixed `ri` from S5.

If no complete independent 2 x 2 exists: **no externally derived lesion null**;
the branch is reported as `not assessed` and the phrase "lesions do not"
remains an **internal** finding.

## 7. Null model, permutations, seeds, MDE80 and leave-one-out (identical to S12)

- the null for route A and for the 2 x 2 `dWT`: 10 000 permutations with a
  **joint sign flip of the complete biological unit or sample columns** (the
  sample structure is preserved; genes are never permuted individually);
- seeds (fixed, before the computation): metadata classification `20260820`,
  permutation `20260825` (atlas), `20260826` (the secondary donor module
  score), `20260827` (the 2 x 2 `iv` permutation);
- reported per study and endpoint: the share, the Wilson 95 % interval, the
  null mean, the null SD, z (z = (effect − null mean)/null SD), the two-sided
  permutation p value, **MDE80 = null mean + 2.80 x null SD** (identical to
  S12), and a leave-one-unit-out per study;
- MDE80 is **not** changed after the result.

## 8. Minimum replication requirements and the strong-carrier rule

- the atlas route: **at least 3 independent biological units per time point**
  for a strong carrier; with fewer than 3 units it is not confirmatory;
- **a strong carrier (atlas)**: at least 3 units per time point **and**
  `share > MDE80` **and** a permutation p of at most 0.05 **and** no single
  leave-one-out step bringing the share to or below the null mean;
- the atlas counts as confirmed only with **at least 2 independent studies** as
  strong carriers; otherwise the route is reported as **not confirmatory** and
  directional tendencies as tendencies;
- for the 2 x 2: the `dWT` endpoint under the same carrier rule; the `iv`
  endpoint only as convergence against chance with p and its expectation (no
  strong-carrier rule — that test is the absence test).

## 9. Missing genes, pseudobulk, technical replicates and several time points

- **Missing genes** remain in the denominator only as "not measurable"; no
  imputation and no symbol-based filling of duplicates (technical base level
  only; multiple hits are not double counted).
- **Pseudobulk:** for scRNA-seq, aggregation within donor and time point
  **before** any module computation.
- **Technical replicates:** averaged within the unit; they never count as
  biological units.
- **Several time points:** the earliest clearly differentiated one is logged
  before any result is visible, identically for all studies of one route.
- Where the time-point decision remains ambiguous: the candidate is descriptive
  only, not confirmatory.

## 10. Stopping criteria

- **AB1:** no two independent studies reach strong carrier status → report the
  atlas as not confirmatory; the 2 x 2 branch separately.
- **AB2:** an analysis needs the gene list, the replicate rule, a threshold, the
  time-point rule, a population or a reference to be relaxed after the result is
  known → stop the session, document the conflict, stop the branch, and
  continue only with non-conflicting work.
- **AB3:** a robust external finding contradicts a statement of the paper →
  report the finding; no autonomous manuscript or submission decision.
- **AB4:** no layer has a sufficient biological replicate structure → report the
  limitation; do not upgrade any level of evidence.

## 11. Planned changes to the figures and the manuscript

- At most 6 main figures; **no Figure 7**. Priority: update **Figure 6 in
  place**:
  - **Panel A:** the donor-resolved atlas forest (one row per study, with
    study-level thresholds; a separate descriptive aggregate only where
    justified).
  - **Panel B:** the donor module-score sensitivity (secondary, explicitly
    labelled; it rescues no missed primary endpoint).
  - **Panel C:** the independent 2 x 2 `iv`, only if an eligible study exists.
  - **Panel D:** the updated internal and external evidence map; untested cells
    explicitly `not assessed`, never `chance`.
- Panel numbers are read by the scripts from CSV files, never set by hand.

## 12. Rules of language (adopted)

- `confirmatory` only where the strong-carrier rule is met.
- Pooled and sensitivity results are `exploratory` unless this preregistration
  says otherwise and the inferential unit is correct.
- Route A is **transcriptomic replication**, never "orthogonal".
- B, C and D are "orthogonal" only where the assay or phenotype is genuinely
  different.
- An external lesion response is never called "at chance" unless an eligible
  external 2 x 2 has been analysed. Without an external 2 x 2, "differentiation
  converges, lesions do not" stays qualified as an **internal** lesion finding.

## 13. Addendum 2026-08-20 (before the first expression download or result) — the SERPINA3 cohort made concrete

For the studies `GSE247491` (chondrogenesis) and `GSE247528` (osteogenesis),
classified in phase 1 as the only eligible donor-resolved cohort, the following
is fixed **before any download of expression data and before any result is
looked at**:

1. **The biological unit:** the donor. The donor coding comes from the GEO SOFT
   metadata (`donor:` in `!Sample_characteristics_ch1`) per GSM. The real
   assignment is permuted against the replicate suffix per study:
   - GSE247491: 2454E = Control_repl1 + KD_repl2; 2802F = KD_repl1 +
     Control_repl3; 8A = Control_repl2 + KD_repl3 (all time points).
   - GSE247528: 8A = both repl1; 2454E = both repl2; 2802F = both repl3.
   There is therefore a complete 2 x 2 cell per donor (D0/D7 x control/KD) in
   both studies.
2. **The primary time point:** the earliest clearly differentiated one is
   **day 7** for both lineages (chondrogenesis and osteogenesis), fixed before
   the result. Day 3 only as a **sensitivity** (never primary, never as a
   separate carrier).
3. **The scale:** per sample `log2(TPM + 1)`; the Salmon abundance files
   (transcript- or gene-level identifiers, after inspecting the file format
   without looking at module gene values) aggregated to gene-level Ensembl
   (sum, nearest hit); module genes mapped through the Ensembl gene identifier.
4. **`dWT` (the atlas endpoint and the 2 x 2 `dWT`):** per gene and per donor,
   `Control_day7 − Control_day0` (a difference on the log2(TPM+1) scale). The
   study endpoint is `share_match` against the fixed `ri`; permutation with
   10 000 whole-donor sign flips, seed `20260825` (atlas, §7).
5. **`iv` (the 2 x 2 lesion endpoint):** per gene and per donor,
   `(KD_day7 − KD_day0) − (Control_day7 − Control_day0)`. The endpoint is the
   number of measurable module genes with a **unanimous** `iv` sign across all
   3 donors (with n = 3 units, unanimity is the only non-arbitrary 90 % level),
   against the permutation expectation (10 000 whole-donor sign flips of the
   `iv` columns, seed `20260827`). In addition, convergence by a majority rule
   (2 of 3) as a sensitivity.
6. **A note for later:** GSE247491 and GSE247528 share the same 3 donors → by
   §3 and §8 they are only **one** independent donor cohort (two programmes).
   The atlas finding needs at least 2 independent strong carriers (§5, §8);
   with one cohort the atlas route can therefore at most be classified as a
   strong programme finding per lineage, but **not** as a two-study
   confirmation. Both programmes are nonetheless run with full numbers.
7. **Missing genes and quality control:** identical to §9; complete cases per
   endpoint; no symbol-based filling and no double counting.

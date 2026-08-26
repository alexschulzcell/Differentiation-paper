> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# Preregistration S9 — a new R2 candidate for a cell-type-specific lesion response

Written and dated **2026-08-19**, **before** the new search query and before
any new download. It builds on the S8 brief (in force) and on the decision of
2026-08-19 to strengthen the work through a load-bearing R2 finding (a
cell-type-specific lesion response).

## 1. Scope and honesty

- **S8 is complete and remains binding.** `PREREG_S8.md` (including S8-AB1 and
  S8-AB2) and the S8 report are **not** revoked. S9 is a **new, prospective,
  preregistered** session; it changes no number of the existing paper and no
  main figure.
- A cell-type-specific lesion response would, if robust, be an **exploratory
  new narrative in need of independent confirmation**. It would not quietly
  push the bulk null statement aside but would be reported as an observation of
  its own, to be tested independently.
- A cell count does not replace biological replicates; cells of the same sample
  are not units.

## 2. What this session is looking for

A clinical single-cell data set that allows a lesion response to be tested
**within stably annotated populations**: diseased or lesioned against control
or healed, with **at least two independent biological sample units per group**
(samples, not cells), and whose cell types can be assigned stably (at least
80 %) between **two external references**.

Candidate indications (not exhaustive): fracture non-union or pseudarthrosis,
osteonecrosis, osteoarthritis, bone tumour or pre-tumour, meniscus, tendon or
ligament (including OPLL), and joint or bone defect healing, provided a
controlled or healed arm exists.

## 3. The single new search query (frozen), dated 2026-08-19

```text
("Homo sapiens"[Organism]) AND
("Expression profiling by high throughput sequencing"[DataSet Type]) AND
("single-cell"[All Fields] OR "scRNA"[All Fields] OR "scRNA-seq"[All Fields] OR "single cell RNA"[All Fields]) AND
(bone[All Fields] OR skeletal[All Fields] OR fracture[All Fields] OR nonunion[All Fields] OR pseudarthrosis[All Fields] OR osteonecrosis[All Fields] OR osteoarthritis[All Fields] OR osteoblast*[All Fields] OR chondrocyte*[All Fields] OR cartilage[All Fields] OR ligament[All Fields] OR meniscus[All Fields] OR tendon[All Fields]) AND
(control[All Fields] OR healthy[All Fields] OR healed[All Fields] OR normal[All Fields] OR patient[All Fields] OR disease[All Fields] OR injury[All Fields] OR repair[All Fields])
```

Run with NCBI E-utilities `esearch` and `esummary` on `gds` up to
`retmax=1000`. All hits are screened, not only the first by relevance;
accessions already screened are marked `A8` and are not counted twice.

## 4. Screening and inclusion rules for download

A data set is a **candidate for download** only if all of the following are
documentable from the GEO metadata (and, where necessary, from the sample
records) **before the download**:

1. a human single-cell RNA sequencing design with an accessible matrix or
   counts and sample metadata;
2. a **diseased or lesioned** arm and a **control or healed** arm;
3. **at least two independent biological samples per group** (donor, patient or
   sample), with cells of the same sample not substituting for units;
4. no sample overlap with the 18 existing data sets and no design already
   excluded or already screened in S8 or its addenda.

A data set is not taken up even if the matrix is large or a result could
support the narrative. Missing sample replicates, a missing controlled or
healed arm, or unclear sample units are logged with the concrete reason.

## 5. Analysis after download (only for candidates taken up)

- Cell-type annotation with **both** external references: (a) the six MSigDB C8
  signatures from `PREREG_S7.md` §3; (b)
  `celldex::HumanPrimaryCellAtlasData` (`label.main`). A cell counts as stable
  if both references place the same label on a common coarse ontology; a
  population is evaluable if it reaches at least 10 % per sample and at least
  100 cells across both groups **and** the cross-reference agreement is at
  least 80 %.
- The R2 statistic: lesion against control per population, with the **sample**
  as the unit of comparison (median, range, a Wilson sample confidence
  interval). NO `iv`, NO 2 x 2 interaction term. "No difference" only with a
  predefined equivalence test and MDE80; otherwise "not decidable".
- NO `dWT`/`iv` main metric; no new gene discovery.

## 6. Stopping criteria

- **S9-AB1:** no hit satisfies the inclusion rules (§4). Report; the aim is
  then not reachable through an R2 component.
- **S9-AB2:** for no candidate is the cell-type annotation stable between the
  two references (at least 80 %). No cell-type claim.
- **S9-AB3:** a planned analysis needs a threshold, gene list, population or
  comparison group to be relaxed after a result is known. **Do not decide
  autonomously; stop the session.**
- **S9-AB4:** a robust finding contradicts a statement of the paper or demands
  a journal decision. **No submission decision without a new decision by the
  authors.**

## 7. Technical

- Scripts carry a reference to this preregistration, the date, the role, the
  sources, the package versions, the seeds and the note "no `dWT`/`iv` main
  metric".
- Archive sizes and SHA-256 are logged before the analysis.
- No main figure; any supplementary figure only through
  `figure_style/publication_style.R`.

## 8. The statement this session starts from

> We are not looking for a substitute for the missing replication and not for a
> new axis. We test whether a lesion response is visible within stably
> annotated populations in a clinical single-cell data set. If no data set
> passes the sample and annotation gates, that is a clean limitation result and
> no reason to rewrite the narrative.

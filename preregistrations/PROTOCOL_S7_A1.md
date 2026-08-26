> Translated from the German original of 2026-08-19. The content, the dates
> and every number are unchanged.

# S7 addendum 1 — the release of the new scRNA GEO screen

Dated **2026-08-19**, before the new search query.

S7-AB2 from `GSE196652` stands. The search is not an attempt to reinterpret
that outcome but a new, prospective test of whether an independent cohort with
real sample replicates is available. The search string and the complete
inclusion rules are given in `PREREG_S7.md` addendum 1.

Before the search query, no new GSE metadata, results or downloads were looked
at. The search may capture all hits up to `retmax=1000`; there is no
result-dependent selection of a single candidate.

The complete hit list and the machine metadata triage are stored in
`s7n1_geo_summaries.csv` and `s7n1_geo_triage.csv`; the candidates checked
directly and their screening decision are stored in the screening table and
added to this protocol. A computation begins only after a protocol of its own
for the screening step and a further addendum fixing the data preparation.

## The result of the screen before the context download

The direct GEO check confirms no further eligible 2 x 2 single-cell data set.
`GSE166824` is already excluded in S1 with `A4` (the differentiated wild type
is missing); `GSE249471` is bulk and contains only the induced day-7 state;
`GSE324998`, `GSE241505` and `GSE255646` have no matching control or
undifferentiated-differentiated structure; `GSE150768` is a mixed human and mouse disease
design without this 2 x 2.

`GSE337700` is the only new data set that remains useful as a separate context:
three clinical fracture non-union and three healed control samples, but no
undifferentiated arm and no perturbation. Following `PREREG_S7.md` addendum 2 it is
therefore downloaded in full but not taken into the reference metric. The
context analysis requires a preregistration of its own.

The data set was then downloaded as `GSE337700_RAW.tar` (1.0 GB), and the eight
RDS files it contains were verified with `tar -tf` and extracted. The RDS
contents were not read before the following preregistration.

## A procedural deviation in the metadata navigation

After the frozen main query was complete, two narrower E-Search filter queries
were run purely to navigate the metadata already stored (one filter with
`osteogenic/chondrogenic/day 0`, one with `mesenchymal stem cell` and
perturbation terms). That was a deviation from the rule of a single query and is
therefore disclosed here. The filter responses were not stored as a new result
pool, no expression value was computed, and no download other than `GSE337700`
was prompted by such a filter. `s7n1_geo_summaries.csv` from the complete
945-hit query remains the binding search pool.

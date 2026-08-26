"""
73_primary_references.py -- the primary studies in the reference apparatus

Purpose  The iScience submission checklist states this explicitly:

           "All datasets, program code, and methods used in your manuscript
            must be appropriately cited in the text and listed in the
            reference section, either in the form of the publications where
            they were first reported or in the form of independent persistent
            identifiers such as the DOI."

         Until 2026-08-24 the primary papers of the reanalysed series lived
         ONLY in Supplementary Table 12. That is not enough: every series
         used must appear in the reference list AND be cited in the text.

         This script generates both from TS12 -- the list block and the
         citation bracket in STAR Methods -- and writes them between markers
         in manuscript/MANUSCRIPT.md. The two therefore cannot drift apart,
         and code/71_check_references.py cross-checks them in both directions.

         Two series are genuinely unpublished (GEO: "Citation missing").
         They are cited via their accession -- exactly the second form the
         checklist allows -- rather than attributed to a paper that does not
         exist.

Input    figures/data/TS12_primary_publications.csv
Output   manuscript/MANUSCRIPT.md (between the markers, in place)
Runtime  seconds
"""
from __future__ import annotations

import os
import pathlib
import re
import sys

import pandas as pd

_env = os.environ.get("PAPER_V2_ROOT")
ROOT = (pathlib.Path(_env) if _env
        else pathlib.Path(__file__).resolve().parents[1])
MS = ROOT / "manuscript" / "MANUSCRIPT.md"
TS12 = ROOT / "figures" / "data" / "TS12_primary_publications.csv"

MARKER_LIST = ("<!-- PRIMARY_SOURCES:START -->", "<!-- PRIMARY_SOURCES:END -->")
MARKER_CITE = ("<!-- PRIMARY_CITATIONS:START -->",
               "<!-- PRIMARY_CITATIONS:END -->")

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:  # noqa: BLE001
    pass

# Publications already in the main reference list; they must not appear
# twice. Keyed by DOI.
ALREADY_LISTED = {
    "10.1038/s41418-022-01035-7": "Hao et al. 2022",   # reference 16
}


def einzeln(df: pd.DataFrame) -> pd.DataFrame:
    """One row per PUBLICATION, not per series.

    Several series share a paper (the two SERPINA3 arms, the two chromatin
    cohorts), and one series can have two candidates (GSE33896, GSE184087).
    The earliest entry per DOI wins -- that is the paper in which the data
    were first reported.
    """
    pub = df[df.doi.notna() & (df.doi.astype(str).str.len() > 4)].copy()
    pub["year"] = pd.to_numeric(pub.year, errors="coerce")
    pub = pub.sort_values(["doi", "year"]).drop_duplicates("doi", keep="first")
    # Which series hang on which paper?
    serien = (df[df.doi.isin(pub.doi)]
              .groupby("doi").accession.apply(lambda s: ", ".join(sorted(set(s)))))
    pub["serien"] = pub.doi.map(serien)
    pub = pub[~pub.doi.isin(ALREADY_LISTED)]
    return pub.sort_values(["year", "first_author"])


def schluessel(r: pd.Series) -> str:
    """The author-year key used to cite in the text -- surname only, uniform
    across all entries ("Kim et al. 2018"), never mixed with initials.

    GEO delivers the first author as "Surname Initials"; surnames can have
    several tokens ("La Manna F"). A token counts as initials when it is
    nothing but capital letters (optionally dotted); any token with lower-case
    letters is part of the surname.
    """
    if pd.isna(r.first_author):
        return f"Anon et al. {int(r.year)}"
    teile = str(r.first_author).split()
    nach = [t for t in teile
            if not re.fullmatch(r"[A-Z]{1,5}\.?", t.strip("."))]
    name = " ".join(nach) if nach else teile[0]
    return f"{name} et al. {int(r.year)}"


def sauber(v: object) -> str:
    """Volume and page numbers without the decimal that pandas attaches to a
    numeric column: from 133.0 to 133, from 2000051.0 to 2000051."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    t = str(v).strip()
    if t.lower() in ("nan", "none", ""):
        return ""
    teile = [x[:-2] if x.endswith(".0") and x[:-2].isdigit() else x
             for x in t.replace("\u2013", "-").split("-")]
    return "-".join(teile)


def eintrag(r: pd.Series) -> str:
    """One reference entry in the Cell Press pattern, without an annotation:
    the accession -> publication mapping lives in Supplementary Table 12."""
    autor = str(r.first_author).strip()
    titel = str(r.title).strip().rstrip(".")
    # A title that ends in "?"/"!" carries its own terminal punctuation;
    # appending a full stop would produce "genes?.".
    if not titel.endswith(("?", "!")):
        titel += "."
    teile = [f"{autor} *et al.*" if autor else "", titel]
    z = str(r.journal).strip()
    if z and z.lower() != "nan":
        band = sauber(r.volume)
        seiten = sauber(r.pages)
        band = f" **{band}**" if band else ""
        seiten = f", {seiten}" if seiten else ""
        teile.append(f"*{z}*{band}{seiten} ({int(r.year)}).")
    else:
        teile.append(f"({int(r.year)}).")
    teile.append(f"doi:{r.doi}")
    return " ".join(t for t in teile if t)


def main() -> int:
    d = pd.read_csv(TS12)
    pub = einzeln(d)
    unpub = d[d.status.astype(str).str.startswith("unpublished")]

    # ---- reference-list block. Final numbering is applied by the package
    # build (64_build_submission.py), following order of first occurrence;
    # here we only count consecutively.
    zeilen = []
    for i, (_, r) in enumerate(pub.iterrows(), 1):
        zeilen.append(f"{i}. {eintrag(r)}")
        zeilen.append("")
    if len(unpub):
        for _, r in unpub.iterrows():
            i += 1
            zeilen.append(f"{i}. {str(r.title).strip().rstrip('.')}. "
                          f"Gene Expression Omnibus, {r.accession}.")
            zeilen.append("")

    # ---- citation bracket for STAR Methods
    keys = [schluessel(r) for _, r in pub.iterrows()]
    akzessionen = ", ".join(unpub.accession)
    zitat = ("The primary publication of every reanalysed series is cited "
             "here [" + "; ".join(keys) + "] and mapped to its accession in "
             f"Supplementary Table 12; the {len(unpub)} series that carry no "
             f"publication ({akzessionen}) are cited by accession.")

    s = MS.read_text(encoding="utf-8")
    for marke, inhalt in ((MARKER_LIST, "\n".join(zeilen)),
                          (MARKER_CITE, zitat)):
        a, e = marke
        if a not in s or e not in s:
            print(f"  ! marker missing: {a} -- block not written")
            continue
        i, j = s.index(a) + len(a), s.index(e)
        s = s[:i] + "\n\n" + inhalt + "\n\n" + s[j:]
    MS.write_text(s, encoding="utf-8")

    print("73_primary_references.py")
    print(f"  {d.accession.nunique()} series -> {len(pub)} references of their own "
          f"({len(ALREADY_LISTED)} already in the list, "
          f"{len(unpub)} unpublished)")
    print(f"  citation bracket with {len(keys)} keys written")
    return 0


if __name__ == "__main__":
    sys.exit(main())

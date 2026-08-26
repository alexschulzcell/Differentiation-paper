# -*- coding: utf-8 -*-
"""
28_geo_primary_publications.py -- the primary publications of the GEO series used.

Every series used needs its primary publication with author, year, journal,
volume, pages and DOI. That can be pulled from GEO, and it is pulled here
rather than copied by hand.

The script queries NCBI E-utilities (`db=gds` for the series, `db=pubmed` for
the linked PMIDs). Series without a linked publication are marked as such --
they then have to be documented by hand or flagged as unpublished.

Output   data_raw/_referenz/geo_primaerarbeiten.csv (with the retrieval date)
Runtime  about a minute; needs network access
"""
from __future__ import annotations

import datetime as dt
import json
import os
import pathlib
import time
import urllib.parse
import urllib.request

import pandas as pd

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
AUS = WURZEL / "data_raw" / "_referenz" / "geo_primaerarbeiten.csv"
AUS.parent.mkdir(parents=True, exist_ok=True)

# The 18 perturbation data sets (Table S1), the seven patient cohorts (S4) and
# the orthogonal chromatin and methylome series.
SERIEN = [
    "GSE102732", "GSE137035", "GSE145235", "GSE190542", "GSE205432",
    "GSE218101", "GSE221128", "GSE226565", "GSE227512", "GSE245585",
    "GSE247528", "GSE251698",
    "GSE160207", "GSE186141", "GSE228522", "GSE22855", "GSE292600",
    "GSE58435", "GSE77758",
    "GSE151311", "GSE151315", "GSE332758", "GSE33896", "GSE129266",
    "GSE129031",
    # Added later: these two series stood as placeholders in the source files
    # (GSE_SERPINA3ch, GSE_MIR181) and were therefore missing here. Both are
    # perturbation data sets among the eighteen and need their primary
    # publication like all the others.
    "GSE247491", "GSE184087",
    # The postnatal growth plate: computed, reported as "not calibratable" --
    # and therefore to be documented exactly like the levels that carry something.
    "GSE288028",
]


def js(u: str):
    return json.load(urllib.request.urlopen(u, timeout=60))


def main() -> None:
    heute = dt.date.today().isoformat()
    treffer, alle_pmids = [], set()
    for g in SERIEN:
        try:
            r = js("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                   + urllib.parse.urlencode({"db": "gds",
                                             "term": "%s[Accession]" % g,
                                             "retmode": "json"}))
            ids = r["esearchresult"]["idlist"]
            if not ids:
                treffer.append((g, "", ""))
                continue
            s = js("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
                   "esummary.fcgi?"
                   + urllib.parse.urlencode({"db": "gds", "id": ids[0],
                                             "retmode": "json"})
                   )["result"][ids[0]]
            pm = [str(x) for x in (s.get("pubmedids") or [])]
            alle_pmids.update(pm)
            treffer.append((g, ";".join(pm), s.get("title", "")))
        except Exception as e:                       # noqa: BLE001
            treffer.append((g, "", "FEHLER: %s" % e))
        time.sleep(0.4)

    lit = {}
    if alle_pmids:
        s = js("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
               + urllib.parse.urlencode({"db": "pubmed",
                                         "id": ",".join(sorted(alle_pmids)),
                                         "retmode": "json"}))["result"]
        for p in sorted(alle_pmids):
            d = s.get(p, {})
            aid = {x.get("idtype"): x.get("value")
                   for x in d.get("articleids", [])}
            lit[p] = {
                "erstautor": (d.get("authors") or [{}])[0].get("name", ""),
                "jahr": (d.get("pubdate", "") or "")[:4],
                "zeitschrift": d.get("source", ""),
                "band": d.get("volume", ""), "seiten": d.get("pages", ""),
                "doi": aid.get("doi", ""), "titel": d.get("title", ""),
            }

    zeilen = []
    for g, pms, titel in treffer:
        if not pms:
            zeilen.append({"gse": g, "pmid": "", "status": "keine verknuepfte "
                           "Publikation in GEO", "geo_titel": titel,
                           "abrufdatum": heute})
            continue
        for p in pms.split(";"):
            zeilen.append(dict({"gse": g, "pmid": p, "status": "ok",
                                "geo_titel": titel, "abrufdatum": heute},
                               **lit.get(p, {})))
    T = pd.DataFrame(zeilen)
    T.to_csv(AUS, index=False, encoding="utf-8")
    ohne = sorted(T.loc[T.status != "ok", "gse"].unique())
    print("%d series, %d publication rows, %d without a publication: %s"
          % (T.gse.nunique(), int((T.status == "ok").sum()), len(ohne),
             ", ".join(ohne)))
    print("-> %s" % AUS)


if __name__ == "__main__":
    main()

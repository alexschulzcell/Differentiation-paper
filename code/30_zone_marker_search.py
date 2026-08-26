# -*- coding: utf-8 -*-
"""A documented PubMed search for zone markers of the human growth plate.

The canonical markers of the four zones are documented from the literature --
independently of our own data and independently of the primary publication of
the postnatal series.

The queries run through NCBI E-utilities (db=pubmed) and are written into the
search log with their hit counts and the screening verdict. The script decides
nothing; it only collects the hit lists for the screening record.
"""
import json
import os
import pathlib
import time

import requests

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
T = WURZEL.parent / "Referenzdaten" / "GSE288028_ChuSciTranslMed2026"
T.mkdir(parents=True, exist_ok=True)
AUS = T / "pubmed_suche_zonenmarker.json"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

ABFRAGEN = {
    # A: the textbook reviews of growth-plate biology
    "A1": 'Kronenberg HM[Author] AND ("Developmental regulation of the growth plate")',
    "A2": 'Mackie EJ[Author] AND "endochondral ossification"',
    "A3": 'Kozhemyakina E[Author] AND ("chondrocyte development" OR "pathway to bone")',
    "A4": '"epiphyseal growth plate"[ti] AND "review"[pt] AND (markers OR zones)',
    # B: single-cell work on the growth plate (mouse and human), published
    #    before 2025, so that the annotation stays independent of the primary
    #    publication of the series being annotated
    "B1": '"growth plate"[tiab] AND "single-cell"[tiab] AND (markers[tiab] OR '
          'zonal[tiab] OR zones[tiab]) AND "2020/01/01"[PDAT] : "2024/12/31"[PDAT]',
    "B2": '"growth plate"[ti] AND ("single cell"[tiab] OR single-cell[tiab])',
    # C: targeted marker questions per zone
    "C1": '"resting zone"[tiab] AND chondrocyte*[tiab] AND (marker*[tiab] OR '
          '"stem cell*"[tiab]) AND (mouse OR murine OR human)',
    "C2": '"prehypertrophic"[tiab] AND (IHH[tiab] OR PTHrP[tiab] OR Indian hedgehog[tiab])',
    "C3": '"hypertrophic"[tiab] AND chondrocyte*[tiab] AND ("type X collagen"[tiab] '
          'OR COL10A1[tiab] OR MMP13[tiab]) AND (marker*[tiab] OR zonal[tiab])',
    "C4": '"proliferative zone"[tiab] AND "growth plate"[tiab] AND '
          '(proliferat*[tiab] AND marker*[tiab])',
}


def main() -> None:
    ergebnis = {}
    for kuerzel, term in ABFRAGEN.items():
        r = requests.get(EUTILS + "esearch.fcgi",
                         params={"db": "pubmed", "term": term,
                                 "retmax": "30", "retmode": "json"},
                         timeout=60)
        d = r.json()["esearchresult"]
        ids = d["idlist"]
        zeilen = []
        if ids:
            time.sleep(0.5)
            r2 = requests.get(EUTILS + "esummary.fcgi",
                              params={"db": "pubmed", "id": ",".join(ids),
                                      "retmode": "json"},
                              timeout=60)
            res = r2.json()["result"]
            for u in res["uids"]:
                e = res[u]
                autoren = ", ".join(a["name"] for a in e.get("authors", [])[:3])
                zeilen.append({
                    "pmid": u,
                    "jahr": e.get("pubdate", "")[:4],
                    "titel": e.get("title", ""),
                    "zeitschrift": e.get("source", ""),
                    "erste_autoren": autoren,
                    "doi": next((x["value"] for x in e.get("articleids", [])
                                 if x["idtype"] == "doi"), ""),
                })
        ergebnis[kuerzel] = {"term": term, "anzahl": int(d.get("count", len(ids))),
                             "treffer": zeilen}
        print(kuerzel, "->", d.get("count"), "hits,", len(zeilen), "loaded")
        time.sleep(0.5)

    with open(AUS, "w", encoding="utf-8") as f:
        json.dump(ergebnis, f, ensure_ascii=False, indent=1)
    print("-> %s" % AUS)


if __name__ == "__main__":
    main()

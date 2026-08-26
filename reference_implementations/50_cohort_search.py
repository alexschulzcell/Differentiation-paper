# -*- coding: utf-8 -*-
"""
50_cohort_search.py -- systematic cohort search for phases M-B and M-C.

Runs the search axes defined in `preregistrations/PRAEREG_M_B.md` §2
against the GEO DataSets database (NCBI E-utilities) and records for each
search the date, query string and hit count. NOTHING is downloaded -- the
script only collects accession numbers and metadata for the screening
table.

Output:
  derived_data/M_patienten/suchlauf.csv     one row per search
  derived_data/M_patienten/treffer_roh.csv  one row per series found
"""
from __future__ import annotations

import pathlib
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import ERGEBNISSE  # noqa: E402

AUS = ERGEBNISSE / "M_patienten"
AUS.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEUTE = date.today().isoformat()
RETMAX = 1000

# Human series only; superseries and pure profile collections drop out at
# screening, not here.
FILTER = '"Homo sapiens"[Organism] AND "gse"[Filter]'

SUCHACHSEN = {
    "OI": '(("osteogenesis imperfecta") OR (COL1A1) OR (COL1A2))',
    "MPS": '(("mucopolysaccharidosis") OR ("mucopolysaccharidoses") OR (ARSB) OR (IDUA) OR (IDS))',
    "PSACH_MED": '(("pseudoachondroplasia") OR ("multiple epiphyseal dysplasia") OR (COMP) OR (MATN3))',
    "FGFR3": '(("achondroplasia") OR ("hypochondroplasia") OR ("thanatophoric") OR (FGFR3))',
    "SHOX": '(("SHOX") OR ("Leri-Weill") OR ("dyschondrosteosis"))',
    "CCD_CMPD": '(("cleidocranial dysplasia") OR ("campomelic dysplasia") OR (RUNX2) OR (SOX9))',
    "FOP": '(("fibrodysplasia ossificans progressiva") OR (ACVR1))',
    "FREITEXT_DYSPLASIE": '("skeletal dysplasia" OR "chondrodysplasia" OR "osteochondrodysplasia")',
    "FREITEXT_KLEINWUCHS": '("short stature" OR "growth failure" OR "dwarfism")',
    "FREITEXT_WACHSTUMSFUGE": '("growth plate" OR "epiphyseal" OR "physis cartilage")',
}


def hole(url: str, versuche: int = 3) -> str:
    for i in range(versuche):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:            # noqa: BLE001
            if i == versuche - 1:
                raise
            print("   ... retrying after error:", e)
            time.sleep(3)
    return ""


def esearch(term: str) -> list[str]:
    url = (f"{EUTILS}/esearch.fcgi?db=gds&retmax={RETMAX}&retmode=json"
           f"&term={urllib.parse.quote(term)}")
    txt = hole(url)
    return re.findall(r'"(\d{6,})"', txt)


def esummary(ids: list[str]) -> list[dict]:
    aus = []
    for i in range(0, len(ids), 100):
        block = ids[i:i + 100]
        url = (f"{EUTILS}/esummary.fcgi?db=gds&retmode=json"
               f"&id={','.join(block)}")
        import json
        d = json.loads(hole(url)).get("result", {})
        for uid in d.get("uids", []):
            r = d[uid]
            aus.append({
                "gse": r.get("accession", ""),
                "titel": (r.get("title", "") or "").replace("\n", " "),
                "zusammenfassung": (r.get("summary", "") or "").replace("\n", " "),
                "n_proben": r.get("n_samples", ""),
                "typ": r.get("gdstype", ""),
                "taxon": r.get("taxon", ""),
                "datum": r.get("PDAT", ""),
            })
        time.sleep(0.4)
    return aus


def main() -> None:
    print("=" * 78)
    print("Cohort search M-B / M-C  --  search date", HEUTE)
    print("=" * 78)

    laeufe, treffer = [], {}
    for name, kern in SUCHACHSEN.items():
        term = f"({kern}) AND {FILTER}"
        ids = esearch(term)
        print(f"{name:24s} {len(ids):5d} hits")
        laeufe.append({"achse": name, "datum": HEUTE, "datenbank": "GEO gds",
                       "suchstring": term, "n_treffer": len(ids),
                       "retmax": RETMAX})
        for r in esummary(ids):
            if not r["gse"].startswith("GSE"):
                continue
            r = dict(r)
            r["achse"] = name
            if r["gse"] in treffer:
                treffer[r["gse"]]["achse"] += "|" + name
            else:
                treffer[r["gse"]] = r
        time.sleep(0.4)

    L = pd.DataFrame(laeufe)
    T = pd.DataFrame(treffer.values()).sort_values("gse")
    L.to_csv(AUS / "suchlauf.csv", index=False)
    T.to_csv(AUS / "treffer_roh.csv", index=False)
    print(f"\n{len(L)} searches, {len(T)} unique series -> {AUS}")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
20_donor_search.py -- the search of phase M-D, by DESIGN rather than by
entity.

Preregistration: `preregistrations/PRAEREG_M_D.md` §5, dated 2026-08-22
before the first download. The search axes of phase M-B were directed at
entities (OI, MPS, FOP, ...). That was exactly the gap: what must be
searched for is the DESIGN -- a patient lesion with isogenic control and a
co-sequenced naive arm.

NOTHING is downloaded. The script collects accession numbers and metadata
for the screening table; the GSM metadata is fetched by `54b_meta.py`.

Output:
  derived_data/M_donoren/suchlauf.csv     one row per search
  derived_data/M_donoren/treffer_roh.csv  one row per series found
"""
from __future__ import annotations

import json
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

AUS = ERGEBNISSE / "M_donoren"
AUS.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEUTE = date.today().isoformat()
RETMAX = 1000

FILTER = '"Homo sapiens"[Organism] AND "gse"[Filter]'

# The three design axes of the preregistration, each crossed with the
# differentiation on which the module is defined.
DIFF = ('("osteogenic" OR "osteogenesis" OR "osteoblast" OR "chondrogenic" '
        'OR "chondrogenesis" OR "chondrocyte" OR "skeletal")')

SUCHACHSEN = {
    "ANLAGE_NAIV": ('("patient-derived" OR "patient derived" OR '
                    '"patient-specific") AND ("day 0" OR "undifferentiated" '
                    'OR "baseline" OR "iPSC stage")'),
    "ANLAGE_ISOGEN": '("isogenic control" OR "isogenic pair" OR "isogenic line") AND ("differentiation" OR "differentiated")',
    "ANLAGE_KORRIGIERT": '("gene-corrected" OR "gene corrected" OR "genetically corrected" OR "corrected iPSC") AND (iPSC OR "iPS cell" OR hiPSC)',
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
    return re.findall(r'"(\d{6,})"', hole(url))


def esummary(ids: list[str]) -> list[dict]:
    aus = []
    for i in range(0, len(ids), 100):
        url = (f"{EUTILS}/esummary.fcgi?db=gds&retmode=json"
               f"&id={','.join(ids[i:i + 100])}")
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
    print("Cohort search M-D  --  by design.  Search date", HEUTE)
    print("=" * 78)

    laeufe, treffer = [], {}
    for name, kern in SUCHACHSEN.items():
        term = f"({kern}) AND {DIFF} AND {FILTER}"
        ids = esearch(term)
        print(f"{name:20s} {len(ids):5d} hits")
        laeufe.append({"achse": name, "datum": HEUTE, "datenbank": "GEO gds",
                       "suchstring": term, "n_treffer": len(ids),
                       "retmax": RETMAX,
                       "an_obergrenze": len(ids) >= RETMAX})
        for r in esummary(ids):
            if not r["gse"].startswith("GSE"):
                continue
            if r["gse"] in treffer:
                treffer[r["gse"]]["achse"] += "|" + name
            else:
                r["achse"] = name
                treffer[r["gse"]] = r
        time.sleep(0.4)

    L = pd.DataFrame(laeufe)
    T = pd.DataFrame(treffer.values()).sort_values("gse")
    L.to_csv(AUS / "suchlauf.csv", index=False)
    T.to_csv(AUS / "treffer_roh.csv", index=False)
    print(f"\n{len(L)} searches, {len(T)} unique series -> {AUS}")


if __name__ == "__main__":
    main()

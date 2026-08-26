# -*- coding: utf-8 -*-
"""Fetch the abstracts of the core sources (for documenting the markers)."""
import os
import pathlib
import time

import requests

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
T = WURZEL.parent / "Referenzdaten" / "GSE288028_ChuSciTranslMed2026"
T.mkdir(parents=True, exist_ok=True)

PMIDS = {
    # Reviews / textbook knowledge
    "12748651": "Kronenberg2003",
    "17659995": "Mackie2008",
    "21642379": "Mackie2011",
    "25715393": "Kozhemyakina2015",
    "31795305": "Hallett2019_review",
    # Zone-specific
    "30401834": "Mizuhashi2018_RZ",
    "40025030": "Otsuru2025_APOE_RZ",
    "34309509": "Hallett2021_WntInhib_RZ",
    "34346115": "Renthal2021_roundcelllayer",
    "8895385": "Lee1996_PTHrP",
    "12050144": "Kobayashi2002_PTHrPIhh",
}

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
aus = []
for pmid, name in PMIDS.items():
    r = requests.get(EUTILS + "efetch.fcgi",
                     params={"db": "pubmed", "id": pmid, "rettype": "abstract",
                             "retmode": "text"}, timeout=60)
    aus.append("##### %s | PMID %s\n%s\n" % (name, pmid, r.text))
    print(name, "ok")
    time.sleep(0.6)

with open(T / "abstracts_kernquellen.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(aus))
print("-> %s" % (T / "abstracts_kernquellen.txt"))

# -*- coding: utf-8 -*-
"""Fetch the open-access full texts of the core sources and search them for zone markers."""
import os
import pathlib
import re
import time

import requests

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
T = WURZEL.parent / "Referenzdaten" / "GSE288028_ChuSciTranslMed2026"
T.mkdir(parents=True, exist_ok=True)

PMCS = {
    "PMC6251707": "Mizuhashi2018_RZ",
    "PMC11873292": "KodamaOtsuru2025_APOE_RZ",
    "PMC4352987": "Kozhemyakina2015",
    "PMC6929081": "Hallett2019_review",
}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

for pmc, name in PMCS.items():
    r = requests.get(EUTILS + "efetch.fcgi", params={"db": "pmc", "id": pmc,
                                                     "retmode": "xml"},
                     timeout=120)
    txt = re.sub(r"<[^>]+>", " ", r.text)
    txt = re.sub(r"\s+", " ", txt)
    p = T / ("%s_%s.xml.txt" % (name, pmc))
    p.write_text(txt, encoding="utf-8")
    print(name, "->", len(txt), "characters")
    time.sleep(1)

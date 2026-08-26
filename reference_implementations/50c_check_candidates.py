# -*- coding: utf-8 -*-
"""
50c_check_candidates.py -- sample metadata of the manual-screening
candidates.

For each accession number passed in, fetches the GSM level from GEO (SOFT,
`targ=gsm&form=text&view=brief`) and stores title, source tissue and
characteristics per sample. Only then is it decidable whether a series
fulfils E1-E5 from `PRAEREG_M_B.md` §2 -- the series title never suffices.

No expression matrix is loaded; metadata only.

Output: data_raw/_meta/<GSE>_proben.csv
"""
from __future__ import annotations

import pathlib
import re
import sys
import time
import urllib.request

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import DATEN, ERGEBNISSE  # noqa: E402

META = DATEN / "_meta"
META.mkdir(parents=True, exist_ok=True)
AUS = ERGEBNISSE / "M_patienten"

URL = ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
       "?acc={gse}&targ=gsm&form=text&view=brief")


def soft(gse: str) -> str:
    p = META / f"{gse}.soft"
    if p.exists() and p.stat().st_size > 200:
        return p.read_text("utf-8", "replace")
    with urllib.request.urlopen(URL.format(gse=gse), timeout=180) as r:
        t = r.read().decode("utf-8", "replace")
    p.write_text(t, encoding="utf-8")
    time.sleep(1.0)
    return t


def parse(gse: str) -> pd.DataFrame:
    zeilen, akt = [], None
    for ln in soft(gse).splitlines():
        if ln.startswith("^SAMPLE"):
            if akt:
                zeilen.append(akt)
            akt = {"gsm": ln.split("=")[1].strip(), "merkmale": []}
        elif akt is None:
            continue
        elif ln.startswith("!Sample_title"):
            akt["titel"] = ln.split("=", 1)[1].strip()
        elif ln.startswith("!Sample_source_name"):
            akt["quelle"] = ln.split("=", 1)[1].strip()
        elif ln.startswith("!Sample_characteristics"):
            akt["merkmale"].append(ln.split("=", 1)[1].strip())
        elif ln.startswith("!Sample_platform_id"):
            akt["plattform"] = ln.split("=", 1)[1].strip()
        elif ln.startswith("!Sample_library_strategy"):
            akt["strategie"] = ln.split("=", 1)[1].strip()
    if akt:
        zeilen.append(akt)
    for z in zeilen:
        z["merkmale"] = " | ".join(z.pop("merkmale"))
    d = pd.DataFrame(zeilen)
    d.insert(0, "gse", gse)
    return d


def main(gses: list[str]) -> None:
    for gse in gses:
        try:
            d = parse(gse)
        except Exception as e:            # noqa: BLE001
            print(f"{gse:11s} ERROR {e}")
            continue
        d.to_csv(META / f"{gse}_proben.csv", index=False)
        pf = d.plattform.nunique() if "plattform" in d else 0
        print(f"\n=== {gse}  {len(d)} samples, {pf} platform(s) ===")
        for _, r in d.head(60).iterrows():
            print("  %-11s %-46s %s" % (r.gsm, str(r.get("titel"))[:46],
                                        str(r.get("merkmale"))[:110]))


if __name__ == "__main__":
    main(sys.argv[1:])

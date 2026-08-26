# -*- coding: utf-8 -*-
"""
_glossary.py -- write column_glossary.csv

The panel data under figures/data/ and every delivered table carry English
column names. The stored intermediates under derived_data/ and results/ keep
the short internal names that the analysis scripts read and write, because
renaming them would mean renaming them in scripts that cannot be re-run
without the raw data. This script writes the bridge between the two, from the
same single mapping that code/_display.py applies when the panel data are
written -- so the glossary cannot drift from the translation.

Columns that do not appear in the glossary are either already English or are
identifiers (sample names, gene symbols, accessions).

Output  column_glossary.csv
"""
from __future__ import annotations

import ast
import collections
import csv
import os
import pathlib
import re
import sys

_env = os.environ.get("PAPER_V2_ROOT")
ROOT = (pathlib.Path(_env) if _env
        else pathlib.Path(__file__).resolve().parents[1])


def mapping() -> dict[str, str]:
    src = (ROOT / "code" / "_display.py").read_text(encoding="utf-8")
    block = re.search(r"^SPALTEN: dict\[str, str \| None\] = \{(.*?)^\}",
                      src, re.S | re.M).group(1)
    return {k: v for k, v in ast.literal_eval("{" + block + "}").items() if v}


def main() -> int:
    m = mapping()
    seen: collections.Counter[str] = collections.Counter()
    for d in ("derived_data", "results"):
        for p in (ROOT / d).rglob("*.csv"):
            if "_ueberholt" in p.parts:
                continue
            try:
                with p.open(encoding="utf-8", newline="") as f:
                    header = next(csv.reader(f))
            except Exception:  # noqa: BLE001
                continue
            for c in header:
                seen[c] += 1
    rows = [(k, m[k], seen[k]) for k in sorted(seen) if k in m]
    out = ROOT / "column_glossary.csv"
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["internal_name", "english_name", "files"])
        w.writerows(rows)
    print(f"column_glossary.csv -- {len(rows)} of {len(seen)} column names "
          f"in derived_data/ and results/")
    return 0


if __name__ == "__main__":
    sys.exit(main())

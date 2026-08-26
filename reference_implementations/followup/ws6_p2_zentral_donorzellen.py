# -*- coding: utf-8 -*-
"""
ws6_p2_zentral_donorzellen.py -- WS6, the central P2 check.

Question: does the fixed 173-gene module also carry in exactly those cells
whose OWN lineage-marker calibration did NOT pass?

WHY not the 18 study-dWT vectors (ws6_p1p2_liniensunabhaengigkeit.py):
that file showed that concordance there lies at 0.73-1.00 in ALL 18
datasets -- but that is **circular**: `ri` is the majority vote of EXACTLY
these 18 dWT vectors (S5_konvergente_gene.csv <-
20f_konvergente_dWT.csv <- median sign over the 18 points, PAPER/
reference_implementations/paper_daten.py lines ~80-92). Testing a dataset against its own
contribution to the majority is no test. The project already named this
trap itself in `reference_implementations/54d_circularity.py` (MODULBILDEND set) and
defused it at donor level there: "only NOT module-forming" (SERPINA3
series) gives S1 z +4.51, while "only module-forming" gives S1 z +0.88
(n.s.).

THIS SCRIPT goes one step further than 54d: it computes the module
concordance on ALL 14 phase-D cells, not only the 7 whose own lineage
calibration passed (`derived_data/M_donoren/zellen.pkl`, field
`meta.eichung_bestanden` -- the same calibration as in `eichung.csv`).
So far, failed cells were excluded from EVERY module analysis (rule: no
number without a passed positive control). Here the question is exactly
reversed: does the module DEVIATE in them?

Statistic: `_module.konkordanz` (directed sign test against `ri`,
baseline-stratified null) -- the same implementation as everywhere in the
project, NOT reinvented.
"""
from __future__ import annotations

import pathlib
import pickle
import sys

import numpy as np
import pandas as pd

HIER = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parents[1] / "reference_implementations"))
from _module import MODUL, konkordanz  # noqa: E402

WURZEL = HIER.parents[1]
AUS = WURZEL / "derived_data" / "followup"
AUS.mkdir(parents=True, exist_ok=True)

RI = MODUL.set_index("symbol")["ri"]
MODULBILDEND = {"GSE218101", "GSE221128", "GSE245585", "LAMA5_USC"}

with open(WURZEL / "derived_data" / "M_donoren" / "zellen.pkl", "rb") as f:
    DATEN = pickle.load(f)

zeilen = []
for studie, d in DATEN.items():
    dwt, basis, meta = d["dwt"], d["basis"], d["meta"]
    for _, row in meta.iterrows():
        z = row["zelle"]
        delta = dwt[z].dropna()
        r = konkordanz(delta.reindex(RI.index).dropna(), RI,
                        hintergrund=delta, schichtung=basis,
                        nziehungen=4000)
        r.update(studie=studie, zelle=z, achse=row.get("achse", ""),
                  eichung_bestanden=bool(row["eichung_bestanden"]),
                  modulbildende_studie=studie in MODULBILDEND)
        zeilen.append(r)
        print("  %-9s %-24s calib=%-5s moduleform=%-5s  concordance %.3f  "
              "z %+6.2f  p %8.4g  n=%s"
              % (studie, z, r["eichung_bestanden"], r["modulbildende_studie"],
                 r.get("konkordanz", np.nan), r.get("konkordanz_z", np.nan),
                 r.get("konkordanz_p", np.nan), r.get("n", "-")))
        pd.DataFrame(zeilen).to_csv(
            AUS / "ws6_p2_modul_je_zelle_alle14.csv", index=False)

T = pd.DataFrame(zeilen)
T = T[T.status == "ok"].copy()

print("\n" + "=" * 78)
print("breakdown by eichung_bestanden (ALL 14 cells)")
print("=" * 78)
for b in [True, False]:
    s = T[T.eichung_bestanden == b]
    print("  eichung_bestanden=%-5s  n=%2d  concordance-z median %.2f  "
          "z>0 & p<0.05: %d/%d"
          % (b, len(s), s.konkordanz_z.median(),
             int(((s.konkordanz_z > 0) & (s.konkordanz_p < 0.05)).sum()), len(s)))

print("\nthe same ONLY for the NOT module-forming studies (SERPINA3 series,")
print("the only truly independent cells):")
U = T[~T.modulbildende_studie]
for b in [True, False]:
    s = U[U.eichung_bestanden == b]
    if len(s) == 0:
        print("  eichung_bestanden=%-5s  n=0 -- no cell in this stratum" % b)
        continue
    print("  eichung_bestanden=%-5s  n=%2d  concordance-z median %.2f  "
          "z>0 & p<0.05: %d/%d  (cells: %s)"
          % (b, len(s), s.konkordanz_z.median(),
             int(((s.konkordanz_z > 0) & (s.konkordanz_p < 0.05)).sum()),
             len(s), ", ".join(s.zelle)))

T.to_csv(AUS / "ws6_p2_modul_je_zelle_alle14.csv", index=False)
print("\n-> %s" % (AUS / "ws6_p2_modul_je_zelle_alle14.csv"))

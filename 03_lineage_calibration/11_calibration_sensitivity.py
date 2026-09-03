# -*- coding: utf-8 -*-
"""
11_calibration_sensitivity.py -- carries the calibration balance of the 18
datasets, or drops it.

Two objections against `10_calibration_18_datasets.py`, both noticed when
checking the coverage and both computed here rather than talked away:

  (E1) **Not calibratable is not failed.** Two datasets (points 5 and 13)
       have only ONE measurable marker of their axis. `kontrast` returns
       "zu wenige Gene" there. Counting them as failures would be wrong --
       the same distinction phase M-B already made with `M6`.

  (E2) **Part of the markers is unreachable in vitro.** `OSTEOGEN` contains
       `SOST`, `DMP1`, `PHEX`, `MEPE` and `PTH1R` -- osteocyte and in-vivo
       markers that a 14- to 21-day MSC or iPSC culture does not reach.
       They pull the set toward zero, independent of whether the experiment
       differentiated osteogenically.

**The honesty status of this file: POST HOC, not preregistered.** It stands
here on the third line because it must also stand there in the report.
`_marker.py` is NOT changed -- the sensitivity analysis uses a subset and
names it as such. The preregistered version remains primary; this one only
says whether the headline count depends on the marker choice.

Output: derived_data/M_kalibrierung/eichung_empfindlichkeit.csv
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import ERGEBNISSE, kontrast, lade_dwt_je_punkt, gencode_karte  # noqa: E402

WURZEL = pathlib.Path(__file__).resolve().parents[1]
AUS = ERGEBNISSE / "M_kalibrierung"
# The per-gene dWT values come from the frozen matrix in this repository
# (see 00_shared/_module.py::lade_dwt_je_punkt). The author's session tree is
# only a fallback and is not part of this repository; point at it with
# SCHERENPAPER_SITZUNGEN if you have it.
GENE20D = (pathlib.Path(os.environ["SCHERENPAPER_SITZUNGEN"])
           / "20_Exploration" / "derived_data"
           if os.environ.get("SCHERENPAPER_SITZUNGEN") else None)

# Maturation markers NOT reachable in vitro within 14-21 days. Rationale per
# gene: osteocyte program (SOST, DMP1, PHEX, MEPE) or receptor of mature bone
# (PTH1R). The list was fixed BEFORE looking at the results of this file and
# is not readjusted.
UNERREICHBAR = ["SOST", "DMP1", "PHEX", "MEPE", "PTH1R"]

OSTEOGEN_INVITRO = [g for g in OSTEOGEN if g not in UNERREICHBAR]
MIN_MARKER = 3          # `kontrast` requires >= 3 per side

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def saetze(satz_a_name: str, invitro: bool) -> tuple[list, list]:
    M = {"OSTEOGEN": OSTEOGEN_INVITRO if invitro else OSTEOGEN,
         "CHONDROGEN": CHONDROGEN, "ADIPOGEN": ADIPOGEN,
         "MYOGEN": MYOGEN, "NAIV": NAIV}
    a = M[satz_a_name]
    b = [g for n, s in M.items() if n != satz_a_name for g in s]
    return a, b


def urteil(r: dict, dwt: pd.Series, a: list, b: list) -> str:
    n_a = len([g for g in a if g in dwt.index])
    n_b = len([g for g in b if g in dwt.index])
    if n_a < MIN_MARKER or n_b < MIN_MARKER:
        return "nicht eichbar"
    if r.get("status") != "ok":
        return "nicht eichbar"
    return ("bestanden" if (r["p"] < 0.05 and r["kontrast"] > 0)
            else "durchgefallen")


def main() -> None:
    log("=" * 78)
    log("Sensitivity of the calibration balance -- POST HOC, not "
        "preregistered")
    log("=" * 78)
    log("Omitted in the in-vitro version: %s" % ", ".join(UNERREICHBAR))
    log("")

    karte = gencode_karte()
    KO = pd.read_csv(WURZEL / "derived_data" / "manuscript" / "f1_kohorte.csv")
    arm = dict(zip(KO.punkt, KO.arm))
    name = dict(zip(KO.punkt, KO.datensatz))

    zeilen = []
    for p, G in lade_dwt_je_punkt(GENE20D):
        G = G.copy()
        G["symbol"] = [karte.get(str(g).split(".")[0]) for g in G.gen]
        G = G[G.symbol.notna() & G.dWT.notna()]
        dwt = G.groupby("symbol").dWT.median()
        satz_a = "OSTEOGEN" if arm.get(p) == "osteogen" else "CHONDROGEN"

        z = dict(punkt=p, datensatz=name.get(p, ""), arm=arm.get(p, ""),
                 satz_a=satz_a)
        for kurz, invitro in (("vorreg", False), ("invitro", True)):
            a, b = saetze(satz_a, invitro)
            r = kontrast(dwt, a, b)
            z["n_%s_a" % kurz] = len([g for g in a if g in dwt.index])
            z["urteil_%s" % kurz] = urteil(r, dwt, a, b)
            z["kontrast_%s" % kurz] = r.get("kontrast", np.nan)
            z["z_%s" % kurz] = r.get("z", np.nan)
            z["p_%s" % kurz] = r.get("p", np.nan)
        zeilen.append(z)
        log("  %2d %-28s %-11s | vorreg n%2d %-14s z %+6.2f | invitro n%2d "
            "%-14s z %+6.2f"
            % (p, str(name.get(p, ""))[:28], z["arm"],
               z["n_vorreg_a"], z["urteil_vorreg"], z["z_vorreg"],
               z["n_invitro_a"], z["urteil_invitro"], z["z_invitro"]))

    T = pd.DataFrame(zeilen)
    T.to_csv(AUS / "eichung_empfindlichkeit.csv", index=False)

    log("")
    log("-" * 78)
    for kurz, titel in (("vorreg", "preregistered marker sets"),
                        ("invitro", "only in-vitro reachable markers")):
        u = T["urteil_%s" % kurz]
        eichbar = T[u != "nicht eichbar"]
        n_best = int((eichbar["urteil_%s" % kurz] == "bestanden").sum())
        log("%s:" % titel)
        log("   calibratable %2d of 18 | passed %2d | failed %2d (%.0f %% "
            "of the calibratable)" % (len(eichbar), n_best, len(eichbar) - n_best,
                                      100 * (len(eichbar) - n_best) / len(eichbar)))
        for a in ("chondrogen", "osteogen"):
            s = eichbar[eichbar.arm == a]
            if len(s):
                nb = int((s["urteil_%s" % kurz] == "bestanden").sum())
                log("      %-11s %d of %d passed" % (a, nb, len(s)))
    log("-" * 78)
    (AUS / "56b_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()

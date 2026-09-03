# -*- coding: utf-8 -*-
"""
paper_daten_donoren.py -- carries the results of phase M-D into derived_data/manuscript
and builds the two map tables for Figure 5.

Nothing is computed here; only selection and renaming and -- for the maps --
assembling one row per level from the numbers already reported in the
protocols. Every row of a map carries the file it comes from.

Complements `paper_daten.py` and `paper_daten_medizin.py`; existing files
remain untouched.
"""
from __future__ import annotations

import pathlib

import pandas as pd

W = pathlib.Path(__file__).resolve().parents[1]
E = W / "derived_data"
D = W / "derived_data" / "manuscript"
T = W / "derived_data" / "reference_tables"


def w(d: pd.DataFrame, name: str, ziel: pathlib.Path = D) -> None:
    d.to_csv(ziel / name, index=False)
    print("  %-38s %4d Zeilen" % (name, len(d)))


# ---------------------------------------------------------------------------
# The levels map (Fig. 5).  One row per measured level, with
#   * calibration: passed / withdrawn / not calibratable
#   * detection limit in the unit of the respective level
#   * observed value for programme and lesion response
#   * status
# The numbers come from the protocols; the column `quelle` says which one.
# None of this is recomputed here.
# ---------------------------------------------------------------------------
LANDKARTE = [
    dict(ebene="perturbation models (18)", kurz="models",
         einheit="genes vs expected", eichung="passed",
         mde="0.35 z", programm=173 / 7.9, laesion=7 / 8.0,
         prog_ueber=True, laes_ueber=False, status="exploratory",
         quelle="f6_internal_summary.csv"),
    dict(ebene="external triangulation (11)", kurz="external",
         einheit="synthesis z", eichung="per cell",
         mde="no single threshold reached", programm=1.21, laesion=0.0,
         prog_ueber=False, laes_ueber=False, status="follow-up",
         quelle="f6_s12_fixed173_summary.csv"),
    dict(ebene="donors within one experiment (7 cells)", kurz="donors",
         einheit="S1 (mean rho)", eichung="dWT itself, 7/14 passed",
         mde="0.344", programm=0.349, laesion=float("nan"),
         prog_ueber=True, laes_ueber=None, status="exploratory",
         quelle="derived_data/M_donoren/statistik.csv"),
    dict(ebene="chromatin (H3K27ac, ATAC)", kurz="chromatin",
         einheit="concordance z", eichung="lineage markers",
         mde="0.588", programm=4.88, laesion=float("nan"),
         prog_ueber=True, laes_ueber=None, status="exploratory",
         quelle="f7_ebenen.csv"),
    dict(ebene="promoter methylome (450K)", kurz="methylome",
         einheit="concordance z", eichung="lineage markers (z +2.94)",
         mde="not reached", programm=1.84, laesion=float("nan"),
         prog_ueber=False, laes_ueber=None, status="exploratory",
         quelle="f7_ebenen.csv"),
    dict(ebene="patient vs control (7 cohorts)", kurz="patients",
         einheit="concordance C", eichung="withdrawn (NAIV unsuitable)",
         mde="0.576-0.689", programm=0.0, laesion=0.0,
         prog_ueber=False, laes_ueber=False, status="exploratory",
         quelle="m3_patienten.csv"),
    dict(ebene="human genetics (7 panels)", kurz="genetics",
         einheit="odds ratio", eichung="OR 17-52 / OR 3.7-5.8",
         mde="OR 1.59", programm=1.00, laesion=1.00,
         prog_ueber=False, laes_ueber=False, status="confirmatory",
         quelle="m4_anker.csv"),
]

# The detection-limit map (Fig. 5B): observed value and limit in ONE common
# unit -- the multiple of the null SD, i.e. z. For levels whose statistic is
# already a z itself, that is the value; for the rest it is the z formed from
# observation, null mean and null SD. The threshold is everywhere the same:
# 2.8 null SDs (MDE80).
GRENZEN = [
    dict(kurz="models", ebene="perturbation models", z=None,
         hinweis="convergence 22x expected", erreichbar=True),
    dict(kurz="external", ebene="external triangulation", z=1.21,
         hinweis="synthesis p 0.0002, no single cell", erreichbar=True),
    dict(kurz="donors", ebene="donors, programme (S1)", z=3.00,
         hinweis="2 of 7 cells over own MDE80", erreichbar=True),
    dict(kurz="donors_iv", ebene="donors, lesion response", z=None,
         hinweis="1 calibrated cell - not measurable", erreichbar=False),
    dict(kurz="donors_syn", ebene="donors, study synthesis (S1)", z=1.66,
         hinweis="MDE80 1.15 on a scale whose maximum is 1.0",
         erreichbar=False),
    dict(kurz="chromatin", ebene="chromatin (H3K27ac)", z=4.88,
         hinweis="only value over MDE80 in the whole study",
         erreichbar=True),
    dict(kurz="methylome", ebene="promoter methylome", z=1.84,
         hinweis="same differentiation, other layer", erreichbar=True),
    dict(kurz="patients", ebene="patient vs control", z=0.0,
         hinweis="no cohort over its own MDE80", erreichbar=True),
    dict(kurz="genetics", ebene="human genetics", z=0.0,
         hinweis="OR 1.00 at a limit of OR 1.59", erreichbar=True),
]


def main() -> None:
    M = E / "M_donoren"

    print("=== Abbildung 3: die donoraufgeloeste Rechnung (Phase M-D) ===")
    w(pd.read_csv(M / "eichung.csv"), "m5_eichung.csv")
    w(pd.read_csv(M / "je_zelle.csv"), "m5_je_zelle.csv")
    w(pd.read_csv(M / "statistik.csv"), "m5_statistik.csv")
    w(pd.read_csv(M / "auslassung.csv"), "m5_auslassung.csv")
    w(pd.read_csv(M / "zirkularitaet.csv"), "m5_zirkularitaet.csv")
    w(pd.read_csv(M / "synthese.csv"), "m5_synthese.csv")
    w(pd.read_csv(M / "zellen_sichtung.csv"), "m5_zellen_sichtung.csv")
    w(pd.read_csv(M / "sichtung_hand.csv"), "m1_sichtung_anlage.csv")
    w(pd.read_csv(M / "selbsttest.csv"), "m5_selbsttest.csv")

    print("\n=== Abbildung 5: die beiden Landkarten ===")
    w(pd.DataFrame(LANDKARTE), "m6_ebenen_landkarte.csv")
    w(pd.DataFrame(GRENZEN), "m6_nachweisgrenzen.csv")

    print("\n=== Ergaenzungstabellen ===")
    w(pd.read_csv(M / "eichung.csv"), "S15_donorzellen_eichung.csv", T)
    w(pd.read_csv(M / "je_zelle.csv"), "S16_donorzellen_statistik.csv", T)
    w(pd.read_csv(M / "sichtung_hand.csv"), "S17_sichtung_anlage.csv", T)
    w(pd.DataFrame(LANDKARTE), "S18_ebenen_landkarte.csv", T)


if __name__ == "__main__":
    main()

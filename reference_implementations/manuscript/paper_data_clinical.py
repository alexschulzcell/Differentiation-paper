# -*- coding: utf-8 -*-
"""
paper_daten_medizin.py -- carries the results of phases M-A, M-B and M-C into
derived_data/manuscript so that the plotting scripts keep reading ONLY from
derived_data/manuscript.

Nothing is computed here, only selection and renaming.

Complements `paper_daten.py`; the existing f* files remain untouched --
in particular `f4_krankheitsanreicherung.csv`, the old negative finding.
"""
from __future__ import annotations

import pathlib

import pandas as pd

W = pathlib.Path(__file__).resolve().parents[2]
E = W / "derived_data"
D = W / "derived_data" / "manuscript"
T = W / "derived_data" / "reference_tables"


def w(d: pd.DataFrame, name: str, ziel: pathlib.Path = D) -> None:
    d.to_csv(ziel / name, index=False)
    print("  %-38s %4d Zeilen" % (name, len(d)))


def main() -> None:
    print("=== Abbildung 1 D/E: die beiden neuen Sichtungen ===")
    S = pd.read_csv(E / "M_patienten" / "sichtung_mechanisch.csv")
    b = (S.groupby(["code"]).size().rename("n").reset_index()
           .sort_values("n", ascending=False))
    b["ebene"] = "Patientenkohorten"
    w(b, "m1_sichtung_patienten.csv")

    C = pd.read_csv(E / "M_diagnosen" / "sichtung.csv")
    c = C.groupby("code").size().rename("n").reset_index().sort_values(
        "n", ascending=False)
    c["ebene"] = "Diagnosen 2x2"
    w(c, "m1_sichtung_diagnosen.csv")

    K = pd.read_csv(E / "M_patienten" / "kohorten_sichtung.csv")
    w(K, "m1_kohorten.csv")

    print("\n=== Abbildung 3: Patientenkonkordanz (Phase M-B) ===")
    R = pd.read_csv(E / "M_patienten" / "streuung.csv")
    sp = ["gse", "entitaet", "gewebe", "satz", "n_patienten", "n_kontrollen",
          "n", "U", "konkordanz", "konkordanz_null", "konkordanz_null_sd",
          "konkordanz_z", "konkordanz_p", "konkordanz_mde80",
          "konkordanz_ki_lo", "konkordanz_ki_hi", "rang_z", "rang_p",
          "kontrast_z", "kontrast_p", "eichung_bestanden"]
    w(R.reindex(columns=sp), "m3_patienten.csv")
    w(pd.read_csv(E / "M_patienten" / "eichung.csv"), "m3_eichung.csv")
    w(pd.read_csv(E / "M_patienten" / "synthese.csv"), "m3_synthese.csv")

    print("\n=== Abbildung 4: humangenetischer Anker (Phase M-A) ===")
    A = pd.read_csv(E / "M_humangenetik" / "anker.csv")
    w(A, "m4_anker.csv")
    w(pd.read_csv(E / "M_humangenetik" / "anker_power.csv"), "m4_power.csv")
    w(pd.read_csv(E / "M_humangenetik" / "eichung_A.csv"), "m4_eichung.csv")

    print("\n=== Ergaenzungstabellen ===")
    w(pd.read_csv(E / "M_patienten" / "kohorten_sichtung.csv"),
      "S10_patientenkohorten.csv", T)
    w(pd.read_csv(E / "M_diagnosen" / "sichtung.csv"),
      "S11_sichtung_diagnosen.csv", T)
    w(A, "S12_humangenetischer_anker.csv", T)
    w(R.reindex(columns=sp), "S13_patientenkonkordanz.csv", T)
    P = pd.read_csv(E / "M_humangenetik" / "panels.csv")
    w(P.groupby(["panel", "quelle", "abruf"]).size().rename("n_gene")
       .reset_index(), "S14_panelfassungen.csv", T)


if __name__ == "__main__":
    main()

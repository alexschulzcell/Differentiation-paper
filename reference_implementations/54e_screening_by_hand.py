# -*- coding: utf-8 -*-
"""Manual screening of the phase-M-D hits -- fixed as a CSV."""
import pandas as pd, pathlib
# --- path parameters (2026-08-23) ----------------------------------------
# Previously hard-coded to ".../Scherenpaper_Folgeprojekt/derived_data/M_donoren"
# -- this directory never existed there; the results live under
# `Paper v2/Ergebnisse`. Overridable via PAPER_V2_ROOT.
import os
_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
AUS = WURZEL / "derived_data" / "M_donoren"

URTEIL = {
 "GSE111163": ("M4", "Skelettmuskel, keine osteogene/chondrogene Achse"),
 "GSE113297": ("A2", "iPSC-abgeleitete kraniofaziale MSC, keine Laesion und keine Kontrollgruppe im 2x2-Sinn"),
 "GSE116928": ("M4", "Somitenuhr, Entwicklungsmodell ohne Laesion"),
 "GSE143453": ("M4", "FSHD2-Myotuben, Muskel; zudem Einzelkernformat (A9)"),
 "GSE143492": ("A9", "Einzelzell-/Einzelkernformat, Muskel"),
 "GSE161025": ("M4", "myogene Spezifizierung, keine skelettale Achse"),
 "GSE169724": ("M4", "FSHD-Muskelfasern"),
 "GSE183525": ("M4", "FOP-Monozyten, entzuendliche Signatur, keine Differenzierungsachse"),
 "GSE189053": ("M4", "DMD, myogene Kulturen"),
 "GSE214626": ("M4", "Myotone Dystrophie, iPSC, keine skelettale Achse"),
 "GSE218101": ("EIN", "MPS VI, 4 Patientenlinien mit isogener Korrektur, iPS -> Tag 14 chondrogen"),
 "GSE224181": ("M4", "Rhabdomyosarkom, Tumorentitaet"),
 "GSE224182": ("M4", "Rhabdomyosarkom, Tumorentitaet"),
 "GSE225148": ("M4", "LGMDR21, Skelettmuskel"),
 "GSE245585": ("EIN", "RB1-Patientenlinie gegen WT, MSC, Osteogenese Tag 0/7/14/21"),
 "GSE26272": ("M3", "ZNF145-Ueberexpression in gesunden MSC, kein Patientendefekt, kein 2x2"),
 "GSE302312": ("M3", "Implantate aus stabilem/hypertrophem Knorpel, keine Laesion"),
 "GSE307443": ("M4", "LAMA2, Muskelstammzellen"),
 "GSE36098": ("M4", "Mesoangioblasten, Muskel"),
 "GSE37521": ("M3", "Dermis- und Fettgewebs-Stromazellen, gesunde Spender, keine Laesion"),
 "GSE58123": ("M4", "familiaeres Krebssyndrom, keine skelettale Achse"),
 "GSE70955": ("M4", "DMD, Muskel"),
}
VORAB = {
 "GSE218101": ("EIN", "vorab benannt; 4 Zellen (Line #1-#4), chondrogen"),
 "GSE221128": ("EIN", "vorab benannt; FOP/resFOP = EINE Linie -> 1 Zelle, ex1-3 sind Replikate"),
 "GSE247491": ("EIN", "vorab benannt; 3 Spender, chondrogen; E2 nicht erfuellt (siRNA) -> nur dWT-tragend"),
 "GSE247528": ("EIN", "vorab benannt; 3 Spender, osteogen; E2 nicht erfuellt (siRNA) -> nur dWT-tragend"),
 "GSE244375": ("A1", "kein naiver Arm: Tag 34 gegen Tag 44, beide differenzierter Knorpel; zudem Klone einer Linie"),
 "GSE148728": ("A1", "kein naiver Arm: nur differenzierte Knorpelpellets; 'biologische Replikate' sind keine Spender"),
 "LAMA5_USC": ("EIN", "eigene Daten; WT1-3/KO9/46/75 sind Klone EINER Linie -> 1 Spender, 2 Achsen"),
}

zeilen = []
for gse, (code, grund) in URTEIL.items():
    zeilen.append(dict(quelle="Suchlauf 2026-08-22", gse=gse,
                       urteil="EIN" if code == "EIN" else "AUS",
                       code="-" if code == "EIN" else code, begruendung=grund))
for gse, (code, grund) in VORAB.items():
    zeilen.append(dict(quelle="vorab benannt (PRAEREG_M_D §5)", gse=gse,
                       urteil="EIN" if code == "EIN" else "AUS",
                       code="-" if code == "EIN" else code, begruendung=grund))
T = pd.DataFrame(zeilen).drop_duplicates(subset=["gse", "quelle"])
T.to_csv(AUS / "sichtung_hand.csv", index=False)
print(T.to_string())
print("\nIN (EIN):", (T.urteil == "EIN").sum(), "| OUT (AUS):", (T.urteil == "AUS").sum())

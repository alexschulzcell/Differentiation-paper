# -*- coding: utf-8 -*-
"""
21_diagnosis_shear.py -- phase M-C: the shear on real diagnoses.

Phase C asks the same as the main part, but the lesion is a **diagnosis**
rather than an engineering intervention. For this it needs the same 2x2
design as the eighteen points:

        naive/undifferentiated   x   differentiated
        control                  x   diagnosis

Figure 1 predicts that most candidates fail at **`A1`** -- no naive arm.
**This failure is itself the result** (`PLAN` §4) and is reported, not
omitted.

All 127 series of the manual screening from `11_cohort_prescreen.py` are audited,
at their GSM metadata, with the same exclusion codes as
`derived_data/reference_tables/S1_sichtung_alle_datensaetze.csv`.

Output: derived_data/M_diagnosen/sichtung.csv
"""
from __future__ import annotations

import pathlib
import re
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import DATEN, ERGEBNISSE  # noqa: E402

import importlib.util  # noqa: E402
spec = importlib.util.spec_from_file_location(
    "p50c", WURZEL / "data_acquisition" / "12_cohort_sample_metadata.py")
p50c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p50c)

AUS = ERGEBNISSE / "M_diagnosen"
AUS.mkdir(parents=True, exist_ok=True)
QUELLE = ERGEBNISSE / "M_patienten"

# A naive arm is recognizable by an undifferentiated starting state.
NAIV = re.compile(
    r"\bday ?0\b|\bd0\b|\bt ?= ?0\b|\b0 ?h\b|undifferentiated|"
    r"\bunstimulated\b|\buninduced\b|before differentiation|baseline|"
    r"\bmonolayer\b|expansion medium|growth medium|basal medi|"
    r"\bhiPSC\b|\biPSC\b(?!.*derived)|\bhESC\b|\bMSC\b|"
    r"mesenchymal (stem|stromal)", re.I)

# A differentiated arm.
DIFF = re.compile(
    r"day ?(?:[1-9]|[1-9]\d)\b|\bd(?:[1-9]|[1-9]\d)\b|week ?[1-9]|"
    r"differentiat|chondrogen|osteogen|adipogen|hypertroph|induced|"
    r"pellet|micromass|mineraliz", re.I)

# A diagnosis as lesion (not an engineering intervention).
DIAGNOSE = re.compile(
    r"disease state|diagnosis|\bpatient\b|genotype|mutation|mutant genes|"
    r"karyotype|affected|proband", re.I)

# Engineering intervention instead of a diagnosis.
EINGRIFF = re.compile(
    r"siRNA|shRNA|knock-?down|knock-?out|CRISPR|overexpress|transfect|"
    r"empty vector|scrambl|control siRNA", re.I)


def audit(gse: str) -> dict:
    try:
        d = p50c.parse(gse)
    except Exception as e:                          # noqa: BLE001
        return dict(gse=gse, urteil="AUS", code="A3",
                    begruendung=f"Metadaten nicht lesbar: {e}")
    if len(d) == 0:
        return dict(gse=gse, urteil="AUS", code="A3", begruendung="keine Proben")
    d.to_csv(DATEN / "_meta" / f"{gse}_proben.csv", index=False)
    txt = (d.get("titel", pd.Series(dtype=str)).astype(str) + " | " +
           d.get("quelle", pd.Series(dtype=str)).astype(str) + " | " +
           d.get("merkmale", pd.Series(dtype=str)).astype(str))
    ganz = " || ".join(txt)

    hat_naiv = bool(txt.str.contains(NAIV).any())
    hat_diff = bool(txt.str.contains(DIFF).any())
    hat_diag = bool(DIAGNOSE.search(ganz))
    hat_eing = bool(EINGRIFF.search(ganz))

    z = dict(gse=gse, n_proben=len(d), naiver_arm=hat_naiv,
             diff_arm=hat_diff, diagnose=hat_diag, eingriff=hat_eing)
    if not hat_diag:
        z.update(urteil="AUS", code="A2",
                 begruendung="keine Diagnose als Laesionsachse")
    elif not (hat_naiv and hat_diff):
        z.update(urteil="AUS", code="A1",
                 begruendung=("kein naiver Arm" if not hat_naiv
                              else "kein differenzierter Arm"))
    else:
        z.update(urteil="PRUEFEN", code="-",
                 begruendung="2x2 aus den Metadaten moeglich -- Handpruefung")
    return z


# Hand verdicts for the seven series that survive the mechanical audit.
# Each was checked against its GSM metadata (`data_raw/_meta/`), none from
# memory.
HANDURTEIL = {
    "GSE218101": ("BEREITS", "-",
                  "MPS VI, Tag 0/14, 4 Patientenlinien gegen isogene Korrektur "
                  "-- vollstaendiges 2x2, aber BEREITS einer der achtzehn "
                  "Punkte (ARSB / MPS VI)"),
    "GSE221128": ("BEREITS", "-",
                  "FOP, iMSC Tag 0/6, FOP gegen resFOP -- vollstaendiges 2x2, "
                  "aber BEREITS einer der achtzehn Punkte (ACVR1 / FOP)"),
    "GSE206213": ("AUS", "A2",
                  "Werner-Syndrom als WRN-Knockdown in hESC modelliert -- "
                  "Engineering-Eingriff, keine Diagnose; Linie statt Patient"),
    "GSE206214": ("AUS", "A2",
                  "dieselbe Studie, ChIP-Arm; WRN-Knockdown in hESC. Der "
                  "Tag-0/14-Arm des Begleitdatensatzes ist MC3T3-E1"),
    "GSE241507": ("AUS", "A7+M2",
                  "humaner OPLL-Arm ohne Differenzierungsachse; die "
                  "Tag-0/14-Achse liegt in MC3T3-E1 -- Maus, immortalisiert"),
    "GSE188698": ("AUS", "M4",
                  "Retinaorganoide, Mueller-Glia -- keine skelettale Entitaet"),
    "GSE222109": ("AUS", "A2",
                  "iPSC gegen Sklerotom, nur Wildtyp -- keine Diagnose als "
                  "Laesionsachse"),
}


def main() -> None:
    K = pd.read_csv(QUELLE / "kandidaten.csv")
    print("=" * 78)
    print("Phase M-C -- A1 audit of the %d manual-screening candidates" % len(K))
    print("=" * 78)
    zeilen = []
    for i, gse in enumerate(K.gse, 1):
        z = audit(gse)
        zeilen.append(z)
        print("%3d/%d  %-11s %-8s %-3s %s"
              % (i, len(K), gse, z["urteil"], z["code"], z["begruendung"][:60]))
    S = pd.DataFrame(zeilen)
    print("\n--- Manual screening of the %d survivors ---"
          % int((S.urteil == "PRUEFEN").sum()))
    for gse, (u, c, b) in HANDURTEIL.items():
        m = S.gse == gse
        if not m.any():
            continue
        S.loc[m, ["urteil", "code", "begruendung"]] = [u, c, b]
        print("   %-11s %-8s %-6s %s" % (gse, u, c, b[:58]))
    offen = S.gse[S.urteil == "PRUEFEN"].tolist()
    if offen:
        print("   NOT HAND-CHECKED: %s" % ", ".join(offen))

    S = S.merge(K[["gse", "titel", "n_proben", "achse"]].rename(
        columns={"n_proben": "n_proben_gds"}), on="gse", how="left")
    S.to_csv(AUS / "sichtung.csv", index=False)

    # New points for the statistic: the 2x2-capable ones minus those already
    # included. Reported separately, never mixed into the confirmatory
    # cohort.
    neu = S[S.urteil == "EIN"]
    neu.to_csv(AUS / "punkte.csv", index=False)

    print("\n" + "=" * 78)
    print(S.code.value_counts().to_string())
    mit_diagnose = S[S.code != "A2"]
    n1 = int((S.code == "A1").sum())
    print("\nOf %d candidates with a diagnosis as lesion axis, %d fail at "
          "A1 -- no naive arm: %.0f %%"
          % (len(mit_diagnose), n1, 100 * n1 / max(1, len(mit_diagnose))))
    print("Complete 2x2: %d series (%s) -- of these already in the cohort: %d"
          % (int((S.urteil.isin(["EIN", "BEREITS"])).sum()),
             ", ".join(S.gse[S.urteil.isin(["EIN", "BEREITS"])]),
             int((S.urteil == "BEREITS").sum())))
    print("NEW points for the statistic: %d" % len(neu))
    print("=" * 78)


if __name__ == "__main__":
    main()

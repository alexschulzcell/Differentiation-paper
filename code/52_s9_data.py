# -*- coding: utf-8 -*-
"""
52_s9_data.py -- panel data of the second cohort: panel F2F and figure S9.

Purpose  Reshapes the outputs of `code/23_lineage_contrast.py` into four panel
         CSV files. S9 carries the result of the second attempt at the lineage
         contrast, and the replication of the decoupling in an independent
         chromatin cohort that turned up while doing it.

         The decision rule was fixed before the run (see the header of
         23_lineage_contrast.py) and sent the lineage contrast to the
         supplement: calibration L passes in 0 of 4 windows in both new
         cohorts, as it already did in GSE332758.

Inputs   results/linienkontrast_{eichung,modultest,diagnose}.csv
         derived_data/B_atac/B_atac_eichung_je_achse.csv  (first cohort)
Outputs  figures/data/S9A_calibration_L_three_cohorts.csv
         figures/data/S9B_marker_sets_lineage_axis.csv
         figures/data/S9C_decoupling_second_cohort.csv
         figures/data/S9D_marker_sets_single_axes.csv
Runtime  seconds
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import _display  # noqa: E402  presentation layer: everything shown is English

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
RES = WURZEL / "results"
ERG = WURZEL / "derived_data"
AUS = WURZEL / "figures" / "data"
AUS.mkdir(parents=True, exist_ok=True)

KOHORTE = {"GSE332758": "GSE332758 (ATAC)",
           "GSE151311": "GSE151311 (ATAC)",
           "GSE151315": "GSE151315 (H3K27ac)"}


def schreib(df: pd.DataFrame, name: str, was: str) -> None:
    # Display values and column headers go out in English.
    df = _display.englisch(df)
    df.to_csv(AUS / (name + ".csv"), index=False)
    print("  %-40s %4d rows x %2d cols  %s"
          % (name, len(df), df.shape[1], was))


def main() -> None:
    print("52_s9_data.py -- panels of supplementary figure S9")

    e_neu = pd.read_csv(RES / "linienkontrast_eichung.csv")
    m_neu = pd.read_csv(RES / "linienkontrast_modultest.csv")
    d_neu = pd.read_csv(RES / "linienkontrast_diagnose.csv")

    # --- S9A  calibration L in all THREE cohorts. The first cohort comes from
    # the existing result file, with the same columns.
    alt = pd.read_csv(ERG / "B_atac" / "B_atac_eichung_je_achse.csv")
    alt = alt[alt.eichung == "L"].copy()
    alt["kohorte"] = "GSE332758"
    sp = ["kohorte", "fenster", "n_a", "n_b", "kontrast", "null_mittel",
          "null_sd", "z", "p", "mde80", "bestanden"]
    neu = e_neu[e_neu.eichung == "L"].copy()
    a = pd.concat([alt[sp], neu[sp]], ignore_index=True)
    a["kohorte_lang"] = a.kohorte.map(KOHORTE)
    a["assay"] = np.where(a.kohorte == "GSE151315", "H3K27ac", "ATAC")
    schreib(a, "S9A_calibration_L_three_cohorts",
            "passed %d of %d" % (int(a.bestanden.sum()), len(a)))

    # --- S9B  marker diagnosis on the difference axis, ALL THREE cohorts.
    # The first cohort is recomputed here with the same calculation, so that
    # the figure does not mix two cohorts from a script with one from the text.
    sys.path.insert(0, str(WURZEL / "reference_implementations"))
    from _marker import ADIPOGEN, NAIV, OSTEOGEN   # noqa: E402
    EPS = 0.05
    alt_diag = []
    for fen in ["P", "T10", "T50", "GB"]:
        M = pd.read_csv(ERG / "B_atac" / ("B_atac_matrix_%s.csv" % fen),
                        index_col=0)
        naiv = M["MSC-0d"]
        lr = pd.DataFrame({c: np.log2((M[c] + EPS) / (naiv + EPS))
                           for c in M.columns if c != "MSC-0d"})
        diffax = (lr[["OB-3d", "OB-5d", "OB-7d"]].mean(axis=1) -
                  lr[["AD-3d", "AD-5d", "AD-7d"]].mean(axis=1)).dropna()
        rng = np.random.default_rng(20260823)
        arr = diffax.values
        for name, satz in (("OSTEO", OSTEOGEN), ("ADIPO", ADIPOGEN),
                           ("NAIV", NAIV)):
            g = [x for x in satz if x in diffax.index]
            null = np.array([arr[rng.choice(len(arr), len(g),
                                            replace=False)].mean()
                             for _ in range(4000)])
            alt_diag.append(dict(kohorte="GSE332758", fenster=fen,
                                 achse="differenz", satz=name, n=len(g),
                                 mittel=float(diffax[g].mean()),
                                 z=float((diffax[g].mean() - null.mean()) /
                                         null.std(ddof=1))))
    b = pd.concat([pd.DataFrame(alt_diag),
                   d_neu[d_neu.achse == "differenz"]], ignore_index=True)
    b["kohorte_lang"] = b.kohorte.map(KOHORTE)
    schreib(b, "S9B_marker_sets_lineage_axis",
            "where each marker set sits on the lineage axis")

    # --- S9C  the decoupling in the second cohort: calibration against module
    ed = e_neu[(e_neu.kohorte == "GSE151315") & (e_neu.eichung == "D")]
    md = m_neu[(m_neu.kohorte == "GSE151315") &
               (m_neu.null == "H1 basisgeschichtet") &
               (m_neu.achse != "differenz")]
    c = md.merge(ed[["fenster", "achse", "z", "kontrast", "mde80",
                     "bestanden"]],
                 on=["fenster", "achse"], how="left",
                 suffixes=("", "_eichung"))
    c = c.rename(columns={"z": "eichung_z", "kontrast": "eichung_kontrast",
                          "mde80": "eichung_mde80",
                          "bestanden": "eichung_bestanden"})
    c = c[["fenster", "achse", "n", "konkordanz", "konkordanz_mde80",
           "konkordanz_z", "konkordanz_p", "ueber_grenze",
           "eichung_kontrast", "eichung_mde80", "eichung_z",
           "eichung_bestanden"]]
    schreib(c, "S9C_decoupling_second_cohort",
            "module above limit %d/%d; calibration passed %d/%d"
            % (int(c.ueber_grenze.sum()), len(c),
               int(c.eichung_bestanden.sum()), len(c)))

    # --- the decomposition in this cohort. Four quantities on one axis, one
    # point per window:
    #   undifferentiated -- is the undifferentiated state left?
    #   osteogenic       -- is the osteogenic lineage reached?
    #   adipogenic       -- does the culture run into the other lineage?
    #   module           -- does the matrix programme run regardless?
    # This decomposition is possible in no other cohort of the study, because
    # only here is there an undifferentiated arm AND two lineages.
    dm = d_neu[(d_neu.kohorte == "GSE151315") &
               (d_neu.achse.isin(["osteogen", "adipogen"]))].copy()
    dm = dm.rename(columns={"satz": "groesse"})[
        ["fenster", "achse", "groesse", "n", "z"]]
    mm = m_neu[(m_neu.kohorte == "GSE151315") &
               (m_neu.null == "H1 basisgeschichtet") &
               (m_neu.achse != "differenz")].copy()
    mm = mm.assign(groesse="Modul")[["fenster", "achse", "groesse", "n",
                                     "konkordanz_z"]]
    mm = mm.rename(columns={"konkordanz_z": "z"})
    f2f = pd.concat([dm, mm], ignore_index=True)
    f2f["achse_lang"] = f2f.achse.map({"osteogen": "osteogenic axis",
                                       "adipogen": "adipogenic axis"})
    f2f["art"] = np.where(f2f.groesse == "Modul", "module", "marker set")
    # This panel used to carry the main figure. Since the decomposition is
    # computed across all 18 data sets (code/25_decomposition_eighteen.py),
    # GSE151315 is the single case in the supplement, and the file name says
    # so.
    schreib(f2f, "S9C_second_cohort_decomposition",
            "GSE151315: naive state left, lineage not reached, module runs")

    # --- S9D  marker diagnosis on the single axes of the second cohort
    d = d_neu[(d_neu.kohorte == "GSE151315") &
              (d_neu.achse.isin(["osteogen", "adipogen"]))].copy()
    schreib(d, "S9D_marker_sets_single_axes",
            "naive state left, lineage not reached")


if __name__ == "__main__":
    main()

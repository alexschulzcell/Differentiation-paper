# -*- coding: utf-8 -*-
"""
10_decomposition_18_datasets.py -- the three-way decomposition on all 18 data sets.

Preregistered in `preregistrations/PRAEREG_F2F.md`, written BEFORE this run.
The decision rule stands there in §2 and once more here, so that it is in the
script header and not only in a second document:

    (i)   the UNDIFFERENTIATED markers fall: z <= -2 against their own null
    (ii)  the lineage markers of the data set's own axis do NOT rise above z = +2
    (iii) the module lies above its own detection limit

    The replication is CONFIRMED when (i) and (iii) hold and (ii) does not.
    When (ii) and (iii) hold -- the lineage is reached and the module runs --
    that is NOT a contradiction but the other case of decoupling, and it is
    reported as such. The decomposition is REFUTED when (i) holds and (iii)
    does not.

    Data sets that do not satisfy (i) are called "not decomposable" and carry
    NEITHER a positive NOR a negative finding (project rule 1).

    Across the cohort: replicated when, among the data sets satisfying (i),
    the majority satisfy (iii) and not (ii). No arithmetic mean over the 18 z
    values and no pooled test statistic -- the 18 are not independent (six
    study units of the donor level are among them). What is reported is a
    count with a Wilson interval.

NO SECOND IMPLEMENTATION. The contrast comes from
`00_shared/_module.py`, the marker sets from
`00_shared/_marker.py`, `dWT` from the frozen `20d_gene_*.csv`,
and the per-data-set module finding unchanged from
`derived_data/followup/ws6_p1p2_modul_je_datensatz.csv`. This script computes
the three marker contrasts and reshapes; nothing else.

THE GENE SPACE
--------------------------------------------------------------------------
The self-test of this script showed that the stored file
`eichung_achtzehn.csv` had been computed with a FILTERED marker set (only
genes with `dWT_n >= 15` of 18 data sets, the MIN_N filter of the internal
gene map). That filter preferentially removes the TERMINAL differentiation
markers and is wrong for a per-experiment positive control; the reasoning is
worked through in `03_lineage_calibration/12_calibration_gene_space.py`.
`03_lineage_calibration/10_calibration_18_datasets.py` is corrected
accordingly (2 of 18 rather than 3 of 18).

The decomposition therefore runs under both gene spaces:

  * `voll`      -- markers and background unfiltered, that is, the same gene
                   space as the corrected calibration. **Primary.** As a
                   self-test it is checked that the own-axis contrast of this
                   run reproduces the stored calibration file exactly.
  * `gefiltert` -- the old MIN_N gene space. **Sensitivity analysis**, so that
                   it stays visible that the verdict does not hang on the
                   filter.

NOTE a tautology that deserves naming: criterion (ii) -- "the lineage markers
of the own axis rise above z = +2" -- is THE SAME test statistic at THE SAME
threshold as the calibration. That the data sets satisfying (ii) are exactly
the data sets passing their calibration is therefore **true by definition and
not a finding**. The finding is (i) and (iii), which are independent of it,
and their joint behaviour. The legend and the text say so.

Inputs    20d_gene_*.csv, derived_data/M_kalibrierung/,
          derived_data/followup/ws6_p1p2_modul_je_datensatz.csv,
          derived_data/manuscript/f1_kohorte.csv
Outputs   results/zerlegung_achtzehn.csv, results/zerlegung_achtzehn_bilanz.csv,
          results/zerlegung_achtzehn_log.txt
Runtime   about two minutes (18 x 3 contrasts x 20 000 draws)
"""
from __future__ import annotations

import importlib.util
import os
import pathlib
import sys

import numpy as np
import pandas as pd

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
sys.path.insert(0, str(WURZEL / "00_shared"))

from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import lade_dwt_je_punkt  # noqa: E402
from _module import kontrast, wilson  # noqa: E402

# The per-gene dWT values come from the frozen matrix in this repository
# (see 00_shared/_module.py::lade_dwt_je_punkt). The author's session tree is
# only a fallback and is not part of this repository; point at it with
# SCHERENPAPER_SITZUNGEN if you have it.
GENE20D = (pathlib.Path(os.environ["SCHERENPAPER_SITZUNGEN"])
           / "20_Exploration" / "derived_data"
           if os.environ.get("SCHERENPAPER_SITZUNGEN") else None)

ERG = WURZEL / "derived_data"
NEU = WURZEL / "derived_data" / "followup"
RES = WURZEL / "results"
RES.mkdir(parents=True, exist_ok=True)

MARKERSAETZE = {"OSTEOGEN": OSTEOGEN, "ADIPOGEN": ADIPOGEN,
                "MYOGEN": MYOGEN, "CHONDROGEN": CHONDROGEN, "NAIV": NAIV}

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def karte_gencode() -> dict:
    """The full gene space: the same GENCODE map as in
    10_calibration_18_datasets.py -- imported, not rebuilt. The file name starts
    with a digit, hence importlib."""
    p = WURZEL / "03_lineage_calibration" / "10_calibration_18_datasets.py"
    spec = importlib.util.spec_from_file_location("_eichung18", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.symbolkarte()


def karte_intern() -> dict:
    """The filtered gene space (MIN_N >= 15 of 18) -- the map with which the
    calibration was mistakenly computed. Sensitivity analysis only."""
    g = pd.read_csv(WURZEL / "derived_data" / "R_intern" /
                    "R_interne_genkarte.csv")
    return dict(zip(g.ensembl, g.symbol))


def zerlege(karte: dict, kartenname: str, arm: dict, name: dict,
            eich: pd.DataFrame, mod: pd.DataFrame) -> pd.DataFrame:
    """The three marker contrasts per data set under one symbol map."""
    zeilen = []
    for pt, G in lade_dwt_je_punkt(GENE20D):
        G = G.copy()
        G["symbol"] = [karte.get(str(g).split(".")[0]) for g in G.gen]
        G = G[G.symbol.notna() & G.dWT.notna()]
        dwt = G.groupby("symbol").dWT.median()

        a = arm.get(pt, "")
        eigen = "OSTEOGEN" if a == "osteogen" else "CHONDROGEN"
        ander = "ADIPOGEN"          # the non-skeletal counter-axis, fixed
        z = {}
        for satz in ("NAIV", eigen, ander):
            b = [g for n, s in MARKERSAETZE.items() if n != satz for g in s]
            z[satz] = kontrast(dwt, MARKERSAETZE[satz], b)
        zeilen.append({
            "genkarte": kartenname,
            "punkt": pt, "datensatz": name.get(pt, ""), "arm": a,
            "n_gene_messbar": len(dwt),
            "eichung_z_abgelegt": float(eich.z.get(pt, np.nan)),
            "eichung_bestanden_abgelegt": bool(eich.bestanden.get(pt, False)),
            "eichung_status_abgelegt": str(eich.status.get(pt, "")),
            "satz_eigen": eigen, "satz_ander": ander,
            "naiv_kontrast": z["NAIV"].get("kontrast", np.nan),
            "naiv_z": z["NAIV"].get("z", np.nan),
            "naiv_mde80": z["NAIV"].get("mde80", np.nan),
            "naiv_status": z["NAIV"].get("status", ""),
            "eigen_kontrast": z[eigen].get("kontrast", np.nan),
            "eigen_z": z[eigen].get("z", np.nan),
            "eigen_status": z[eigen].get("status", ""),
            "ander_kontrast": z[ander].get("kontrast", np.nan),
            "ander_z": z[ander].get("z", np.nan),
            "ander_status": z[ander].get("status", ""),
            "modul_konkordanz": float(mod.konkordanz.get(pt, np.nan)),
            "modul_z": float(mod.konkordanz_z.get(pt, np.nan)),
            "modul_mde80": float(mod.konkordanz_mde80.get(pt, np.nan)),
            "modul_n": float(mod.n.get(pt, np.nan)),
        })
        log("  %2d  %-30s naive z %+6.2f | %-10s z %+6.2f | adipo z %+6.2f "
            "| module z %+6.2f"
            % (pt, str(name.get(pt, ""))[:30], zeilen[-1]["naiv_z"],
               eigen.lower(), zeilen[-1]["eigen_z"], zeilen[-1]["ander_z"],
               zeilen[-1]["modul_z"]))

    T = pd.DataFrame(zeilen)
    # ---- apply the preregistered rule --------------------------------------
    T["i_naiv_faellt"] = T.naiv_z <= -2
    T["ii_linie_erreicht"] = T.eigen_z > 2
    T["iii_module_ueber_grenze"] = T.modul_konkordanz > T.modul_mde80
    T["urteil"] = np.where(
        ~T.i_naiv_faellt, "nicht zerlegbar",
        np.where(T.iii_module_ueber_grenze,
                 np.where(T.ii_linie_erreicht,
                          "Entkopplung, anderer Fall (Linie erreicht)",
                          "Zerlegung bestaetigt"),
                 "Zerlegung widerlegt"))
    return T.sort_values("punkt")


def bilanz(T: pd.DataFrame, kartenname: str) -> dict:
    zerlegbar = T[T.i_naiv_faellt]
    n_z = len(zerlegbar)
    n_best = int((zerlegbar.urteil == "Zerlegung bestaetigt").sum())
    n_ander = int(zerlegbar.urteil.str.startswith("Entkopplung").sum())
    n_wider = int((zerlegbar.urteil == "Zerlegung widerlegt").sum())
    lo, hi = wilson(n_best, n_z) if n_z else (np.nan, np.nan)
    log("")
    log("-" * 78)
    log("GENE MAP: %s" % kartenname)
    log("(i) naive markers fall (z <= -2): %d of %d points -- the remaining %d"
        % (n_z, len(T), len(T) - n_z))
    log("    are called 'not decomposable' and carry no finding.")
    log("Among the %d decomposable points:" % n_z)
    log("    decomposition confirmed     %2d  (Wilson %.3f-%.3f)"
        % (n_best, lo, hi))
    log("    decoupling, the other case  %2d" % n_ander)
    log("    decomposition refuted       %2d" % n_wider)
    log("    -> %s"
        % ("REPLICATED" if n_z and n_best > n_z / 2 else "NOT replicated"))
    log("Module above its own limit, all 18: %d"
        % int(T.iii_module_ueber_grenze.sum()))
    log("Breakdown by stored calibration status:")
    for st, s in T.groupby(T.eichung_bestanden_abgelegt):
        log("    calibration %-13s n %2d | naive falls %2d | lineage reached %2d "
            "| module above limit %2d"
            % ("passed" if st else "failed", len(s),
               int(s.i_naiv_faellt.sum()), int(s.ii_linie_erreicht.sum()),
               int(s.iii_module_ueber_grenze.sum())))
    log("-" * 78)
    return {"genkarte": kartenname, "n_punkte": len(T), "n_zerlegbar_i": n_z,
            "n_bestaetigt": n_best, "n_anderer_fall": n_ander,
            "n_widerlegt": n_wider,
            "anteil_bestaetigt": (n_best / n_z) if n_z else np.nan,
            "wilson_lo": lo, "wilson_hi": hi,
            "n_module_ueber_grenze": int(T.iii_module_ueber_grenze.sum()),
            "repliziert": bool(n_z and n_best > n_z / 2)}


def main() -> None:
    log("=" * 78)
    log("The three-way decomposition on the 18 perturbation data sets")
    log("Preregistered: preregistrations/PRAEREG_F2F.md, 2026-08-24")
    log("=" * 78)

    KO = pd.read_csv(WURZEL / "derived_data" / "manuscript" / "f1_kohorte.csv")
    arm = dict(zip(KO.punkt, KO.arm))
    name = dict(zip(KO.punkt, KO.datensatz))
    eich = pd.read_csv(ERG / "M_kalibrierung" /
                       "eichung_achtzehn.csv").set_index("punkt")
    mod = pd.read_csv(NEU /
                      "ws6_p1p2_modul_je_datensatz.csv").set_index("punkt")

    teile, bilanzen = [], []
    for kn, kf in (("voll", karte_gencode), ("gefiltert", karte_intern)):
        log("")
        log("### Gene map: %s" % kn)
        T = zerlege(kf(), kn, arm, name, eich, mod)
        if kn == "voll":
            # Self-test: the own-axis contrast MUST reproduce the stored
            # calibration number exactly. If it does not, one of the two
            # computations is not what it is taken to be.
            d = (T.eigen_z - T.eichung_z_abgelegt).abs().max()
            log("  self-test against eichung_achtzehn.csv: max |Delta z| = %.2e"
                % d)
            assert not (d > 1e-9), ("the own-axis contrast misses the stored "
                                    "calibration: %.3g" % d)
        teile.append(T)
        bilanzen.append(bilanz(T, kn))

    A = pd.concat(teile, ignore_index=True)
    A.to_csv(RES / "zerlegung_achtzehn.csv", index=False)
    B = pd.DataFrame(bilanzen)
    B.to_csv(RES / "zerlegung_achtzehn_bilanz.csv", index=False)

    log("")
    log("Verdict identical under both maps: %s"
        % ("YES" if len(set(B.repliziert)) == 1 else "NO"))
    log("-> %s" % (RES / "zerlegung_achtzehn.csv"))
    (RES / "zerlegung_achtzehn_log.txt").write_text(
        chr(10).join(LOG) + chr(10), encoding="utf-8")


if __name__ == "__main__":
    main()

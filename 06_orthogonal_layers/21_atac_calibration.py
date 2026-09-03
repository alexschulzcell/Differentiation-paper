# -*- coding: utf-8 -*-
"""
21_atac_calibration.py -- layer B, final version.

Why this script replaces 23_
----------------------------
23_ calibrated the layer against the wrong question. What was calibrated was
LINEAGE SEPARATION (osteogenic vs adipogenic), but what was tested included
the DIFFERENTIATION AXIS (differentiated vs naive). A layer that resolves
differentiation cleanly but cannot keep the two lineages apart fails this
calibration although it serves the actual question. Each axis needs its own
calibration:

  CALIBRATION D (differentiation axis)  two-set contrast
        `lineage markers minus naive/proliferation markers`
        within each axis. Licenses the module tests on
        `osteogen` and `adipogen`.

  CALIBRATION L (lineage axis)  two-set contrast
        `osteogenic minus adipogenic markers` on the difference axis.
        Licenses only the module test on `differenz`.

Additionally two hardenings of the null, because the result on the
differentiation axis is positive and a positive finding must withstand more
than a negative one:

  H1  BASELINE ACCESSIBILITY. Module genes could simply be the genes that
      are closed in the naive state and hence can only open up. The null
      therefore draws decile-wise from the same MSC-0d accessibility.
  H2  REGION WIDTH. In the gene-body window the signal tracks gene length.
      The null draws decile-wise from the same window width.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
from _marker import ADIPOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import ERGEBNISSE, MODUL, konkordanz, kontrast, wilson  # noqa: E402

AUS = ERGEBNISSE / "B_atac"
LOG: list[str] = []
EPS = 0.05

FENSTER = {
    "P":   "promoter, TSS -2000/+500",
    "T10": "TSS +- 10 kb",
    "T50": "TSS +- 50 kb",
    "GB":  "gene body",
}


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def achsen(M: pd.DataFrame) -> pd.DataFrame:
    naiv = M["MSC-0d"]
    lr = pd.DataFrame({c: np.log2((M[c] + EPS) / (naiv + EPS))
                       for c in M.columns if c != "MSC-0d"})
    D = pd.DataFrame({
        "osteogen": lr[["OB-3d", "OB-5d", "OB-7d"]].mean(axis=1),
        "adipogen": lr[["AD-3d", "AD-5d", "AD-7d"]].mean(axis=1),
    })
    D["differenz"] = D.osteogen - D.adipogen
    D["basis"] = naiv
    return D


def main() -> None:
    log("=" * 78)
    log("Layer B -- chromatin, GSE332758: calibration per axis, then module test")
    log("=" * 78)

    daten = {f: achsen(pd.read_csv(AUS / ("B_atac_matrix_%s.csv" % f), index_col=0))
             for f in FENSTER}

    # ------------------------------------------------------ CALIBRATION D
    log("\n--- CALIBRATION D: does the layer resolve differentiation? --------")
    log("Contrast: lineage markers minus naive/proliferation markers, per axis.")
    log("%-5s %-9s %5s %5s %8s %8s %8s %7s %9s  %s"
        % ("win", "axis", "n_lin", "n_naiv", "contr.", "null", "MDE80",
           "z", "p", "verdict"))
    eich = []
    for f in FENSTER:
        for a, marker in (("osteogen", OSTEOGEN), ("adipogen", ADIPOGEN)):
            r = kontrast(daten[f][a], marker, NAIV)
            ok = r.get("status") == "ok" and r["kontrast"] >= r["mde80"]
            log("%-5s %-9s %5d %5d %+8.3f %+8.3f %+8.3f %+7.2f %9.4g  %s"
                % (f, a, r.get("n_a", 0), r.get("n_b", 0), r.get("kontrast", np.nan),
                   r.get("null_mittel", np.nan), r.get("mde80", np.nan),
                   r.get("z", np.nan), r.get("p", np.nan),
                   "PASSED" if ok else "failed"))
            eich.append(dict(eichung="D", fenster=f, achse=a, bestanden=ok,
                             **{k: v for k, v in r.items() if k != "status"}))

    # ------------------------------------------------------ CALIBRATION L
    log("\n--- CALIBRATION L: does the layer keep the two lineages apart? ----")
    log("Contrast: osteogenic minus adipogenic markers on the difference axis.")
    for f in FENSTER:
        r = kontrast(daten[f]["differenz"], OSTEOGEN, ADIPOGEN)
        ok = r.get("status") == "ok" and r["kontrast"] >= r["mde80"]
        log("%-5s %-9s %5d %5d %+8.3f %+8.3f %+8.3f %+7.2f %9.4g  %s"
            % (f, "differenz", r.get("n_a", 0), r.get("n_b", 0),
               r.get("kontrast", np.nan), r.get("null_mittel", np.nan),
               r.get("mde80", np.nan), r.get("z", np.nan), r.get("p", np.nan),
               "PASSED" if ok else "failed"))
        eich.append(dict(eichung="L", fenster=f, achse="differenz", bestanden=ok,
                         **{k: v for k, v in r.items() if k != "status"}))
    E = pd.DataFrame(eich)
    E.to_csv(AUS / "B_atac_eichung_je_achse.csv", index=False)

    # ------------------------------------------------------ MODULE TEST
    log("\n--- MODULE TEST ---------------------------------------------------")
    log("Expected sign +ri (open chromatin and transcription run in the")
    log("same direction). Null: background draw; H1 additionally decile-")
    log("stratified by baseline accessibility in the naive state.")
    sym_ri = dict(zip(MODUL.symbol, MODUL.ri))
    zeilen = []
    for f in FENSTER:
        D = daten[f]
        log("\nWindow %s (%s)" % (f, FENSTER[f]))
        for a in ("osteogen", "adipogen", "differenz"):
            lic = E[(E.fenster == f) & (E.achse == a)]
            bestanden = bool(lic.bestanden.iloc[0]) if len(lic) else False
            hg = D[a].dropna()
            mod = hg[hg.index.isin(sym_ri)]
            erw = pd.Series({s: sym_ri[s] for s in mod.index})
            for hname, sch in (("Hintergrund", None),
                               ("H1 basisgeschichtet", D["basis"])):
                r = konkordanz(mod, erw, hintergrund=hg, schichtung=sch)
                k = int((np.sign(mod.values) == erw.reindex(mod.index).values).sum())
                lo, hi = wilson(k, len(mod))
                schwelle = r["konkordanz"] >= r["konkordanz_mde80"]
                urteil = ("Ebene fuer diese Achse nicht geeicht" if not bestanden
                          else ("SCHWELLE ERREICHT" if schwelle else "unter Schwelle"))
                log("  %-10s %-20s n %3d | C %.3f [%.3f-%.3f] | Null %.3f+-%.3f "
                    "| z %+5.2f | p %.4g | %s"
                    % (a, hname, r["n"], r["konkordanz"], lo, hi,
                       r["konkordanz_null"], r["konkordanz_null_sd"],
                       r["konkordanz_z"], r["konkordanz_p"], urteil))
                zeilen.append(dict(ebene="ATAC", datensatz="GSE332758", fenster=f,
                                   achse=a, null=hname, geeicht=bestanden,
                                   schwelle_erreicht=bool(schwelle), k=k,
                                   wilson_lo=lo, wilson_hi=hi, urteil=urteil,
                                   **{x: y for x, y in r.items() if x != "status"}))

    Z = pd.DataFrame(zeilen)
    Z.to_csv(AUS / "B_atac_modultest_final.csv", index=False)

    # ------------------------------------------------- gene values for figure
    D = daten["T50"].copy()
    D["im_module"] = D.index.isin(sym_ri)
    D["ri"] = [sym_ri.get(s, np.nan) for s in D.index]
    D["markerklasse"] = ["osteogen" if s in OSTEOGEN else
                         "adipogen" if s in ADIPOGEN else
                         "naiv" if s in NAIV else "" for s in D.index]
    D.reset_index().rename(columns={"index": "symbol"}).to_csv(
        AUS / "B_atac_genwerte_T50.csv", index=False)

    # ------------------------------------------------------- summary
    log("\n--- Summary of layer B --------------------------------------------")
    g = Z[(Z.geeicht) & (Z.null == "H1 basisgeschichtet")]
    if len(g):
        for _, r in g.iterrows():
            log("  %-4s %-10s C %.3f (null %.3f) z %+.2f p %.4g -> %s"
                % (r.fenster, r.achse, r.konkordanz, r.konkordanz_null,
                   r.konkordanz_z, r.konkordanz_p, r.urteil))
    else:
        log("  No window passed a calibration.")

    (AUS / "B_atac_final_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

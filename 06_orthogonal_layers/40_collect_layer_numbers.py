# -*- coding: utf-8 -*-
"""
40_collect_layer_numbers.py -- collects the numbers of the orthogonal layers into
derived_data/manuscript so the figure script only reads from PAPER/.

The same separation as for Figures 1-4: computing and drawing are separate
runs, and the drawing script knows no raw data.

Written:
    f7_ebenen.csv        one row per orthogonal layer: design, number of
                         biological units, calibration passed, module test
    f7_robustheit.csv    the module finding across eight preparations
    f7_dwt_gegen_iv.csv  the paired comparison per layer
    f7_modul_gegen_marker.csv
    f7_streu_atac.csv    gene values of the replicated chromatin layer for
                         the scatter plot
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import ERGEBNISSE, MODUL, WURZEL  # noqa: E402

ZIEL = WURZEL / "derived_data" / "manuscript"
ZIEL.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


# Design of the layers -- maintained by hand, because it comes from the
# metadata of the datasets, not from a computation.
EBENEN = [
    dict(ebene="Chromatin", kurz="ATAC", datensatz="GSE224251",
         achse="osteogen", zelltyp="hMSC", einheiten=3,
         einheit="biologische Replikate", messung="ATAC-seq, Peakzaehlungen",
         schluessel="ATAC GSE224251 osteogen (n=3)"),
    dict(ebene="Chromatin", kurz="ATAC", datensatz="GSE332758",
         achse="osteogen", zelltyp="MSC-Linie", einheiten=1,
         einheit="biologische Linie", messung="ATAC-seq, BigWig",
         schluessel="ATAC GSE332758 osteogen (n=1)"),
    dict(ebene="Chromatin", kurz="ATAC", datensatz="GSE332758",
         achse="adipogen", zelltyp="MSC-Linie", einheiten=1,
         einheit="biologische Linie", messung="ATAC-seq, BigWig",
         schluessel="ATAC GSE332758 adipogen (n=1)"),
    dict(ebene="Chromatin", kurz="H3K27ac", datensatz="GSE129031",
         achse="chondrogen", zelltyp="BM-MSC", einheiten=2,
         einheit="Donorlinien", messung="H3K27ac ChIP-seq",
         schluessel="H3K27ac GSE129031 chondrogen (n=2)"),
    dict(ebene="Methylierung", kurz="DNAm 450K", datensatz="GSE129266",
         achse="chondrogen", zelltyp="BM-MSC", einheiten=2,
         einheit="gepaarte Donoren", messung="Illumina 450K",
         schluessel="DNAm GSE129266 chondrogen (n=2)"),
    dict(ebene="Methylierung", kurz="DNAm 27K", datensatz="GSE33896",
         achse="osteogen", zelltyp="hASC", einheiten=3,
         einheit="gepaarte Donoren", messung="Illumina 27K",
         schluessel="DNAm GSE33896 osteogen (n=3)"),
    dict(ebene="Methylierung", kurz="DNAm 27K", datensatz="GSE33896",
         achse="myogen", zelltyp="hASC", einheiten=3,
         einheit="gepaarte Donoren", messung="Illumina 27K",
         schluessel="DNAm GSE33896 myogen (n=3)"),
]


def main() -> None:
    log("=" * 78)
    log("Figure 5 -- collecting data (filename carries the old count 7)")
    log("=" * 78)

    Z = pd.read_csv(ERGEBNISSE / "Z_integration" / "Z_integration_statistik.csv")
    P = pd.read_csv(ERGEBNISSE / "Z_integration" / "Z_integration_partiell.csv")
    D = pd.read_csv(ERGEBNISSE / "Z_integration" / "Z_integration_dezile.csv")

    # ---- calibration per layer: ONE statistic for all six layers.
    # Calibration uses the same directed statistic that also measures the
    # module, applied to the canonical lineage marker set of the respective
    # axis (from 31_). A layer that does not find its own textbook marker
    # set can say nothing about the module.
    MM = pd.read_csv(ERGEBNISSE / "Z_integration" / "Z_modul_gegen_marker.csv")
    mk = MM[MM.satzart == "Marker"].set_index("ebene")
    md = MM[MM.satzart == "Modul"].set_index("ebene")
    gk = MM[MM.satzart == "Groessenkontrolle"].set_index("ebene")
    E = pd.DataFrame(EBENEN)
    E["marker_z"] = E.schluessel.map(mk.z)
    E["marker_p"] = E.schluessel.map(mk.p)
    E["marker_n"] = E.schluessel.map(mk.n)
    E["modul_z"] = E.schluessel.map(md.z)
    E["modul_p"] = E.schluessel.map(md.p)
    E["modul_n"] = E.schluessel.map(md.n)
    E["teilmengen_schlagen_marker"] = E.schluessel.map(1 - gk.p)

    eich = {}
    p = ERGEBNISSE / "B_atac" / "B2_GSE224251_eichung_hart.csv"
    if p.exists():
        t = pd.read_csv(p)
        r = t[t.schluessel == "nur Promotorpeaks | Peakfilter | CPM+Quantil"]
        if len(r):
            eich["ATAC GSE224251 osteogen (n=3)"] = (
                float(r.eichung2_z.iloc[0]), bool(r.eichung2_ok.iloc[0]))
    p = ERGEBNISSE / "B_atac" / "B_atac_eichung_je_achse.csv"
    if p.exists():
        t = pd.read_csv(p)
        for a, s in (("osteogen", "ATAC GSE332758 osteogen (n=1)"),
                     ("adipogen", "ATAC GSE332758 adipogen (n=1)")):
            r = t[(t.eichung == "D") & (t.fenster == "T50") & (t.achse == a)]
            if len(r):
                eich[s] = (float(r.z.iloc[0]), bool(r.bestanden.iloc[0]))
    p = ERGEBNISSE / "A_dnam" / "A_dnam450_GSE129266_eichung.csv"
    if p.exists():
        t = pd.read_csv(p)
        # expected direction is negative -> flip the sign for display
        eich["DNAm GSE129266 chondrogen (n=2)"] = (
            -float(t.z.iloc[0]), bool(t.bestanden.iloc[0]))

    E["eichung_z"] = E.schluessel.map(lambda s: eich.get(s, (np.nan, None))[0])
    E["eichung_ok"] = E.schluessel.map(lambda s: eich.get(s, (np.nan, None))[1])

    for dosis, marke in (("dWT_kons", "dwt"), ("iv_kons", "iv")):
        sub = Z[Z.dosis == dosis].set_index("ebene")
        E["z_" + marke] = E.schluessel.map(sub.z)
        E["p_" + marke] = E.schluessel.map(sub.p)
        E["zb_" + marke] = E.schluessel.map(sub.z_basisgeschichtet)
        E["pb_" + marke] = E.schluessel.map(sub.p_basisgeschichtet)
        sp = P[P.dosis == dosis].set_index("ebene")
        E["rho_" + marke] = E.schluessel.map(sp.rho_roh)
        E["rhop_" + marke] = E.schluessel.map(sp.rho_partiell)

    E.to_csv(ZIEL / "f7_ebenen.csv", index=False)
    log("\nf7_ebenen.csv -- %d layers" % len(E))
    log("%-34s %5s %8s %9s %9s %9s"
        % ("layer", "n", "calib z", "zb(dWT)", "zb(iv)", "rho_p(dWT)"))
    for _, r in E.iterrows():
        log("%-34s %5d %8s %+9.2f %+9.2f %+10.4f"
            % (r.schluessel, r.einheiten,
               ("%+.2f" % r.eichung_z) if pd.notna(r.eichung_z) else "-",
               r.zb_dwt, r.zb_iv, r.rhop_dwt))

    # f7_dosiskurve.csv belongs to 32_dose_response.py (directed and
    # baseline-controlled) and is explicitly NOT written here.
    log("\nDecile curve from 30_: %d rows, %d layers (stays in derived_data/)"
        % (len(D), D.ebene.nunique()))

    E[["schluessel", "ebene", "kurz", "datensatz", "achse", "einheiten",
       "eichung_z", "eichung_ok", "z_dwt", "p_dwt", "zb_dwt", "pb_dwt",
       "z_iv", "p_iv", "zb_iv", "pb_iv"]].to_csv(
        ZIEL / "f7_dwt_gegen_iv.csv", index=False)

    p = ERGEBNISSE / "Z_integration" / "Z_modul_gegen_marker.csv"
    if p.exists():
        M = pd.read_csv(p)
        M.to_csv(ZIEL / "f7_modul_gegen_marker.csv", index=False)
        log("\nf7_modul_gegen_marker.csv -- %d rows" % len(M))
    else:
        log("\nWARNING: Z_modul_gegen_marker.csv still missing.")

    # ---- scatter plot of the replicated chromatin layer
    p = ERGEBNISSE / "B_atac" / "B2_GSE224251_genwerte_hart.csv"
    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    K = K[K.symbol.notna()].drop_duplicates("symbol").set_index("symbol")
    if p.exists():
        G = pd.read_csv(p, index_col=0)
        gem = G.index.intersection(K.index)
        S = pd.DataFrame({
            "symbol": gem,
            "dWT_kons": K.dWT_kons.reindex(gem).values,
            "iv_kons": K.iv_kons.reindex(gem).values,
            "basis": K.basis_med.reindex(gem).values,
            "atac": G.osteogen.reindex(gem).values,
            "im_module": K.im_module.reindex(gem).values,
            "ri": K.ri.reindex(gem).values,
        })
        S.to_csv(ZIEL / "f7_streu_atac.csv", index=False)
        log("\nf7_streu_atac.csv -- %d genes, of which %d in the module"
            % (len(S), int(S.im_module.sum())))

    (ERGEBNISSE / "Z_integration" / "f7_daten_log.txt").write_text(
        "\n".join(LOG) + "\n", encoding="utf-8")
    print("\nwritten to", ZIEL)


def _lauf() -> None:
    main()
    robustheit()


def robustheit() -> None:
    """f7_robustheit.csv -- the module finding across all eight preparations
    of the replicated chromatin layer. Shows that it is not an artifact of a
    pipeline decision."""
    p = ERGEBNISSE / "B_atac" / "B2_GSE224251_modultest_hart.csv"
    if not p.exists():
        return
    R = pd.read_csv(p)
    R["peaks_kurz"] = R.peaks.map({"alle Peaks": "alle Peaks",
                                   "nur Promotorpeaks": "Promotor"})
    R[["schluessel", "peaks", "peaks_kurz", "filter", "normierung", "n",
       "konkordanz", "konkordanz_null", "konkordanz_null_sd", "konkordanz_z",
       "konkordanz_p", "rang", "rang_z", "rang_p",
       "wilson_lo", "wilson_hi"]].to_csv(ZIEL / "f7_robustheit.csv", index=False)
    print("f7_robustheit.csv -- %d preparations, continuous z from %+.2f to %+.2f"
          % (len(R), R.rang_z.min(), R.rang_z.max()))


if __name__ == "__main__":
    _lauf()

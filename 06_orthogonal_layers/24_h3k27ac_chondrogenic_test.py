# -*- coding: utf-8 -*-
"""
24_h3k27ac_chondrogenic_test.py -- layer B3: calibration and module test on H3K27ac.

Input: derived_data/B_atac/B3_GSE129031_matrix_<Fenster>.csv from 27_.

Axis: `CHON` (day-14 chondrocyte) against `MSC` (naive), formed per donor
line and then averaged over the two lines. Two biological units --
descriptive at cohort level, inferential only at gene level.

Sign rule, fixed before computation: H3K27ac marks active promoters and
enhancers and runs in the SAME direction as transcription. Expected sign is
`+ri`, as for accessibility and unlike methylation.

The purpose of this layer is the LAYER COMPARISON WITHIN ONE AXIS:
GSE129266 measures the methylome of the same in-vitro chondrogenesis. If
the module shows up here but not there, the difference lies in the
measurement layer, not in the differentiation axis.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
from _marker import CHONDROGEN, NAIV  # noqa: E402
from _module import ERGEBNISSE, MODUL, konkordanz, kontrast, wilson  # noqa: E402

AUS = ERGEBNISSE / "B_atac"
LOG: list[str] = []
EPS = 0.02
FENSTER = {"P": "promoter, TSS -2000/+500", "T10": "TSS +- 10 kb",
           "T50": "TSS +- 50 kb"}
LINIEN = ["8A", "2454e"]


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def achse(M: pd.DataFrame) -> pd.DataFrame:
    lr = {}
    for L in LINIEN:
        lr[L] = np.log2((M["CHON_%s" % L] + EPS) / (M["MSC_%s" % L] + EPS))
    LR = pd.DataFrame(lr)
    D = pd.DataFrame({"chondrogen": LR.mean(axis=1)})
    vz = np.sign(LR)
    D["linienkonsistent"] = (vz.abs().sum(axis=1) == 2) & (vz.sum(axis=1).abs() == 2)
    D["basis"] = M[["MSC_%s" % L for L in LINIEN]].mean(axis=1)
    for L in LINIEN:
        D["linie_" + L] = LR[L]
    return D


def main() -> None:
    log("=" * 78)
    log("Layer B3 -- H3K27ac, GSE129031, MSC vs chondrocyte (2 donor lines)")
    log("=" * 78)
    log("Purpose: the same layer comparison WITHIN the chondrogenic axis,")
    log("against the methylome of the same differentiation (GSE129266).")
    log("Donor labels of the two series do NOT match -- it is the same system")
    log("and lab, but not verifiably the same donors.")

    sym_ri = dict(zip(MODUL.symbol, MODUL.ri))
    eich, zeilen = [], []
    daten = {}

    for f, besch in FENSTER.items():
        M = pd.read_csv(AUS / ("B3_GSE129031_H3K27ac_matrix_%s.csv" % f), index_col=0)
        D = achse(M)
        daten[f] = D

        # INPUT CONTROL. The input measures mappability, copy number and
        # fragmentation propensity. An H3K27ac difference that also appears
        # in the input is not a difference of the mark. The axis is
        # therefore additionally computed on the mark/input ratio per sample.
        pi_ = AUS / ("B3_GSE129031_input_matrix_%s.csv" % f)
        D_in = None
        if pi_.exists():
            I = pd.read_csv(pi_, index_col=0)
            gem = M.index.intersection(I.index)
            R = pd.DataFrame({c: (M.loc[gem, c] + EPS) / (I.loc[gem, c] + EPS)
                              for c in M.columns})
            D_in = achse(R)
            daten[f + "_input"] = D_in
        log("\n" + "=" * 78)
        log("Window %s (%s) -- %d genes" % (f, besch, len(D)))
        r_lin = stats.pearsonr(D["linie_8A"], D["linie_2454e"])
        log("Agreement of the two donor lines: r = %.3f" % r_lin.statistic)
        log("line-consistent genes: %d (%.1f %%)"
            % (D.linienkonsistent.sum(), 100 * D.linienkonsistent.mean()))

        # ------------------------------------------------------- calibration
        hg = D.chondrogen.dropna()
        e = kontrast(hg, CHONDROGEN, NAIV)
        ok_e = e.get("status") == "ok" and e["kontrast"] >= e["mde80"]
        log("\nCALIBRATION (chondrogenic minus naive markers, %d vs %d genes):"
            % (e.get("n_a", 0), e.get("n_b", 0)))
        log("  contrast %+.3f | null %+.3f +- %.3f | MDE80 %+.3f | z %+.2f | p %.4g | %s"
            % (e["kontrast"], e["null_mittel"], e["null_sd"], e["mde80"],
               e["z"], e["p"], "PASSED" if ok_e else "failed"))
        eich.append(dict(datensatz="GSE129031", fenster=f, bestanden=ok_e,
                         **{k: v for k, v in e.items() if k != "status"}))

        # ------------------------------------------------------ module test
        mod = hg[hg.index.isin(sym_ri)]
        erw = pd.Series({s: sym_ri[s] for s in mod.index})
        log("\nMODULE TEST (expected sign +ri), %d of 173 genes measurable:"
            % len(mod))
        for hname, sch in (("Hintergrund", None), ("basisgeschichtet", D["basis"])):
            res = konkordanz(mod, erw, hintergrund=hg, schichtung=sch)
            k = int((np.sign(mod.values) == erw.reindex(mod.index).values).sum())
            lo, hi = wilson(k, len(mod))
            schwelle = res["konkordanz"] >= res["konkordanz_mde80"]
            log("  %-18s C %.3f [%.3f-%.3f] | null %.3f+-%.3f | MDE80 %.3f | "
                "z %+5.2f | p %.4g | %s"
                % (hname, res["konkordanz"], lo, hi, res["konkordanz_null"],
                   res["konkordanz_null_sd"], res["konkordanz_mde80"],
                   res["konkordanz_z"], res["konkordanz_p"],
                   "THRESHOLD REACHED" if schwelle else "below threshold"))
            log("  %-18s continuous %+.4f | z %+5.2f | p %.4g"
                % ("", res["rang"], res["rang_z"], res["rang_p"]))
            zeilen.append(dict(ebene="H3K27ac", datensatz="GSE129031", fenster=f,
                               achse="chondrogen", null=hname, geeicht=ok_e,
                               k=k, wilson_lo=lo, wilson_hi=hi,
                               schwelle_erreicht=bool(schwelle),
                               **{x: y for x, y in res.items() if x != "status"}))

        # ---------------------------------------------- input control
        if D_in is not None:
            hg_i = D_in.chondrogen.dropna()
            e_i = kontrast(hg_i, CHONDROGEN, NAIV)
            mod_i = hg_i[hg_i.index.isin(sym_ri)]
            erw_i = pd.Series({s2: sym_ri[s2] for s2 in mod_i.index})
            r_i = konkordanz(mod_i, erw_i, hintergrund=hg_i,
                             schichtung=D_in["basis"])
            log("  [input control, mark/input] calibration z %+.2f | "
                "module C %.3f | null %.3f | z %+5.2f | p %.4g | continuous z %+.2f"
                % (e_i.get("z", np.nan), r_i["konkordanz"], r_i["konkordanz_null"],
                   r_i["konkordanz_z"], r_i["konkordanz_p"], r_i["rang_z"]))
            zeilen.append(dict(ebene="H3K27ac/Input", datensatz="GSE129031",
                               fenster=f, achse="chondrogen",
                               null="basisgeschichtet",
                               geeicht=bool(e_i.get("status") == "ok"
                                            and e_i["kontrast"] >= e_i["mde80"]),
                               k=np.nan, wilson_lo=np.nan, wilson_hi=np.nan,
                               schwelle_erreicht=bool(
                                   r_i["konkordanz"] >= r_i["konkordanz_mde80"]),
                               **{x: y for x, y in r_i.items() if x != "status"}))

        # per donor line separately
        for L in LINIEN:
            hg2 = D["linie_" + L].dropna()
            mod2 = hg2[hg2.index.isin(sym_ri)]
            erw2 = pd.Series({s: sym_ri[s] for s in mod2.index})
            r2 = konkordanz(mod2, erw2, hintergrund=hg2)
            log("  [line %-6s] C %.3f | null %.3f | z %+5.2f | p %.4g"
                % (L, r2["konkordanz"], r2["konkordanz_null"],
                   r2["konkordanz_z"], r2["konkordanz_p"]))
            zeilen.append(dict(ebene="H3K27ac", datensatz="GSE129031", fenster=f,
                               achse="Linie " + L, null="Hintergrund",
                               geeicht=ok_e, k=np.nan, wilson_lo=np.nan,
                               wilson_hi=np.nan,
                               schwelle_erreicht=bool(
                                   r2["konkordanz"] >= r2["konkordanz_mde80"]),
                               **{x: y for x, y in r2.items() if x != "status"}))

    pd.DataFrame(eich).to_csv(AUS / "B3_GSE129031_eichung.csv", index=False)
    pd.DataFrame(zeilen).to_csv(AUS / "B3_GSE129031_modultest.csv", index=False)

    # reference window for the integration: T10 (H3K27ac lies promoter-proximal
    # AND distal; T10 is the compromise and is chosen BEFORE the result,
    # because it matches the window of the accessibility layer).
    D = daten["T10"].copy()
    D["im_module"] = D.index.isin(sym_ri)
    D["ri"] = [sym_ri.get(s, np.nan) for s in D.index]
    D.index.name = "symbol"
    D.to_csv(AUS / "B3_GSE129031_genwerte_T10.csv")

    log("\n" + "=" * 78)
    log("OVERVIEW")
    log("=" * 78)
    Z = pd.DataFrame(zeilen)
    Z = Z[(Z.achse == "chondrogen") & (Z.null == "basisgeschichtet")]
    log("%-5s %8s %9s %9s %10s" % ("win", "C", "z_con", "z_cont", "p_cont"))
    for _, r in Z.iterrows():
        log("%-5s %8.3f %+9.2f %+9.2f %10.4g"
            % (r.fenster, r.konkordanz, r.konkordanz_z, r.rang_z, r.rang_p))

    (AUS / "B3_GSE129031_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

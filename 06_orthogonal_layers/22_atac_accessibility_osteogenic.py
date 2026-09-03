# -*- coding: utf-8 -*-
"""
22_atac_accessibility_osteogenic.py -- layer B (main), hardened version.

Two problems of the first version (25_) are fixed here:

  (1) DEPTH IMBALANCE. Library sizes scatter by a factor of eight
      (312 279 to 2 509 594 counts in peaks), and the two thinnest samples
      BOTH lie in the osteogenic flat arm. Pure CPM normalization leaves
      this inequality in the null distribution: with a thin library many
      genes sit at zero, `log2(CPM + 1)` compresses them, and the
      difference `osteogen minus naiv` acquires a depth-dependent skew.
      Remedy: (a) peak filter for sufficient total coverage,
      (b) quantile normalization of the per-gene log matrix, forcing every
      sample onto the same marginal distribution.

  (2) CALIBRATION. The marker set `naiv` in `_marker.py` mixes surface
      markers (THY1, ENG, NT5E) and proliferation markers (MKI67, TOP2A).
      Neither must close in chromatin during a seven-day osteogenic
      induction. The calibration is therefore additionally computed in the
      form that is testable on this layer at all: osteogenic markers
      against the BACKGROUND (enrichment test), instead of against a
      second, possibly unsuitable marker set.

Both normalizations and both calibrations are reported side by side. There
is no selection of whichever yields the friendlier result.
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
from _marker import NAIV, OSTEOGEN  # noqa: E402
from _module import DATEN, ERGEBNISSE, MODUL, konkordanz, kontrast, wilson  # noqa: E402

AUS = ERGEBNISSE / "B_atac"
LOG: list[str] = []

PROBEN = {
    "SF": ["SF.1_S1", "SF.2_S2", "SF.3_S3"],
    "SP": ["SP.1_S4", "SP.2_S5", "SP.3_S6"],
    "OF": ["OF.1_S7", "OF.2_S8", "OF.3_S9"],
    "OP": ["OP.1_S10", "OP.2_S11", "OP.3_S12"],
}
SPALTEN = [s for gr in PROBEN.values() for s in gr]
MIN_ZAEHL = 20      # minimum total of a peak over all twelve samples


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def quantilnorm(L: pd.DataFrame) -> pd.DataFrame:
    """Force every column onto the common marginal distribution."""
    R = L.rank(axis=0, method="average")
    ziel = np.sort(L.values, axis=0).mean(axis=1)
    out = pd.DataFrame(index=L.index, columns=L.columns, dtype=float)
    n = len(L)
    for c in L.columns:
        pos = np.clip((R[c].values - 1), 0, n - 1)
        lo = np.floor(pos).astype(int)
        hi = np.ceil(pos).astype(int)
        frac = pos - lo
        out[c] = ziel[lo] * (1 - frac) + ziel[hi] * frac
    return out


def achse(L: pd.DataFrame) -> pd.DataFrame:
    d_f = L[PROBEN["OF"]].mean(axis=1) - L[PROBEN["SF"]].mean(axis=1)
    d_p = L[PROBEN["OP"]].mean(axis=1) - L[PROBEN["SP"]].mean(axis=1)
    D = pd.DataFrame({"flach": d_f, "saeule": d_p})
    D["osteogen"] = D[["flach", "saeule"]].mean(axis=1)
    D["basis"] = L[PROBEN["SF"] + PROBEN["SP"]].mean(axis=1)
    return D


def main() -> None:
    log("=" * 78)
    log("Layer B (main, hardened) -- ATAC GSE224251, hMSC naive vs osteogenic")
    log("=" * 78)

    T = pd.read_csv(DATEN / "GSE224251" / "GSE224251_count_table.csv.gz",
                    low_memory=False)
    T = T[T.symbol.notna() & (T.symbol.astype(str) != "nan")].copy()
    T["promotor"] = T.annotation.astype(str).str.startswith("Promoter")
    T["summe"] = T[SPALTEN].sum(axis=1)

    log("Peaks with gene annotation: %d" % len(T))
    log("of these with total coverage >= %d: %d (%.1f %%)"
        % (MIN_ZAEHL, (T.summe >= MIN_ZAEHL).sum(),
           100 * (T.summe >= MIN_ZAEHL).mean()))
    bib = T[SPALTEN].sum(axis=0)
    log("Library sizes: %.0f to %.0f (factor %.1f)"
        % (bib.min(), bib.max(), bib.max() / bib.min()))

    sym_ri = dict(zip(MODUL.symbol, MODUL.ri))
    zeilen, eichzeilen = [], []
    daten = {}

    for pname, pmaske in (("alle Peaks", pd.Series(True, index=T.index)),
                          ("nur Promotorpeaks", T.promotor)):
        for fname, fmaske in (("ungefiltert", pd.Series(True, index=T.index)),
                              ("Peakfilter", T.summe >= MIN_ZAEHL)):
            sub = T[pmaske & fmaske]
            M = sub.groupby("symbol")[SPALTEN].sum()
            L0 = np.log2(M.div(bib, axis=1) * 1e6 + 1.0)
            for nname, L in (("CPM", L0), ("CPM+Quantil", quantilnorm(L0))):
                schluessel = "%s | %s | %s" % (pname, fname, nname)
                D = achse(L)
                daten[schluessel] = D

                # ---- calibration 1: markers against markers
                e1 = kontrast(D.osteogen, OSTEOGEN, NAIV)
                ok1 = e1.get("status") == "ok" and e1["kontrast"] >= e1["mde80"]
                # ---- calibration 2: osteogenic markers against the background
                hg = D.osteogen.dropna()
                mk = [g for g in OSTEOGEN if g in hg.index]
                rng = np.random.default_rng(20260821)
                beob = hg[mk].mean()
                null = np.array([hg.values[rng.choice(len(hg), len(mk), False)].mean()
                                 for _ in range(20000)])
                z2 = (beob - null.mean()) / null.std(ddof=1)
                p2 = min(1.0, 2 * (1 + (null >= beob).sum()) / (1 + len(null)))
                ok2 = beob >= null.mean() + 2.8 * null.std(ddof=1)

                # ---- module test
                mod = hg[hg.index.isin(sym_ri)]
                erw = pd.Series({s: sym_ri[s] for s in mod.index})
                res = konkordanz(mod, erw, hintergrund=hg, schichtung=D["basis"])
                k = int((np.sign(mod.values) == erw.reindex(mod.index).values).sum())
                lo, hi = wilson(k, len(mod))

                log("\n--- %s ---" % schluessel)
                log("  genes %d | agreement of the surfaces r = %.3f"
                    % (len(D), stats.pearsonr(D.flach, D.saeule).statistic))
                log("  calibration 1 (osteo minus naiv):     %+7.3f | z %+5.2f | p %.4g | %s"
                    % (e1.get("kontrast", np.nan), e1.get("z", np.nan),
                       e1.get("p", np.nan), "passed" if ok1 else "failed"))
                log("  calibration 2 (osteo vs background): %+7.3f | z %+5.2f | p %.4g | %s"
                    % (beob - null.mean(), z2, p2,
                       "passed" if ok2 else "failed"))
                log("  MODULE n %3d | C %.3f [%.3f-%.3f] | null %.3f+-%.3f | "
                    "z %+5.2f | p %.4g"
                    % (res["n"], res["konkordanz"], lo, hi, res["konkordanz_null"],
                       res["konkordanz_null_sd"], res["konkordanz_z"],
                       res["konkordanz_p"]))
                log("         continuous %+.4f | z %+5.2f | p %.4g"
                    % (res["rang"], res["rang_z"], res["rang_p"]))

                eichzeilen.append(dict(schluessel=schluessel, peaks=pname,
                                       filter=fname, normierung=nname,
                                       eichung1_kontrast=e1.get("kontrast"),
                                       eichung1_z=e1.get("z"), eichung1_p=e1.get("p"),
                                       eichung1_ok=ok1,
                                       eichung2_diff=float(beob - null.mean()),
                                       eichung2_z=float(z2), eichung2_p=float(p2),
                                       eichung2_ok=bool(ok2)))
                zeilen.append(dict(ebene="ATAC", datensatz="GSE224251",
                                   schluessel=schluessel, peaks=pname, filter=fname,
                                   normierung=nname, k=k, wilson_lo=lo, wilson_hi=hi,
                                   **{x: y for x, y in res.items() if x != "status"}))

    Eich = pd.DataFrame(eichzeilen)
    Z = pd.DataFrame(zeilen)
    Eich.to_csv(AUS / "B2_GSE224251_eichung_hart.csv", index=False)
    Z.to_csv(AUS / "B2_GSE224251_modultest_hart.csv", index=False)

    log("\n" + "=" * 78)
    log("OVERVIEW over all eight preparations")
    log("=" * 78)
    log("%-42s %7s %7s %8s %7s" % ("preparation", "C", "z_con", "z_cont", "p_cont"))
    for _, r in Z.iterrows():
        log("%-42s %7.3f %+7.2f %+8.2f %7.4g"
            % (r.schluessel, r.konkordanz, r.konkordanz_z, r.rang_z, r.rang_p))
    log("\nRange of the continuous z over all preparations: %+.2f to %+.2f"
        % (Z.rang_z.min(), Z.rang_z.max()))
    log("All eight preparations show the same sign: %s"
        % ("yes" if (np.sign(Z.rang_z) == np.sign(Z.rang_z.iloc[0])).all() else "no"))

    best = "nur Promotorpeaks | Peakfilter | CPM+Quantil"
    D = daten[best].copy()
    D["im_module"] = D.index.isin(sym_ri)
    D["ri"] = [sym_ri.get(s, np.nan) for s in D.index]
    D["markerklasse"] = ["osteogen" if s in OSTEOGEN else
                         "naiv" if s in NAIV else "" for s in D.index]
    D.index.name = "symbol"
    D.to_csv(AUS / "B2_GSE224251_genwerte_hart.csv")
    log("\ngene values of the reference preparation written: %s" % best)

    (AUS / "B2_GSE224251_hart_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

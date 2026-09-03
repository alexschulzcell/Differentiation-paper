# -*- coding: utf-8 -*-
"""
32_dose_response.py -- the dose-response in the form that carries the claim.

The decile curve in 30_ does not carry the claim: it puts the SIGNED
convergence on the x axis and the raw layer difference on the y axis.
Because the orthogonal layer has a global offset (during differentiation
more chromatin opens than closes), all points then lie above zero, and the
curve rises equally for `dWT` and `iv`.

The right form asks for the STRENGTH of convergence, not its direction:

    x  = decile of the convergence strength |kons|, from "no common sign"
         to "identical in all eighteen datasets"
    y  = mean DIRECTED layer value, sign(kons) * (rank(delta) - 0.5),
         centered on the layer background

This removes the global offset of the layer (it affects both signs alike),
and y = 0 means "the layer knows nothing". A rising curve means: the more
the eighteen transcriptome datasets agree on a gene, the more reliably it
moves on the foreign layer in exactly the predicted direction.

Additionally, per decile a spread band from 2 000 bootstrap draws of the
genes is computed -- the genes are the unit, not the samples.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import ERGEBNISSE, SEED, WURZEL  # noqa: E402

import importlib.util  # noqa: E402
spec = importlib.util.spec_from_file_location(
    "i30", pathlib.Path(__file__).resolve().parent / "30_convergence_dose_integration.py")
i30 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(i30)

AUS = ERGEBNISSE / "Z_integration"
ZIEL = WURZEL / "derived_data" / "manuscript"
NBOOT = 2000
NDEZ = 8
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def main() -> None:
    log("=" * 78)
    log("Dose-response: convergence strength against directed layer movement")
    log("=" * 78)

    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    K = K[K.symbol.notna()].drop_duplicates("symbol").set_index("symbol")
    E = i30.lade_ebenen()

    rng = np.random.default_rng(SEED)
    zeilen = []
    for ename, delta in E.items():
        gem = delta.index.intersection(K.index)
        if len(gem) < 500:
            continue
        d = delta.reindex(gem)
        r = pd.Series(stats.rankdata(d.values) / len(d) - 0.5, index=gem)

        # BASELINE CONTROL. The layer value tracks the baseline level: highly
        # expressed genes have more open chromatin and more H3K27ac. Because
        # the internal convergence genes start systematically LOW (the
        # confounder from Figure 2), a dose-response could arise from that
        # alone. The computation therefore uses the RESIDUAL of the layer
        # rank after linear adjustment on the baseline rank -- the same
        # control that 30_ performs via the stratified null, here in a form
        # that allows drawing a curve.
        b = K.basis_med.reindex(gem)
        ok = b.notna() & r.notna()
        br = pd.Series(np.nan, index=gem)
        br[ok] = stats.rankdata(b[ok].values) / int(ok.sum()) - 0.5
        gp = np.polyfit(br[ok].values, r[ok].values, 1)
        r_res = pd.Series(np.nan, index=gem)
        r_res[ok] = r[ok].values - np.polyval(gp, br[ok].values)
        log("\n--- %s ---" % ename)
        log("%6s %8s %7s %10s %10s %10s"
            % ("dose", "decile", "n", "|kons|", "directed", "95 % CI"))
        for dosis in ("dWT_kons", "iv_kons"):
            kons = K[dosis].reindex(gem)
            T = pd.DataFrame({"kons": kons, "r": r, "r_res": r_res}).dropna()
            T["staerke"] = T.kons.abs()
            T["gerichtet_roh"] = np.sign(T.kons) * T.r
            T["gerichtet"] = np.sign(T.kons) * T.r_res   # baseline-controlled
            T = T[T.kons != 0]
            T["dezil"] = pd.qcut(T.staerke.rank(method="first"), NDEZ, labels=False)
            for bb, g in T.groupby("dezil"):
                v = g.gerichtet.values
                bs = np.array([v[rng.integers(0, len(v), len(v))].mean()
                               for _ in range(NBOOT)])
                lo, hi = np.quantile(bs, [0.025, 0.975])
                zeilen.append(dict(ebene=ename, dosis=dosis, dezil=int(bb),
                                   n=len(v), staerke=float(g.staerke.median()),
                                   staerke_min=float(g.staerke.min()),
                                   staerke_max=float(g.staerke.max()),
                                   gerichtet=float(v.mean()),
                                   gerichtet_roh=float(g.gerichtet_roh.mean()),
                                   lo=float(lo), hi=float(hi)))
                log("%6s %8d %7d %10.3f %+10.4f  %+.4f..%+.4f"
                    % (dosis.replace("_kons", ""), bb, len(v),
                       g.staerke.median(), v.mean(), lo, hi))

    D = pd.DataFrame(zeilen)
    D.to_csv(AUS / "Z_dosiskurve_gerichtet.csv", index=False)
    D.to_csv(ZIEL / "f7_dosiskurve.csv", index=False)

    # --------------------------------------------- trend test per layer/dose
    log("\n" + "=" * 78)
    log("Trend across the deciles (Spearman between decile rank and y)")
    log("=" * 78)
    log("%-36s %-6s %8s %10s %12s"
        % ("layer", "dose", "rho", "p", "y top"))
    tr = []
    for (ename, dosis), g in D.groupby(["ebene", "dosis"]):
        s = stats.spearmanr(g.dezil, g.gerichtet)
        log("%-36s %-6s %+8.3f %10.4g %+12.4f"
            % (ename, dosis.replace("_kons", ""), s.statistic, s.pvalue,
               g.gerichtet.iloc[-1]))
        tr.append(dict(ebene=ename, dosis=dosis, trend_rho=s.statistic,
                       trend_p=s.pvalue, y_oberstes=g.gerichtet.iloc[-1],
                       lo_oberstes=g.lo.iloc[-1], hi_oberstes=g.hi.iloc[-1]))
    pd.DataFrame(tr).to_csv(AUS / "Z_dosiskurve_trend.csv", index=False)
    pd.DataFrame(tr).to_csv(ZIEL / "f7_dosistrend.csv", index=False)

    (AUS / "Z_dosiskurve_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\nwritten to", AUS, "and", ZIEL)


if __name__ == "__main__":
    main()

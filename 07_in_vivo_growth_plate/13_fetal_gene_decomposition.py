# -*- coding: utf-8 -*-
"""
13_fetal_gene_decomposition.py -- is the in vivo trend broadly carried?

Only 76 of 164 module genes (46.3 %, CI 38.7-54.0 %) are individually
concordant in vivo, and that interval covers chance. A clean split into
"confirmed in vivo" and "in vitro only" cannot be derived from that number;
what it needs is a ranking of individual genes by |Delta|.

This is that ranking.

-----------------------------------------------------------------------------
THE DECISION RULE -- written here before the first number was computed, and
not changed afterwards.
-----------------------------------------------------------------------------
Question: is the aggregate trend (rho 0.456, z +4.80, limit rho 0.274,
`results/invivo_spendertest.csv`) carried by a MINORITY of large, aligned
excursions, or by a broad majority of small ones?

Test: the same trend test, with the same null (permutation of zone labels
WITHIN the specimen, `07_in_vivo_growth_plate/12_fetal_donor_trend_test.py`, imported unchanged), on
a module from which the **10 % of genes with the largest |Delta|** have been
removed. `Delta` is `delta_hyper_minus_mes` from
`derived_data/followup/ws4_p2_gen_konkordanz.csv`; 10 % means ceil(0.10 * n)
among the genes with a finite Delta. Genes are removed from BOTH direction
sets together, by |Delta|, without regard to whether a gene is concordant --
otherwise the selection would already be an answer.

  * If the trend STAYS above its own limit, it is BROADLY CARRIED, and the
    legend of Figure 3C says so.
  * If it FALLS below the limit, Figure 3C is reported as carried by a few
    genes and the discussion is weakened accordingly.

The limit is determined AFRESH in the reduced run, from that run's own null --
a reduced module has its own detection limit (project rule 1).

In addition, DESCRIPTIVE and explicitly not part of the rule: the same
computation for 5, 20 and 30 % of genes removed, and the distribution of Delta
by direction. They describe; they decide nothing.

-----------------------------------------------------------------------------
NO SECOND IMPLEMENTATION. The contrast comes from
`00_shared/_module.py`, the trend test from
`07_in_vivo_growth_plate/12_fetal_donor_trend_test.py`, the pseudobulk matrix from
`07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py`. As a SELF-TEST the module contrast of the
full set is compared against
`derived_data/followup/ws4_modulwert_je_probe.csv`; if it does not match, the
script stops.

Inputs   results/invivo_pseudobulk.csv.gz,
         derived_data/followup/ws4_p2_gen_konkordanz.csv,
         derived_data/followup/ws4_modulwert_je_probe.csv
Outputs  results/invivo_genzerlegung.csv, results/invivo_genzerlegung_log.txt
Runtime  a few minutes
"""
from __future__ import annotations

import importlib.util
import math
import os
import pathlib
import sys

import numpy as np
import pandas as pd

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
sys.path.insert(0, str(WURZEL / "00_shared"))
from _module import MODUL, kontrast  # noqa: E402

NEU = WURZEL / "derived_data" / "followup"
RES = WURZEL / "results"

NZIEHUNGEN = 20000
SEED = 20260822          # as in 10_fetal_limb_atlas_pseudobulk_build.py
ANTEILE = [0.10, 0.05, 0.20, 0.30]   # 0.10 is the rule, the rest descriptive

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def trendtest():
    """The donor-stratified trend test from 07_in_vivo_growth_plate/12_fetal_donor_trend_test.py --
    imported, not rebuilt."""
    p = WURZEL / "07_in_vivo_growth_plate" / "12_fetal_donor_trend_test.py"
    spec = importlib.util.spec_from_file_location("_spendertest", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.test


def modulkontrast(PB: pd.DataFrame, gene: list) -> pd.DataFrame:
    """Module contrast per (zone, sample) on a given gene set."""
    up = [g for g in MODUL.loc[MODUL.ri > 0, "symbol"] if g in gene]
    dn = [g for g in MODUL.loc[MODUL.ri < 0, "symbol"] if g in gene]
    gen_spalten = [c for c in PB.columns
                   if c not in ("zone", "probe", "n_zellen")]
    zeilen = []
    for _, r in PB.iterrows():
        werte = r[gen_spalten]
        res = kontrast(werte, up, dn, nziehungen=NZIEHUNGEN, seed=SEED)
        res.update({"zone": r.zone, "probe": r.probe})
        zeilen.append(res)
    return pd.DataFrame(zeilen)


def main() -> None:
    log("=" * 78)
    log("Is the in vivo trend broadly carried? Decision rule in the header.")
    log("=" * 78)

    PB = pd.read_csv(RES / "invivo_pseudobulk.csv.gz")
    GK = pd.read_csv(NEU / "ws4_p2_gen_konkordanz.csv")
    test = trendtest()

    alle = MODUL["symbol"].tolist()
    im_atlas = [g for g in alle if g in PB.columns]
    log("Module genes measurable in the atlas: %d of %d"
        % (len(im_atlas), len(alle)))

    # ---- self-test: the full set must reproduce the stored values ----------
    voll = modulkontrast(PB, im_atlas)
    ws4 = pd.read_csv(NEU / "ws4_modulwert_je_probe.csv")
    j = voll.merge(ws4[["zone", "probe", "kontrast"]], on=["zone", "probe"],
                   suffixes=("", "_ws4"))
    d = float((j.kontrast - j.kontrast_ws4).abs().max())
    log("self-test against ws4_modulwert_je_probe.csv: max |Delta| = %.3e "
        "over %d points" % (d, len(j)))
    # Threshold 1e-6 rather than 0: the pseudobulk values are float32 means
    # passed through a text file; the remaining deviation is 3e-8, that is,
    # a matter of representation and not of the computation.
    if d > 1e-6:
        raise SystemExit("The pseudobulk matrix is not the one of WS4 "
                         "(max |Delta| %.3g) -- abort." % d)

    t_voll = test(voll, "Modul voll (%d Gene)" % len(im_atlas))
    log("Full set: rho %.4f | z %+.3f | limit rho %.4f | above limit %s"
        % (t_voll["rho"], t_voll["z"], t_voll["mde80_rho"],
           t_voll["ueber_mde80"]))

    # ---- the distribution of Delta -----------------------------------------
    g = GK.dropna(subset=["delta_hyper_minus_mes"]).copy()
    g = g[g.symbol.isin(im_atlas)]
    g["abs_delta"] = g.delta_hyper_minus_mes.abs()
    g = g.sort_values("abs_delta", ascending=False)
    ges = float(g.abs_delta.sum())
    log("")
    log("Distribution of |Delta| over %d measurable module genes:" % len(g))
    for k in (10, 17, 25, 50):
        log("   the %2d largest carry %5.1f %% of the sum of |Delta|"
            % (k, 100 * g.abs_delta.head(k).sum() / ges))
    log("   median |Delta| %.4f | 90th percentile %.4f | maximum %.4f"
        % (g.abs_delta.median(), g.abs_delta.quantile(0.9), g.abs_delta.max()))
    kon = GK.konkordant.astype("string").str.lower()
    log("   individually concordant: %d of %d (%.1f %%)"
        % (int((kon == "true").sum()), int(kon.isin(["true", "false"]).sum()),
           100 * (kon == "true").sum() / max(1, kon.isin(["true", "false"]).sum())))
    auf = g[g.ri > 0].delta_hyper_minus_mes
    ab = g[g.ri < 0].delta_hyper_minus_mes
    log("   median Delta at ri=+1: %+.4f (n %d) | at ri=-1: %+.4f (n %d)"
        % (auf.median(), len(auf), ab.median(), len(ab)))

    # ---- the rule -----------------------------------------------------------
    zeilen = [dict(t_voll, anteil_entfernt=0.0, n_entfernt=0,
                   rolle="Ausgangswert")]
    for a in ANTEILE:
        k = math.ceil(a * len(g))
        raus = set(g.symbol.head(k))
        rest = [x for x in im_atlas if x not in raus]
        kt = modulkontrast(PB, rest)
        t = test(kt, "Modul ohne die %.0f %% groessten |Delta| (%d Gene)"
                 % (100 * a, len(rest)))
        t.update(anteil_entfernt=a, n_entfernt=k,
                 rolle="Regel" if abs(a - 0.10) < 1e-9 else "deskriptiv")
        zeilen.append(t)
        log("without the %2.0f %% largest |Delta| (%3d genes out, %3d left): "
            "rho %.4f | z %+.3f | limit %.4f | above limit %s   [%s]"
            % (100 * a, k, len(rest), t["rho"], t["z"], t["mde80_rho"],
               t["ueber_mde80"], t["rolle"]))

    T = pd.DataFrame(zeilen)
    T.to_csv(RES / "invivo_genzerlegung.csv", index=False)

    regel = T[T.anteil_entfernt == 0.10].iloc[0]
    log("")
    log("-" * 78)
    log("DECISION per the rule in the script header:")
    if bool(regel.ueber_mde80):
        log("  The trend stays above its own limit "
            "(rho %.4f > %.4f, z %+.2f)." % (regel.rho, regel.mde80_rho,
                                             regel.z))
        log("  -> BROADLY CARRIED. The caption of Fig. 3C says so.")
    else:
        log("  The trend falls below its own limit "
            "(rho %.4f <= %.4f, z %+.2f)." % (regel.rho, regel.mde80_rho,
                                              regel.z))
        log("  -> CARRIED BY FEW GENES. Fig. 3C is reported as such and")
        log("     the discussion §3 is weakened.")
    log("-" * 78)

    (RES / "invivo_genzerlegung_log.txt").write_text(
        chr(10).join(LOG) + chr(10), encoding="utf-8")
    log("-> %s" % (RES / "invivo_genzerlegung.csv"))


if __name__ == "__main__":
    main()

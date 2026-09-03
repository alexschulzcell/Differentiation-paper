# -*- coding: utf-8 -*-
"""
11_methylation_sensitivity.py -- layer A, follow-up check of the null finding.

10_methylation_osteogenic_27k.py finds no direction of the fixed 173-gene module on the
methylation layer. Before this is carried as a biological statement, three
technical explanations must be ruled out:

  (S1) DYNAMIC RANGE. Promoter CpG islands of active genes are
       constitutively unmethylated (beta ~ 0). A probe at the floor cannot
       show direction. Repeat on probes with dynamic range only
       (0.10 <= beta <= 0.90 in at least one sample).
  (S2) ISLAND STATUS. Outside islands lies the larger part of
       differentiation-dependent methylation dynamics. Repeat separately
       for island / non-island.
  (S3) GLOBAL DRIFT. The osteogenic-vs-myogenic comparison in 10_ is
       carried partly by the different global drifts of the two axes. The
       drift-controlled test draws the axis difference
       (delta_osteo - delta_myo) and confronts it with a background draw of
       the same difference.

Additionally the built-in positive control of the layer: known osteogenic
key genes (RUNX2, SP7, ALPL, BGLAP, IBSP, COL1A1) and myogenic ones
(MYOD1, MYOG, DES, MYF5) -- does the layer show any directed
differentiation response at all?
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import ERGEBNISSE, MODUL, NZIEHUNGEN, SEED, konkordanz  # noqa: E402
import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location(
    "a10", pathlib.Path(__file__).resolve().parent / "10_methylation_osteogenic_27k.py")
a10 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a10)

AUS = ERGEBNISSE / "A_dnam"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


OSTEO_MARKER = ["RUNX2", "SP7", "ALPL", "BGLAP", "IBSP", "COL1A1", "SPP1", "DLX5"]
MYO_MARKER = ["MYOD1", "MYOG", "DES", "MYF5", "MYH3", "TNNT1"]


def main() -> None:
    log("=" * 78)
    log("Layer A -- sensitivity and positive control (GSE33896)")
    log("=" * 78)

    man = a10.lade_manifest()
    beta, meta = a10.lade_beta()
    behalte = meta[meta.arm.isin(["naiv", "osteo", "myo"]) & meta.donor.isin([1, 2, 3])]
    beta = beta[behalte.gsm.tolist()].apply(pd.to_numeric, errors="coerce")
    man = man[man.sonde.isin(beta.index)].copy()
    man["abstand_tss"] = pd.to_numeric(man.abstand_tss, errors="coerce")
    man = man[man.abstand_tss.abs() <= 1500]

    # probe differences per axis, donor-paired
    def sondendiff(arm: str) -> pd.Series:
        sp = {}
        for d in (1, 2, 3):
            n = behalte[(behalte.donor == d) & (behalte.arm == "naiv")].gsm.iloc[0]
            t = behalte[(behalte.donor == d) & (behalte.arm == arm)].gsm.iloc[0]
            sp[d] = beta[t] - beta[n]
        return pd.DataFrame(sp).mean(axis=1)

    d_ost = sondendiff("osteo")
    d_myo = sondendiff("myo")

    bmin = beta.min(axis=1)
    bmax = beta.max(axis=1)
    dynamisch = (bmax >= 0.10) & (bmin <= 0.90)
    insel = man.set_index("sonde").cpg_insel.astype(str).str.upper().eq("TRUE")

    log("\nProbe balance (promoter-proximal, |d(TSS)| <= 1500 bp): %d" % len(man))
    log("  with dynamic range (0.10 <= beta <= 0.90 anywhere): %d (%.1f %%)"
        % (dynamisch.reindex(man.sonde).sum(),
           100 * dynamisch.reindex(man.sonde).mean()))
    log("  on CpG island: %d (%.1f %%)" % (insel.sum(), 100 * insel.mean()))

    # ----------------------------------------- positive control of the layer
    log("\n--- Positive control: does the layer know differentiation? -------")
    log("%-10s %8s %10s %10s" % ("gene", "probes", "d_osteo", "d_myo"))
    kontr = []
    for g in OSTEO_MARKER + MYO_MARKER:
        s = man.sonde[man.symbol == g]
        if len(s) == 0:
            log("%-10s %8s %10s %10s" % (g, "-", "n/a", "n/a"))
            continue
        do, dm = d_ost.reindex(s).mean(), d_myo.reindex(s).mean()
        log("%-10s %8d %+10.4f %+10.4f" % (g, len(s), do, dm))
        kontr.append(dict(gen=g, klasse="osteogen" if g in OSTEO_MARKER else "myogen",
                          n_sonden=len(s), d_osteo=do, d_myo=dm))
    pd.DataFrame(kontr).to_csv(AUS / "A_dnam_positivkontrolle.csv", index=False)
    K = pd.DataFrame(kontr)
    if len(K):
        log("\nExpectation: osteogenic markers hypomethylated in the osteogenic axis,")
        log("myogenic markers hypomethylated in the myogenic axis.")
        for kl, sp_eigen, sp_fremd in (("osteogen", "d_osteo", "d_myo"),
                                       ("myogen", "d_myo", "d_osteo")):
            sub = K[K.klasse == kl]
            log("  %-9s markers (n %d): own axis %+.4f, other axis %+.4f"
                % (kl, len(sub), sub[sp_eigen].mean(), sub[sp_fremd].mean()))

    # ------------------------------------------------- sensitivities
    sym_module = dict(zip(MODUL.symbol, MODUL.ri))
    zeilen = []

    def lauf(name: str, sondenmenge: pd.Index, delta: pd.Series, achse: str) -> None:
        mm = man[man.sonde.isin(sondenmenge)]
        gen = mm.assign(d=mm.sonde.map(delta)).groupby("symbol").d.mean().dropna()
        mod = gen[gen.index.isin(sym_module)]
        if len(mod) < 8:
            log("  %-38s too few module genes (%d)" % (name, len(mod)))
            return
        erw = pd.Series({s: -sym_module[s] for s in mod.index})
        r = konkordanz(mod, erw, hintergrund=gen)
        log("  %-38s n %3d  C %.3f (Null %.3f+-%.3f)  z %+5.2f  p %.3g"
            % (name, r["n"], r["konkordanz"], r["konkordanz_null"],
               r["konkordanz_null_sd"], r["konkordanz_z"], r["konkordanz_p"]))
        zeilen.append(dict(variante=name, achse=achse,
                           **{k: v for k, v in r.items() if k != "status"}))

    for achse, delta in (("osteo", d_ost), ("myo", d_myo)):
        log("\n--- Sensitivities, axis %s ---------------------------------" % achse)
        alle = pd.Index(man.sonde)
        lauf("alle promotornahen Sonden", alle, delta, achse)
        lauf("S1 nur Sonden mit Dynamik",
             pd.Index(man.sonde[dynamisch.reindex(man.sonde).fillna(False).values]),
             delta, achse)
        lauf("S2a nur CpG-Insel",
             pd.Index(man.sonde[insel.reindex(man.sonde).fillna(False).values]),
             delta, achse)
        lauf("S2b nur ausserhalb der Insel",
             pd.Index(man.sonde[~insel.reindex(man.sonde).fillna(True).values]),
             delta, achse)

    # ------------------------------- S3 drift-controlled axis difference
    log("\n--- S3 drift-controlled: osteogenic minus myogenic axis -----------")
    log("Question: does the module move more strongly in the expected")
    log("direction in the osteogenic than in the myogenic axis, beyond global drift?")
    dd = d_ost - d_myo
    gen_dd = man.assign(d=man.sonde.map(dd)).groupby("symbol").d.mean().dropna()
    mod_dd = gen_dd[gen_dd.index.isin(sym_module)]
    erw = pd.Series({s: -sym_module[s] for s in mod_dd.index})
    r = konkordanz(mod_dd, erw, hintergrund=gen_dd)
    log("n %d | concordance %.3f against null %.3f +- %.3f | z %+.2f | p %.4g"
        % (r["n"], r["konkordanz"], r["konkordanz_null"], r["konkordanz_null_sd"],
           r["konkordanz_z"], r["konkordanz_p"]))
    log("continuous: %+.4f, z %+.2f, p %.4g" % (r["rang"], r["rang_z"], r["rang_p"]))
    log("MDE80 threshold %.3f -- reached: %s"
        % (r["konkordanz_mde80"],
           "YES" if r["konkordanz"] >= r["konkordanz_mde80"] else "no"))
    zeilen.append(dict(variante="S3 osteo minus myo (driftkontrolliert)",
                       achse="differenz",
                       **{k: v for k, v in r.items() if k != "status"}))

    pd.DataFrame(zeilen).to_csv(AUS / "A_dnam_sensitivitaeten.csv", index=False)
    (AUS / "A_dnam_sensitivitaet_log.txt").write_text("\n".join(LOG) + "\n",
                                                      encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

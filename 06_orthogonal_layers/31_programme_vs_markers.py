# -*- coding: utf-8 -*-
"""
31_programme_vs_markers.py -- the direct competition on foreign ground.

Observation from 26_ and 30_: on the replicated chromatin layer
(GSE224251) the data-driven 173-gene module moves directionally and robust
to baseline, while the canonical osteogenic marker set does not. This is
not a side note -- it is a testable claim, and it is tested directly here.

Design: on each orthogonal layer the same directed statistic is computed
for three gene sets, against the same null and with the same weights:

    MODULE     the 173 convergent genes, direction `ri`
    MARKER     the canonical lineage marker set, direction "up"
    RANDOM     sets of equal size from the background (the null itself)

The question is not which set is "better", but: which of the two is
recoverable on a measurement layer that neither of them produced. A marker
set stems from decades of single-gene work, the module from eighteen
transcriptome datasets. If only the module carries, that says something
about how differentiation is organized at the chromatin level -- namely
distributed and not in a few lead genes.

Additionally the size question is settled: marker sets are small (10-18
genes), the module is large (173). To make the two comparable, the module
is drawn DOWN to marker size -- 2 000 random subsets of the module with
size equal to the marker set. Reported is the share of these subsets that
beats the marker set.
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
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, OSTEOGEN  # noqa: E402
from _module import ERGEBNISSE, MODUL, SEED  # noqa: E402

import importlib.util  # noqa: E402
spec = importlib.util.spec_from_file_location(
    "i30", pathlib.Path(__file__).resolve().parent / "30_convergence_dose_integration.py")
i30 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(i30)

AUS = ERGEBNISSE / "Z_integration"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []
NZIEH = 20000
NTEIL = 2000

# Which marker set belongs to which layer, and in which direction?
# Positive layer difference = "matches a gene running up" (see 30_).
MARKERZUORDNUNG = {
    "ATAC GSE224251 osteogen (n=3)": ("osteogene Marker", OSTEOGEN),
    "ATAC GSE332758 osteogen (n=1)": ("osteogene Marker", OSTEOGEN),
    "ATAC GSE332758 adipogen (n=1)": ("adipogene Marker", ADIPOGEN),
    "DNAm GSE33896 osteogen (n=3)": ("osteogene Marker", OSTEOGEN),
    "DNAm GSE33896 myogen (n=3)": ("myogene Marker", MYOGEN),
    "DNAm GSE129266 chondrogen (n=2)": ("chondrogene Marker", CHONDROGEN),
    "H3K27ac GSE129031 chondrogen (n=2)": ("chondrogene Marker", CHONDROGEN),
}


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def gerichtet(delta: pd.Series, gene, richtung, seed: int = SEED) -> dict:
    """Mean signed rank of a gene set against the background draw of equal
    size and equal sign set."""
    d = delta.dropna()
    r = pd.Series(stats.rankdata(d.values) / len(d) - 0.5, index=d.index)
    g = [x for x in gene if x in r.index]
    if len(g) < 5:
        return {"n": len(g), "status": "too few"}
    vz = np.array([richtung[x] if isinstance(richtung, dict) else richtung
                   for x in g], dtype=float)
    beob = float((vz * r.reindex(g).values).mean())
    rng = np.random.default_rng(seed)
    rv = r.values
    null = np.empty(NZIEH)
    for i in range(NZIEH):
        idx = rng.choice(len(rv), size=len(g), replace=False)
        null[i] = (rng.permutation(vz) * rv[idx]).mean()
    sd = null.std(ddof=1)
    return {"n": len(g), "statistik": beob, "null": float(null.mean()),
            "null_sd": float(sd), "z": float((beob - null.mean()) / sd),
            "p": float(min(1.0, 2 * (1 + (null >= beob).sum()) / (1 + len(null)))),
            "status": "ok"}


def main() -> None:
    log("=" * 78)
    log("Module against marker panel -- same layer, same null, same statistic")
    log("=" * 78)

    E = i30.lade_ebenen()
    ri = dict(zip(MODUL.symbol, MODUL.ri))
    zeilen = []

    log("\n%-34s %-20s %5s %9s %8s %9s"
        % ("layer", "gene set", "n", "statistic", "z", "p"))
    for ename, delta in E.items():
        mname, mgene = MARKERZUORDNUNG.get(ename, (None, None))
        r_mod = gerichtet(delta, list(ri), ri)
        r_mar = gerichtet(delta, mgene, +1) if mgene else {"status": "none"}
        for satz, nm, r in (("Modul", "173 konvergente Gene", r_mod),
                            ("Marker", mname or "-", r_mar)):
            if r.get("status") != "ok":
                log("%-34s %-20s %s" % (ename, nm, r.get("status")))
                continue
            log("%-34s %-20s %5d %+9.4f %+8.2f %9.4g"
                % (ename, nm, r["n"], r["statistik"], r["z"], r["p"]))
            zeilen.append(dict(ebene=ename, satzart=satz, satz=nm,
                               **{k: v for k, v in r.items() if k != "status"}))

        # ------------------------- size control: module at marker size
        if r_mar.get("status") == "ok" and r_mod.get("status") == "ok":
            d = delta.dropna()
            rr = pd.Series(stats.rankdata(d.values) / len(d) - 0.5, index=d.index)
            mg = [x for x in ri if x in rr.index]
            k = r_mar["n"]
            rng = np.random.default_rng(SEED + 1)
            werte = np.empty(NTEIL)
            for i in range(NTEIL):
                aus = rng.choice(len(mg), size=min(k, len(mg)), replace=False)
                sub = [mg[j] for j in aus]
                vz = np.array([ri[x] for x in sub], dtype=float)
                werte[i] = (vz * rr.reindex(sub).values).mean()
            anteil = float((werte > r_mar["statistik"]).mean())
            log("%-34s %-20s module subsets of size %d beat the "
                "marker set in %.1f %% of cases (median %+.4f vs %+.4f)"
                % ("", "-> size control", k, 100 * anteil,
                   float(np.median(werte)), r_mar["statistik"]))
            zeilen.append(dict(ebene=ename, satzart="Groessenkontrolle",
                               satz="Modulteilmengen n=%d" % k, n=k,
                               statistik=float(np.median(werte)),
                               null=r_mar["statistik"], null_sd=np.nan,
                               z=np.nan, p=1 - anteil))
        log("")

    Z = pd.DataFrame(zeilen)
    Z.to_csv(AUS / "Z_modul_gegen_marker.csv", index=False)

    log("=" * 78)
    log("SUMMARY")
    log("=" * 78)
    log("%-34s %10s %10s %12s" % ("layer", "z(module)", "z(marker)", "difference"))
    for ename in Z.ebene.unique():
        a = Z[(Z.ebene == ename) & (Z.satzart == "Modul")]
        b = Z[(Z.ebene == ename) & (Z.satzart == "Marker")]
        if len(a) and len(b):
            log("%-34s %+10.2f %+10.2f %+12.2f"
                % (ename, a.z.iloc[0], b.z.iloc[0], a.z.iloc[0] - b.z.iloc[0]))
    log("")
    log("A module that carries on a foreign layer, where the canonical marker")
    log("set of the same lineage does not, says something about the")
    log("organization of differentiation: it is distributed and not confined")
    log("to the few genes the field reads it off.")

    (AUS / "Z_modul_gegen_marker_log.txt").write_text("\n".join(LOG) + "\n",
                                                      encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

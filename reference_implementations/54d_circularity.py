# -*- coding: utf-8 -*-
"""
54d_circularity.py -- the circularity control of phase M-D.

**POST HOC, not preregistered.** This stands in the first line because it
must also stand in the first line of the report.

The occasion. The fixed 173-gene module was formed from **eighteen**
perturbation experiments. Four of the phase-D studies belong to these
eighteen: GSE218101, GSE221128, GSE245585 and the own LAMA5-USC series.
Their `dWT` vectors co-determined the module -- high agreement with `ri` is
partly built in there. Only the SERPINA3 series (GSE247491/GSE247528) did
**not** participate in module formation; they have, however, already been
analyzed per study as Fig. S3C (never donor-resolved).

The same ladder is therefore computed once more on two subsets:

  (a) only cells that did **not** co-define the module
  (b) only cells that did co-define it

The statistic S3a (directed concordance against `ri`) is affected most by
circularity, S1 and S2 least -- they ask about the agreement of donors
**with each other** and do not use `ri` at all. Exactly why S1 is the
primary statistic.

Output: derived_data/M_donoren/zirkularitaet.csv
"""
from __future__ import annotations

import pathlib
import pickle
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import ERGEBNISSE, MODUL, NZIEHUNGEN, SEED_D, leiter  # noqa: E402

AUS = ERGEBNISSE / "M_donoren"

# The eighteen points from which the module was formed -- as far as they
# reappear in phase D. Source: reference_implementations/manuscript/methods/03_metric.R.
MODULBILDEND = {"GSE218101", "GSE221128", "GSE245585", "LAMA5_USC"}


def main() -> None:
    with open(AUS / "zellen.pkl", "rb") as f:
        daten = pickle.load(f)
    E = pd.read_csv(AUS / "eichung.csv")
    g = E[E.bestanden]

    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
    from importlib import import_module
    lc = import_module("54c_leiter") if False else None   # noqa: F841

    def sammle(zellen):
        teile, basen = [], []
        for gse, d in daten.items():
            sp = [c for c in d["dwt"].columns if c in zellen]
            if sp:
                teile.append(d["dwt"][sp].T)
                basen.append(d["basis"])
        X = pd.concat(teile, axis=0, join="inner")
        basis = pd.concat(basen, axis=1).mean(axis=1).reindex(X.columns).dropna()
        X = X[basis.index]
        X = X.loc[:, X.notna().all(axis=0)]
        return X, basis.reindex(X.columns)

    ri = MODUL.set_index("symbol").ri
    zeilen = []
    teilmengen = {
        "alle geeichten Zellen": list(g.zelle),
        "nur NICHT modulbildend": list(g[~g.studie.isin(MODULBILDEND)].zelle),
        "nur modulbildend": list(g[g.studie.isin(MODULBILDEND)].zelle),
    }
    print("=" * 78)
    print("Circularity control phase M-D  --  POST HOC, "
          "not preregistered")
    print("=" * 78)
    for name, zellen in teilmengen.items():
        if len(zellen) < 2:
            print("\n%-26s %d cell(s) -- no statistic possible"
                  % (name, len(zellen)))
            zeilen.append(dict(teilmenge=name, n_zellen=len(zellen),
                               status="too few cells"))
            continue
        X, basis = sammle(zellen)
        sp = list(g.set_index("zelle").spender.reindex(X.index))
        gene = [x for x in MODUL.symbol if x in X.columns]
        r = leiter(X, gene, ri, X, basis, spender=sp,
                   nziehungen=NZIEHUNGEN, seed=SEED_D)
        r.pop("je_zelle", None)
        r.update(teilmenge=name, n_studien=len(set(
            g.set_index("zelle").studie.reindex(X.index))))
        zeilen.append(r)
        print("\n%-26s %d cells, %d donors, %d studies"
              % (name, r["n_zellen"], len(set(sp)), r["n_studien"]))
        for s in ("S1", "S2", "S3a", "S3b"):
            print("   %-4s %+7.4f | null %+7.4f +- %6.4f | z %+6.2f | "
                  "p %8.4g | MDE80 %+7.4f %s"
                  % (s, r[f"{s}_beobachtet"], r[f"{s}_null_mittel"],
                     r[f"{s}_null_sd"], r[f"{s}_z"], r[f"{s}_p"],
                     r[f"{s}_mde80"],
                     "ABOVE" if r[f"{s}_beobachtet"] > r[f"{s}_mde80"] else ""))
    pd.DataFrame(zeilen).to_csv(AUS / "zirkularitaet.csv", index=False)
    print("\n->", AUS / "zirkularitaet.csv")


if __name__ == "__main__":
    main()

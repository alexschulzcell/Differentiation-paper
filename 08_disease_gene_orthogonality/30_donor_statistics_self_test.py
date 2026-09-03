# -*- coding: utf-8 -*-
"""
30_donor_statistics_self_test.py -- self-test of the statistical ladder from `_module.leiter`.

Preregistration `PRAEREG_M_D.md` §12: before the first real computation the
new statistics are checked against the known null rate. Two questions, and
only these two:

  (1) **Is the null calibrated?**  If the cell vectors are drawn THEMSELVES
      from the background -- i.e. without any shared structure -- the rate
      of cases with z > 2 must be around 2.3 % (one-sided, normal
      approximation), and the rate |z| > 2 around 5 %. If it is clearly
      higher, the null is too narrow and every finding of this phase would
      be an artifact.

  (2) **Does the ladder find what it should find?**  If a shared program of
      known strength is implanted, S1, S2 and S3 must see it.

The script computes NOTHING on real data and changes no result file except
its own.

Output: derived_data/M_donoren/selbsttest.csv
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
from _module import ERGEBNISSE, SEED_D, leiter  # noqa: E402

AUS = ERGEBNISSE / "M_donoren"
AUS.mkdir(parents=True, exist_ok=True)

N_ZELLEN = 6          # as the order of magnitude of phase D
N_GENE = 170          # as the measurable module genes
N_HG = 4000           # background genes
N_ZIEH = 500          # draws per computation (self-test, not the main run)
N_WDH = 200           # repetitions


def eine_runde(rng: np.random.Generator, effekt: float) -> dict:
    """One artificial study: background + optional shared program."""
    # baseline expression per gene, so that stratification has something to do
    basis = pd.Series(rng.normal(0, 1.5, N_HG),
                      index=[f"G{i}" for i in range(N_HG)])
    # background: cell x gene, with gene-dependent spread (like real data)
    skala = np.exp(rng.normal(0, 0.4, N_HG))
    HG = pd.DataFrame(rng.normal(0, 1, (N_ZELLEN, N_HG)) * skala,
                      columns=basis.index,
                      index=[f"Z{i}" for i in range(N_ZELLEN)])
    gene = list(basis.index[:N_GENE])
    ri = pd.Series(rng.choice([-1, 1], N_HG), index=basis.index)
    if effekt > 0:
        # shared program: the same directed offset in every cell
        HG.loc[:, gene] = HG.loc[:, gene].to_numpy() + effekt * ri[gene].to_numpy()[None, :]
    X = HG[gene]
    return leiter(X, gene, ri, HG, basis, nziehungen=N_ZIEH,
                  seed=int(rng.integers(1, 2**31)))


def main() -> None:
    print("=" * 78)
    print("Self-test of the statistical ladder  --  seed", SEED_D)
    print("=" * 78)
    zeilen = []
    for effekt, wdh in ((0.0, N_WDH), (0.35, 40)):
        rng = np.random.default_rng(SEED_D + int(effekt * 1000))
        werte = {s: [] for s in ("S1", "S2", "S3a", "S3b")}
        for _ in range(wdh):
            r = eine_runde(rng, effekt)
            for s in werte:
                werte[s].append(r[f"{s}_z"])
        for s, v in werte.items():
            v = np.array(v, dtype=float)
            v = v[np.isfinite(v)]
            zeilen.append({
                "effekt": effekt, "kennzahl": s, "n_wiederholungen": len(v),
                "mittlerer_z": float(v.mean()), "sd_z": float(v.std(ddof=1)),
                "anteil_z_ueber_2": float((v > 2).mean()),
                "anteil_abs_z_ueber_2": float((np.abs(v) > 2).mean()),
            })
            print("Effect %.2f | %-4s | mean z %+6.2f | sd %.2f | "
                  "z>2 %5.1f %% | |z|>2 %5.1f %%"
                  % (effekt, s, v.mean(), v.std(ddof=1),
                     100 * (v > 2).mean(), 100 * (np.abs(v) > 2).mean()))
    T = pd.DataFrame(zeilen)
    T.to_csv(AUS / "selbsttest.csv", index=False)
    print("\n->", AUS / "selbsttest.csv")
    print("Expectation: at effect 0.00, z>2 is around 2.3 % and |z|>2 around "
          "5 %. At effect 0.35 all four statistics must fire.")


if __name__ == "__main__":
    main()

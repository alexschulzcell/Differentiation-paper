# -*- coding: utf-8 -*-
"""
_hardening.py -- shared foundation for the held-out validation and robustness
analyses of Figure 2G-I.

Reproduces the frozen 173-gene module-derivation rule
(reference_implementations/manuscript/methods/20f_convergence_dwt.py) in a form
that can be applied to ANY subset of the 18 datasets, and provides the
project's concordance test (the concordance branch of
reference_implementations/_module.py `konkordanz`, identical statistic).

Everything runs in Ensembl-ID space. Reads only repository files:
  derived_data/reference_tables/20d_dWT_matrix.csv.gz   long: gen, punkt, dWT
  derived_data/reference_tables/dataset_study_map.csv   point -> gse study map
  derived_data/reference_tables/S5_konvergente_gene.csv the frozen 173 genes
  derived_data/followup/ws1_genkarte_erweitert.csv      per-gene covariates
"""
from __future__ import annotations
import os
import pathlib
import numpy as np
import pandas as pd

SEED = 20260821

_env = os.environ.get("PAPER_V2_ROOT")
ROOT = (pathlib.Path(_env) if _env
        else pathlib.Path(__file__).resolve().parents[1])
TAB = ROOT / "derived_data" / "reference_tables"


def paths() -> dict:
    return {
        "dwt_matrix": TAB / "20d_dWT_matrix.csv.gz",
        "aufteilung": TAB / "dataset_study_map.csv",
        "s5": TAB / "S5_konvergente_gene.csv",
        "genkarte": ROOT / "derived_data" / "followup" / "ws1_genkarte_erweitert.csv",
    }


def load_dwt_wide() -> pd.DataFrame:
    """Genes x 18 points matrix of dWT (index = Ensembl gen)."""
    long = pd.read_csv(paths()["dwt_matrix"], compression="gzip", dtype={"gen": str})
    wide = long.pivot(index="gen", columns="punkt", values="dWT")
    return wide.reindex(columns=sorted(wide.columns))


def study_map() -> pd.DataFrame:
    """One row per point: punkt, datensatz, gse, arm, klasse."""
    a = pd.read_csv(paths()["aufteilung"], dtype={"punkt": int})
    cols = ["punkt", "datensatz", "gse", "arm", "klasse"]
    return a[cols].drop_duplicates("punkt").sort_values("punkt").reset_index(drop=True)


def derive_module(wide: pd.DataFrame, points: list[int],
                  min_frac: float = 14 / 18, vn: float = 0.90) -> pd.DataFrame:
    """Apply the frozen convergence rule to an arbitrary subset of points.

    Returns DataFrame index=gen with columns n, v, ri (ri in {-1,+1}).
    min_frac scales the >=14/18 universe threshold to the subset size.
    """
    sub = wide[points]
    min_da = int(np.ceil(min_frac * len(points)))
    available = sub.notna().sum(axis=1)
    universe = sub.loc[available >= min_da]
    medians = universe.median(axis=0)
    signs = np.sign(universe.subtract(medians, axis="columns")).fillna(0)
    positive = (signs > 0).sum(axis=1)
    n = (signs != 0).sum(axis=1)
    v = np.maximum(positive, n - positive)
    ri = np.where(positive >= n - positive, 1, -1)
    out = pd.DataFrame({"n": n.astype(int), "v": v.astype(int), "ri": ri.astype(int)},
                       index=universe.index)
    out = out.loc[out["n"] > 0]
    out = out.loc[out["v"] / out["n"] >= vn]
    return out


def concordance(delta: pd.Series, ri: pd.Series, background: pd.Series,
                nziehungen: int = 5000, seed: int = SEED) -> dict:
    """Direction test: share of set genes with sign(delta)==ri versus a
    size-matched, sign-preserving background draw. Identical logic to
    reference_implementations/_module.py konkordanz (concordance branch)."""
    d = pd.Series(delta).dropna()
    e = pd.Series(ri).reindex(d.index)
    ok = e.notna() & (d != 0)
    d, e = d[ok], e[ok].astype(int)
    k = len(d)
    if k < 8:
        return {"n": k, "status": "too few"}
    hg = pd.Series(background).dropna()
    hg = hg[hg != 0]
    beob = float((np.sign(d.values) == e.values).mean())
    rng = np.random.default_rng(seed)
    hgv, ev = hg.values, e.values
    nc = np.empty(nziehungen)
    for i in range(nziehungen):
        idx = rng.choice(len(hgv), size=k, replace=False)
        nc[i] = (np.sign(hgv[idx]) == rng.permutation(ev)).mean()
    m, s = float(nc.mean()), float(nc.std(ddof=1))
    z = (beob - m) / s if s > 0 else np.nan
    p = (1 + int((nc >= beob).sum())) / (1 + len(nc))
    return {"n": k, "concordance": beob, "null_mean": m, "null_sd": s,
            "z": float(z), "p": float(min(1.0, 2 * p)),
            "mde80": m + 2.8 * s, "above_mde80": bool(beob > m + 2.8 * s),
            "status": "ok"}

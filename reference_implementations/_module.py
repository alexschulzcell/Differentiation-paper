# -*- coding: utf-8 -*-
"""
_module.py -- shared foundation of all orthogonal layers.

Provides:
  * MODUL   : the fixed 173 convergent genes with their direction `ri`
              (from derived_data/reference_tables/S5_konvergente_gene.csv; not recomputed)
  * konkordanz(...) : the direction test of a module on any orthogonal
              measurement layer, against a background-drawn null.

Logic of the direction test
---------------------------
On each orthogonal layer there is per gene a difference `delta_g`
(differentiated minus naive). The module predicts a sign `s_g` for every
gene. What is tested is the concordance

    C = share of module genes with sign(delta_g) == s_g .

The null draws gene sets of equal size from the genes measurable on this
layer and assigns them **the same set of signs** that the module carries.
The null is therefore robust against two plausible artifacts:

  (1) global drift of the layer (e.g. global hypomethylation during
      differentiation) -- it affects module and null alike;
  (2) the sign imbalance of the module itself (129 up, 44 down) -- it is
      carried along rather than averaged out.

In addition a continuous variant is reported: the mean signed rank of
`delta_g`. It is insensitive to the scale of the respective layer and to
single outliers.

There is exactly one implementation of this test in this project, and it is
this one.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
from scipy import stats

WURZEL = pathlib.Path(__file__).resolve().parents[1]
TABELLEN = WURZEL / "derived_data" / "reference_tables"
ERGEBNISSE = WURZEL / "derived_data"
DATEN = WURZEL / "data_raw"

SEED = 20260821          # fixed; the same root for all orthogonal layers
NZIEHUNGEN = 20000


def lade_module() -> pd.DataFrame:
    """The 173 convergent genes with Ensembl ID, symbol and direction `ri`."""
    m = pd.read_csv(TABELLEN / "S5_konvergente_gene.csv")
    m = m.rename(columns={"gen": "ensembl"})
    m["symbol"] = m["symbol"].astype(str).str.strip()
    assert len(m) == 173, "The module is fixed: 173 genes."
    assert set(m.ri.unique()) <= {-1, 1}
    return m[["ensembl", "symbol", "ri", "n", "v", "med"]]


MODUL = lade_module()


def konkordanz(delta: pd.Series,
               erwartet: pd.Series,
               hintergrund: pd.Series | None = None,
               schichtung: pd.Series | None = None,
               nziehungen: int = NZIEHUNGEN,
               seed: int = SEED) -> dict:
    """Direction test of a module on an orthogonal layer.

    Parameters
    ----------
    delta        : series `gene -> delta` for the module genes measurable on
                   this layer (index = gene key).
    erwartet     : series `gene -> +1/-1`, same index set as `delta`.
    hintergrund  : series `gene -> delta` over ALL genes measurable on this
                   layer. Default: `delta` itself (then the null is only the
                   sign null, which is weaker).
    """
    d = pd.Series(delta).dropna()
    e = pd.Series(erwartet).reindex(d.index)
    ok = e.notna() & (d != 0)
    d, e = d[ok], e[ok].astype(int)
    k = len(d)
    if k < 8:
        return {"n": k, "status": "too few measurable genes"}

    hg = pd.Series(hintergrund if hintergrund is not None else d).dropna()
    hg = hg[hg != 0]

    beob_c = float((np.sign(d.values) == e.values).mean())

    # continuous variant: signed rank within the background
    rang_hg = pd.Series(stats.rankdata(hg.values) / len(hg), index=hg.index)
    gem = d.index.intersection(rang_hg.index)
    beob_r = float((e.reindex(gem).values * (rang_hg.reindex(gem).values - 0.5)).mean())

    rng = np.random.default_rng(seed)
    hgv = hg.values
    hgr = rang_hg.values
    nc = np.empty(nziehungen)
    nr = np.empty(nziehungen)
    ev = e.values

    # Stratified draw: `schichtung` is a nuisance variable per gene (e.g.
    # the baseline accessibility in the naive state). The null then draws for
    # each module gene from THE SAME decile stratum. This rules out the
    # result arising because module genes start out systematically different.
    if schichtung is not None:
        sch = pd.Series(schichtung).reindex(hg.index).dropna()
        hg2 = hg.reindex(sch.index)
        hgv, hgr = hg2.values, rang_hg.reindex(sch.index).values
        dez = pd.qcut(sch.rank(method="first"), 10, labels=False).values
        pos = {b: np.flatnonzero(dez == b) for b in range(10)}
        d_mod = pd.Series(schichtung).reindex(d.index)
        # assign module genes to their strata (via the stratum boundaries)
        raender = np.quantile(sch.values, np.linspace(0, 1, 11)[1:-1])
        b_mod = np.searchsorted(raender, d_mod.values)
        b_mod = np.where(np.isnan(d_mod.values), -1, b_mod)
        gueltig = b_mod >= 0
        ev_g = ev[gueltig]
        b_g = b_mod[gueltig]
        for i in range(nziehungen):
            idx = np.array([rng.choice(pos[b]) for b in b_g])
            vz = rng.permutation(ev_g)
            nc[i] = (np.sign(hgv[idx]) == vz).mean()
            nr[i] = (vz * (hgr[idx] - 0.5)).mean()
    else:
        for i in range(nziehungen):
            idx = rng.choice(len(hgv), size=k, replace=False)
            vz = rng.permutation(ev)
            nc[i] = (np.sign(hgv[idx]) == vz).mean()
            nr[i] = (vz * (hgr[idx] - 0.5)).mean()

    def z_und_p(beob, null):
        sd = null.std(ddof=1)
        z = (beob - null.mean()) / sd if sd > 0 else np.nan
        # empirical, one-sided upper, with edge correction
        p_ein = (1 + int((null >= beob).sum())) / (1 + len(null))
        return float(z), float(min(1.0, 2 * p_ein)), float(null.mean()), float(sd)

    z_c, p_c, m_c, s_c = z_und_p(beob_c, nc)
    z_r, p_r, m_r, s_r = z_und_p(beob_r, nr)

    return {
        "n": k,
        "konkordanz": beob_c,
        "konkordanz_null": m_c,
        "konkordanz_null_sd": s_c,
        "konkordanz_z": z_c,
        "konkordanz_p": p_c,
        "konkordanz_mde80": m_c + 2.8 * s_c,     # project threshold
        "rang": beob_r,
        "rang_null": m_r,
        "rang_null_sd": s_r,
        "rang_z": z_r,
        "rang_p": p_r,
        "status": "ok",
    }


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson interval for a proportion."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    mitte = (p + z * z / (2 * n)) / (1 + z * z / n)
    halb = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, mitte - halb), min(1.0, mitte + halb))


if __name__ == "__main__":
    print(MODUL.head())
    print("Module:", len(MODUL), "genes |", int((MODUL.ri > 0).sum()), "up,",
          int((MODUL.ri < 0).sum()), "down")


# ---------------------------------------------------------------------------
# The two-set contrast -- the same statistic as in the paper core,
# transferred to an orthogonal layer.
# ---------------------------------------------------------------------------
def kontrast(werte: "pd.Series", satz_a, satz_b,
             nziehungen: int = NZIEHUNGEN, seed: int = SEED) -> dict:
    """Difference of two gene sets against the same background-drawn null.

    `werte` is a series `gene -> value` over the entire measurable
    background. What is tested is mean(A) - mean(B) against a null that
    draws two random sets of the same sizes from the same background.
    This statistic is insensitive to any global offset of the layer --
    exactly the argument that saves the per-gene z scale in the paper core.
    """
    w = pd.Series(werte).dropna()
    a = [g for g in satz_a if g in w.index]
    b = [g for g in satz_b if g in w.index]
    if len(a) < 3 or len(b) < 3:
        return {"n_a": len(a), "n_b": len(b), "status": "too few genes"}
    beob = float(w[a].mean() - w[b].mean())
    rng = np.random.default_rng(seed)
    v = w.values
    null = np.empty(nziehungen)
    for i in range(nziehungen):
        idx = rng.choice(len(v), size=len(a) + len(b), replace=False)
        null[i] = v[idx[:len(a)]].mean() - v[idx[len(a):]].mean()
    sd = null.std(ddof=1)
    z = (beob - null.mean()) / sd if sd > 0 else np.nan
    p_ein = (1 + int((null >= beob).sum())) / (1 + len(null))
    return {"n_a": len(a), "n_b": len(b), "kontrast": beob,
            "null_mittel": float(null.mean()), "null_sd": float(sd),
            "z": float(z), "p": float(min(1.0, 2 * p_ein)),
            "mde80": float(null.mean() + 2.8 * sd), "status": "ok"}


# ===========================================================================
# Phase M-D -- the statistical ladder of the donor-resolved shear.
#
# Preregistration: `preregistrations/PRAEREG_M_D.md` §7/§8, dated
# 2026-08-22 before the first download. These functions are an EXTENSION of
# that file, not a second implementation: the null is the same
# baseline-stratified background draw as in `konkordanz`; only the statistic
# computed on the drawn set differs.
#
# Data situation. Each study provides a matrix X (rows = cells/donors,
# columns = genes). The rows are either the per-donor `dWT` vectors
# (program) or the per-donor `iv` vectors (lesion response), formed as in
# `03_metric.R`.
#
# The four statistics (all three are reported; S1 is primary):
#
#   S1  mean pairwise Spearman correlation between different donors,
#       restricted to the set genes. Per cell additionally the mean over
#       its partners -- this is the quantity that gets its own MDE80.
#   S2  share of variance carried by the first principal component.
#       IMPORTANT: **not** centered over donors. Exactly the shared share
#       that centering over donors would remove is the quantity of interest
#       (addendum 1 of the preregistration).
#   S3a directed sign concordance against `ri` -- the statistic from M-B.
#   S3b direction-free agreement: mean_g |mean_p sign(x_pg)|.
#
# The null draws gene sets of equal size from the same deciles of baseline
# expression and computes the same four statistics on them.
# ===========================================================================

SEED_D = 20260822


def _spearman_paare(M: np.ndarray) -> np.ndarray:
    """All pairwise Spearman rhos of the rows of M (cells x genes)."""
    R = np.apply_along_axis(stats.rankdata, 1, M)
    R = R - R.mean(axis=1, keepdims=True)
    n = np.sqrt((R * R).sum(axis=1))
    n[n == 0] = np.nan
    return (R @ R.T) / np.outer(n, n)


def _pc1_anteil(M: np.ndarray) -> float:
    """Share of the first principal component in the total sum of squares.

    Without centering over donors -- see header comment.
    """
    if M.shape[0] < 2:
        return np.nan
    s = np.linalg.svd(M, compute_uv=False)
    ges = float((s ** 2).sum())
    return float(s[0] ** 2 / ges) if ges > 0 else np.nan


def _kennzahlen(M: np.ndarray, ri: np.ndarray,
                paare: np.ndarray) -> dict:
    """S1, S2, S3a, S3b on a matrix M (cells x genes) with direction `ri`.

    `paare` is a boolean mask (cells x cells): which pairs are admissible,
    i.e. from different donors? The diagonal is always off.
    """
    n = M.shape[0]
    rho = _spearman_paare(M)
    ob = np.triu(paare, 1)
    with np.errstate(invalid="ignore"):
        s1_je_zelle = np.array([np.nanmean(rho[i][paare[i]])
                                if paare[i].any() else np.nan
                                for i in range(n)])
        s1 = float(np.nanmean(rho[ob])) if ob.any() else np.nan
    vz = np.sign(M)
    s3a_je_zelle = (vz == ri[None, :]).mean(axis=1)
    return {
        "S1": s1,
        "S1_je_zelle": s1_je_zelle,
        "S2": _pc1_anteil(M),
        "S3a": float(s3a_je_zelle.mean()),
        "S3a_je_zelle": s3a_je_zelle,
        "S3b": float(np.abs(vz.mean(axis=0)).mean()),
    }


def leiter(X: "pd.DataFrame",
           gene: list,
           ri: "pd.Series",
           hintergrund: "pd.DataFrame",
           basis: "pd.Series",
           spender: list | None = None,
           nziehungen: int = NZIEHUNGEN,
           seed: int = SEED_D) -> dict:
    """The preregistered statistical ladder with baseline-stratified null.

    Parameters
    ----------
    X            DataFrame cell x gene -- the per-donor vectors.
    gene         the set genes (their order determines `ri`).
    ri           series gene -> +1/-1, the direction predicted by the module.
    hintergrund  DataFrame cell x gene over ALL measurable genes, same rows
                 as X. The null is drawn from it.
    basis        series gene -> baseline expression (naive control arm),
                 averaged over cells. Stratification variable of the null.
    spender      donor identifier per row. Pairs of the same donor are
                 excluded. Default: every row is its own donor.

    Returns: observed statistics, null mean, null SD, z, p and MDE80 per
    statistic -- and the same per cell for S1 and S3a.
    """
    gene = [g for g in gene if g in X.columns]
    if len(gene) < 8 or X.shape[0] < 2:
        return {"status": "too few genes or cells",
                "n_gene": len(gene), "n_zellen": int(X.shape[0])}
    M = np.nan_to_num(X[gene].to_numpy(dtype=float), nan=0.0)
    ev = ri.reindex(gene).to_numpy(dtype=float)

    sp = np.asarray(list(spender) if spender is not None else list(X.index))
    paare = (sp[:, None] != sp[None, :])

    beob = _kennzahlen(M, ev, paare)

    # ---- baseline-stratified background draw ------------------------------
    H = hintergrund.reindex(index=X.index)
    b = pd.Series(basis).reindex(H.columns).dropna()
    HM = np.nan_to_num(H[b.index].to_numpy(dtype=float), nan=0.0)
    dez = pd.qcut(b.rank(method="first"), 10, labels=False).to_numpy()
    pos = {d: np.flatnonzero(dez == d) for d in range(10)}
    raender = np.quantile(b.to_numpy(), np.linspace(0, 1, 11)[1:-1])
    b_set = np.searchsorted(raender, pd.Series(basis).reindex(gene).to_numpy())

    rng = np.random.default_rng(seed)
    null = {s: np.empty(nziehungen) for s in ("S1", "S2", "S3a", "S3b")}
    null_s1_zelle = np.empty((nziehungen, M.shape[0]))
    null_s3a_zelle = np.empty((nziehungen, M.shape[0]))
    for i in range(nziehungen):
        idx = np.array([rng.choice(pos[d]) for d in b_set])
        vz = rng.permutation(ev)
        r = _kennzahlen(HM[:, idx], vz, paare)
        for s in null:
            null[s][i] = r[s]
        null_s1_zelle[i] = r["S1_je_zelle"]
        null_s3a_zelle[i] = r["S3a_je_zelle"]

    def fasse(beobachtet, nullwerte):
        nw = np.asarray(nullwerte, dtype=float)
        nw = nw[np.isfinite(nw)]
        if not np.isfinite(beobachtet) or len(nw) < 100:
            return {"beobachtet": float(beobachtet), "null_mittel": np.nan,
                    "null_sd": np.nan, "z": np.nan, "p": np.nan,
                    "mde80": np.nan}
        m, s = float(nw.mean()), float(nw.std(ddof=1))
        z = (beobachtet - m) / s if s > 0 else np.nan
        p_ein = (1 + int((nw >= beobachtet).sum())) / (1 + len(nw))
        return {"beobachtet": float(beobachtet), "null_mittel": m,
                "null_sd": s, "z": float(z),
                "p": float(min(1.0, 2 * p_ein)), "mde80": m + 2.8 * s}

    aus = {"status": "ok", "n_gene": len(gene), "n_zellen": int(M.shape[0]),
           "n_paare": int(np.triu(paare, 1).sum()),
           "seed": seed, "nziehungen": nziehungen}
    for s in ("S1", "S2", "S3a", "S3b"):
        for schl, v in fasse(beob[s], null[s]).items():
            aus[f"{s}_{schl}"] = v
    j = pd.DataFrame({
        "zelle": list(X.index),
        "S1": beob["S1_je_zelle"],
        "S1_null_mittel": np.nanmean(null_s1_zelle, axis=0),
        "S1_null_sd": np.nanstd(null_s1_zelle, axis=0, ddof=1),
        "S3a": beob["S3a_je_zelle"],
        "S3a_null_mittel": np.nanmean(null_s3a_zelle, axis=0),
        "S3a_null_sd": np.nanstd(null_s3a_zelle, axis=0, ddof=1),
    })
    for s in ("S1", "S3a"):
        j[f"{s}_z"] = (j[s] - j[f"{s}_null_mittel"]) / j[f"{s}_null_sd"]
        j[f"{s}_mde80"] = j[f"{s}_null_mittel"] + 2.8 * j[f"{s}_null_sd"]
        j[f"{s}_ueber_mde80"] = j[s] > j[f"{s}_mde80"]
    aus["je_zelle"] = j
    return aus


def synthese_flip(matrizen: dict,
                  kennzahl: str = "S1",
                  nziehungen: int = 10000,
                  seed: int = SEED_D) -> dict:
    """Study synthesis against the joint donor-flip null (Fig. S3B).

    `matrizen` is a mapping
        study -> (M, ri, spender)
    with M = cells x set genes, ri = direction vector, spender = donor
    identifier per row. The null flips the sign of **entire donors** per
    draw (all cells of the same donor together, never separately) and
    recomputes the same statistic.

    Reported is the statistic averaged over studies against the joint null
    -- not an arithmetic mean over cohort z values.

    WARNING, one degenerate case that deserves naming: **S2 is invariant
    under donor flips.** Flipping row signs yields the new matrix D*M with D
    orthogonal-diagonal; the singular values and hence the PC1 share do not
    change. The flip null therefore has zero spread for S2 and is no null.
    In this case `entartet = True` is set and z/p are NaN -- **no** number
    is reported rather than an apparent one. The same holds for any study
    with only one donor: there S1 is undefined and the flip null has no
    effect.
    """
    rng = np.random.default_rng(seed)

    def kz(M, ri, sp):
        sp = np.asarray(sp)
        return _kennzahlen(M, ri, sp[:, None] != sp[None, :])[kennzahl]

    beob = {n: kz(M, ri, sp) for n, (M, ri, sp) in matrizen.items()}
    beob_mittel = float(np.nanmean(list(beob.values())))

    null = np.empty(nziehungen)
    je_studie = {n: np.empty(nziehungen) for n in matrizen}
    for i in range(nziehungen):
        werte = []
        for n, (M, ri, sp) in matrizen.items():
            sp_a = np.asarray(sp)
            uniq = pd.unique(sp_a)
            f = dict(zip(uniq, rng.choice([-1.0, 1.0], size=len(uniq))))
            fl = np.array([f[s] for s in sp_a])[:, None]
            v = kz(M * fl, ri, sp_a)
            je_studie[n][i] = v
            werte.append(v)
        null[i] = float(np.nanmean(werte))
    m, s = float(np.nanmean(null)), float(np.nanstd(null, ddof=1))
    entartet = not (s > 1e-12)
    p_ein = (1 + int((null >= beob_mittel).sum())) / (1 + len(null))
    zeilen = []
    for n in matrizen:
        nm = float(np.nanmean(je_studie[n]))
        ns = float(np.nanstd(je_studie[n], ddof=1))
        ent = not (ns > 1e-12)
        zeilen.append({"studie": n, "kennzahl": kennzahl,
                       "beobachtet": beob[n], "null_mittel": nm,
                       "null_sd": ns, "entartet": ent,
                       "z": np.nan if ent else (beob[n] - nm) / ns,
                       "mde80": np.nan if ent else nm + 2.8 * ns,
                       "ueber_mde80": (False if ent
                                       else bool(beob[n] > nm + 2.8 * ns))})
    return {"kennzahl": kennzahl, "n_studien": len(matrizen),
            "beobachtet": beob_mittel, "null_mittel": m, "null_sd": s,
            "entartet": entartet,
            "z": np.nan if entartet else (beob_mittel - m) / s,
            "p": np.nan if entartet else float(min(1.0, 2 * p_ein)),
            "mde80": np.nan if entartet else m + 2.8 * s,
            "seed": seed, "nziehungen": nziehungen,
            "je_studie": pd.DataFrame(zeilen)}

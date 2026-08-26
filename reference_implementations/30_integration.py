# -*- coding: utf-8 -*-
"""
30_integration.py -- the convergence: does internal convergence predict
movement on foreign measurement layers?

The step beyond the paper so far
--------------------------------
Until now the external test was a yes/no question posed to a fixed 173-gene
set. Here it becomes DOSE-RESPONSE. For each of the roughly 11 500
evaluable genes there is an internal continuous convergence strength:

    dWT_kons  = (2 v/n - 1) * sign   for the differentiation response
    iv_kons   = the same for the lesion response

The question: does the directed movement on an orthogonal layer increase
with internal convergence strength? And -- this is the actual test of the
paper title -- does this hold for `dWT` but NOT for `iv`?

With this, the sentence "Differentiation converges, lesions do not"
migrates from the RNA layer, where it was found, to measurement layers that
did not find it. That is the difference between an internal observation and
a finding.

Three analyses per layer:
  (A) rank correlation between `kons` and the directed layer difference.
  (B) decile curve: median of the layer difference per decile of `kons`.
  (C) the paired comparison `dWT` against `iv` on the same gene set and the
      same layer, with a joint null.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import ERGEBNISSE, SEED  # noqa: E402

AUS = ERGEBNISSE / "Z_integration"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []
NZIEH = 20000


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def lade_ebenen() -> dict[str, pd.Series]:
    """Directed layer difference per gene symbol.

    Sign convention: POSITIVE means "matches a gene that runs HIGH in
    differentiation". On the methylation layer the sign of the beta
    difference is flipped accordingly.
    """
    E: dict[str, pd.Series] = {}

    p = ERGEBNISSE / "B_atac" / "B2_GSE224251_genwerte_hart.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        E["ATAC GSE224251 osteogen (n=3)"] = t.osteogen

    p = ERGEBNISSE / "B_atac" / "B_atac_genwerte_T50.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        E["ATAC GSE332758 osteogen (n=1)"] = t.osteogen
        E["ATAC GSE332758 adipogen (n=1)"] = t.adipogen

    p = ERGEBNISSE / "B_atac" / "B3_GSE129031_genwerte_T10.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        E["H3K27ac GSE129031 chondrogen (n=2)"] = t.chondrogen

    p = ERGEBNISSE / "A_dnam" / "A_dnam_GSE33896_alle_gene_osteo.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        E["DNAm GSE33896 osteogen (n=3)"] = -t.delta      # sign flipped
    p = ERGEBNISSE / "A_dnam" / "A_dnam_GSE33896_alle_gene_myo.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        E["DNAm GSE33896 myogen (n=3)"] = -t.delta

    p = ERGEBNISSE / "A_dnam" / "A_dnam450_GSE129266_gene.csv"
    if p.exists():
        t = pd.read_csv(p, index_col=0)
        E["DNAm GSE129266 chondrogen (n=2)"] = -t.delta

    return {k: v[~v.index.duplicated()].dropna() for k, v in E.items()}


def dezilkurve(kons: pd.Series, delta: pd.Series, n_dez: int = 10) -> pd.DataFrame:
    d = pd.DataFrame({"kons": kons, "delta": delta}).dropna()
    d["dezil"] = pd.qcut(d.kons.rank(method="first"), n_dez, labels=False)
    g = d.groupby("dezil").agg(n=("delta", "size"),
                               kons_mitte=("kons", "median"),
                               delta_median=("delta", "median"),
                               delta_mittel=("delta", "mean"),
                               delta_se=("delta", lambda x: x.std(ddof=1) / np.sqrt(len(x))))
    return g.reset_index()


def gerichtete_konkordanz(kons: pd.Series, delta: pd.Series,
                          basis: pd.Series | None = None,
                          seed: int = SEED) -> dict:
    """Weighted directed agreement against a permutation null.

    Statistic: the mean signed rank of the layer difference, weighted by
    the convergence strength |kons|.

    `basis` switches the null from a free sign permutation to a STRATIFIED
    one: signs are swapped only within deciles of the baseline level. This
    rules out the statistic arising because highly convergent genes start
    systematically low and the orthogonal layer only sees the starting
    level -- the confounder Figure 2 of this paper warns about.
    """
    tab = {"kons": kons, "delta": delta}
    if basis is not None:
        tab["basis"] = basis
    d = pd.DataFrame(tab).dropna()
    d = d[(d.kons != 0) & (d.delta != 0)]
    if len(d) < 50:
        return {"n": len(d), "status": "too few"}
    w = d.kons.abs().values
    s = np.sign(d.kons.values)
    r = stats.rankdata(d.delta.values) / len(d) - 0.5
    beob = float(np.average(s * r, weights=w))
    rng = np.random.default_rng(seed)
    if basis is None:
        null = np.array([np.average(rng.permutation(s) * r, weights=w)
                         for _ in range(NZIEH)])
    else:
        dez = pd.qcut(d.basis.rank(method="first"), 10, labels=False).values
        bloecke = [np.flatnonzero(dez == b) for b in range(10)]
        null = np.empty(NZIEH)
        sp = s.copy()
        for i in range(NZIEH):
            for idx in bloecke:
                sp[idx] = rng.permutation(s[idx])
            null[i] = np.average(sp * r, weights=w)
    z = (beob - null.mean()) / null.std(ddof=1)
    p = min(1.0, 2 * (1 + (null >= beob).sum()) / (1 + len(null)))
    return {"n": len(d), "statistik": beob, "null": float(null.mean()),
            "null_sd": float(null.std(ddof=1)), "z": float(z), "p": float(p),
            "status": "ok"}


def main() -> None:
    log("=" * 78)
    log("Integration -- interne Konvergenz als Dosis, orthogonale Ebene als Wirkung")
    log("=" * 78)

    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    K = K[K.symbol.notna()].drop_duplicates("symbol").set_index("symbol")
    log("Internal gene map: %d genes with symbol" % len(K))

    E = lade_ebenen()
    log("Orthogonal layers: %d" % len(E))
    for k, v in E.items():
        log("  %-36s %6d genes | overlap with internal: %d"
            % (k, len(v), len(v.index.intersection(K.index))))

    zeilen, kurven = [], []
    log("\n" + "=" * 78)
    log("(A) and (C): rank correlation and the paired dWT-vs-iv comparison")
    log("=" * 78)
    log("%-36s %-9s %6s %9s %8s %8s %9s"
        % ("layer", "dose", "n", "statistic", "z", "p", "rho_S"))

    for ename, delta in E.items():
        gem = delta.index.intersection(K.index)
        if len(gem) < 200:
            log("%-36s too little overlap (%d)" % (ename, len(gem)))
            continue
        d = delta.reindex(gem)
        for dosis in ("dWT_kons", "iv_kons"):
            kons = K[dosis].reindex(gem)
            basis = K.basis_med.reindex(gem)
            r = gerichtete_konkordanz(kons, d)
            rb = gerichtete_konkordanz(kons, d, basis=basis)
            rho = stats.spearmanr(kons, d, nan_policy="omit")
            log("%-36s %-9s %6d %+9.4f %+8.2f %8.4g %+9.4f  | basisgesch. z %+6.2f p %.4g"
                % (ename, dosis.replace("_kons", ""), r.get("n", 0),
                   r.get("statistik", np.nan), r.get("z", np.nan),
                   r.get("p", np.nan), rho.statistic,
                   rb.get("z", np.nan), rb.get("p", np.nan)))
            zeilen.append(dict(ebene=ename, dosis=dosis, spearman_rho=rho.statistic,
                               spearman_p=rho.pvalue,
                               z_basisgeschichtet=rb.get("z"),
                               p_basisgeschichtet=rb.get("p"),
                               **{k: v for k, v in r.items() if k != "status"}))
            kur = dezilkurve(kons, d)
            kur["ebene"] = ename
            kur["dosis"] = dosis
            kurven.append(kur)
        log("")

    Z = pd.DataFrame(zeilen)
    Z.to_csv(AUS / "Z_integration_statistik.csv", index=False)
    pd.concat(kurven, ignore_index=True).to_csv(AUS / "Z_integration_dezile.csv",
                                                index=False)

    # ------------------------------------------------------ the core claim
    log("=" * 78)
    log("THE CORE COMPARISON: same layer, same genes, two doses")
    log("=" * 78)
    log("Both columns baseline-stratified -- the free null is in the CSV.")
    log("%-36s %10s %10s %12s" % ("layer", "z(dWT)", "z(iv)", "difference"))
    for ename in Z.ebene.unique():
        a = Z[(Z.ebene == ename) & (Z.dosis == "dWT_kons")]
        b = Z[(Z.ebene == ename) & (Z.dosis == "iv_kons")]
        if len(a) and len(b):
            za, zb = a.z_basisgeschichtet.iloc[0], b.z_basisgeschichtet.iloc[0]
            log("%-36s %+10.2f %+10.2f %+12.2f" % (ename, za, zb, za - zb))
    log("")
    log("Reading: where z(dWT) clearly exceeds z(iv), the sentence")
    log('"Differentiation converges, lesions do not" is repeated on a')
    log("measurement layer that did not produce it.")

    # ------------------------------------- confounder baseline level
    log("\n--- Confounder control: the baseline level -----------------------")
    log("The per-gene z scale forces cor(basis, dWT) < 0. If the orthogonal")
    log("layer only saw the baseline level, the same finding would arise")
    log("without any convergence. The partial correlation against `basis_med`")
    log("is therefore additionally computed.")
    log("%-36s %-9s %10s %10s" % ("layer", "dose", "rho raw", "rho partial"))
    part = []
    for ename, delta in E.items():
        gem = delta.index.intersection(K.index)
        if len(gem) < 200:
            continue
        d = delta.reindex(gem)
        b = K.basis_med.reindex(gem)
        for dosis in ("dWT_kons", "iv_kons"):
            kons = K[dosis].reindex(gem)
            ok = kons.notna() & d.notna() & b.notna()
            rk = stats.rankdata
            x, y, z0 = rk(kons[ok]), rk(d[ok]), rk(b[ok])
            rxy = np.corrcoef(x, y)[0, 1]
            rxz = np.corrcoef(x, z0)[0, 1]
            ryz = np.corrcoef(y, z0)[0, 1]
            rp = (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
            log("%-36s %-9s %+10.4f %+10.4f"
                % (ename, dosis.replace("_kons", ""), rxy, rp))
            part.append(dict(ebene=ename, dosis=dosis, rho_roh=rxy,
                             rho_partiell=rp, rho_kons_basis=rxz,
                             rho_delta_basis=ryz, n=int(ok.sum())))
    pd.DataFrame(part).to_csv(AUS / "Z_integration_partiell.csv", index=False)

    (AUS / "Z_integration_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

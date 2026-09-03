# -*- coding: utf-8 -*-
"""
80_lineage_independence_crossarm.py -- WS6, checks P1 and P2.

Question: does the fixed 173-gene module run independently of the (often
failed) lineage-marker calibration of the 18 perturbation datasets?

Data source (NOT recomputed, verbatim from the frozen tables):
  * per dataset `dWT` per gene  -> 20_Exploration/derived_data/20d_gene_XX.csv
    (the same source that 10_calibration_18_datasets.py uses for lineage calibration)
  * the lineage calibration per dataset -> derived_data/M_kalibrierung/eichung_achtzehn.csv
  * the module (173 genes, direction ri) -> 00_shared/_module.py (MODUL)
  * the statistic -> 00_shared/_module.py `konkordanz()` (sign test against
    background-drawn null), the same implementation as everywhere in the project.

P1: does the module run in the same direction in both represented arms
    (osteogenic, chondrogenic)? (ADIPOGEN/MYOGEN are NOT represented among
    the 18 -- this is noted explicitly, not concealed.)
P2: does the module also carry in exactly those datasets whose OWN
    lineage-marker calibration failed?
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pandas as pd

HIER = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parents[1] / "00_shared"))
from _module import MODUL, konkordanz  # noqa: E402

WURZEL = HIER.parents[1]
ANTRAG = WURZEL.parent
AUS = WURZEL / "derived_data" / "followup"
AUS.mkdir(parents=True, exist_ok=True)

# Not part of this repository: the author's session tree, holding the
# per-run intermediates of the original exploration. Point at it with
# SCHERENPAPER_SITZUNGEN. This script is marked "needs raw data".
SITZUNGEN = (pathlib.Path(os.environ["SCHERENPAPER_SITZUNGEN"])
             if os.environ.get("SCHERENPAPER_SITZUNGEN") else ANTRAG)
GENE20D = (SITZUNGEN /
           "20_Exploration" / "derived_data")
EICHUNG = pd.read_csv(WURZEL / "derived_data" / "M_kalibrierung" / "eichung_achtzehn.csv")
KO = pd.read_csv(WURZEL / "derived_data" / "manuscript" / "f1_kohorte.csv")
arm_map = dict(zip(KO.punkt, KO.arm))
name_map = dict(zip(KO.punkt, KO.datensatz))

RI = MODUL.set_index("symbol")["ri"]


def gencode_karte() -> dict:
    """The same Ensembl->symbol map as in 10_calibration_18_datasets.py."""
    import gzip
    import re
    DATEN = WURZEL / "data_raw"
    for p in sorted((DATEN / "_referenz").glob("*.gtf*")):
        m = {}
        op = gzip.open if p.suffix == ".gz" else open
        with op(p, "rt", encoding="utf-8", errors="replace") as f:
            for ln in f:
                if ln[0] == "#" or "\tgene\t" not in ln:
                    continue
                g = re.search(r'gene_id "([^".]+)', ln)
                s = re.search(r'gene_name "([^"]+)', ln)
                if g and s:
                    m[g.group(1)] = s.group(1)
        if m:
            return m
    raise RuntimeError("no Gencode reference found")


def main() -> None:
    karte = gencode_karte()
    zeilen = []
    basis_je_punkt = {}
    dwt_je_punkt = {}

    for f in sorted(GENE20D.glob("20d_gene_*.csv")):
        G = pd.read_csv(f)
        p = int(G.punkt.iloc[0])
        G["symbol"] = [karte.get(str(g).split(".")[0]) for g in G.gen]
        G = G[G.symbol.notna() & G.dWT.notna()]
        dwt = G.groupby("symbol").dWT.median()
        basis = G.groupby("symbol").basis.median()
        dwt_je_punkt[p] = dwt
        basis_je_punkt[p] = basis

        r = konkordanz(dwt.reindex(RI.index).dropna(), RI,
                        hintergrund=dwt, nziehungen=3000)
        r.update(punkt=p, datensatz=name_map.get(p, ""), arm=arm_map.get(p, ""))
        zeilen.append(r)
        pd.DataFrame(zeilen).to_csv(AUS / "ws6_p1p2_modul_je_datensatz_zwischenstand.csv", index=False)
        print("  point %2d  %-28s %-11s  concordance %.3f  z %+6.2f  p %8.4g  (n=%s)"
              % (p, str(name_map.get(p, ""))[:28], arm_map.get(p, ""),
                 r.get("konkordanz", np.nan), r.get("konkordanz_z", np.nan),
                 r.get("konkordanz_p", np.nan), r.get("n", "-")))

    T = pd.DataFrame(zeilen)
    T = T.merge(EICHUNG[["punkt", "bestanden", "n_gene_messbar"]],
                on="punkt", how="left", suffixes=("", "_eich"))
    T.to_csv(AUS / "ws6_p1p2_modul_je_datensatz.csv", index=False)

    print("\n" + "=" * 78)
    print("P2 -- does the module carry in those datasets whose OWN")
    print("      lineage-marker calibration did NOT pass?")
    print("=" * 78)
    ok = T[T.status == "ok"].copy()
    durchgefallen = ok[~ok.bestanden]
    bestanden = ok[ok.bestanden]
    print("  passed (n=%d): concordance-z median %.2f, share p<0.05 directed: %d/%d"
          % (len(bestanden), bestanden.konkordanz_z.median(),
             int((bestanden.konkordanz_p < 0.05).sum()), len(bestanden)))
    print("  failed (n=%d): concordance-z median %.2f, share p<0.05 directed: %d/%d"
          % (len(durchgefallen), durchgefallen.konkordanz_z.median(),
             int((durchgefallen.konkordanz_p < 0.05).sum()), len(durchgefallen)))

    # pooled test ONLY over the failed datasets: a joint module-dWT over
    # these datasets (median per gene over the failed points), against the
    # same baseline-stratified null.
    punkte_fail = durchgefallen.punkt.tolist()
    punkte_pass = bestanden.punkt.tolist()

    def gepoolt(punkte, label):
        M = pd.concat([dwt_je_punkt[p].rename(p) for p in punkte if p in dwt_je_punkt], axis=1)
        pooled = M.median(axis=1, skipna=True)
        B = pd.concat([basis_je_punkt[p].rename(p) for p in punkte if p in basis_je_punkt], axis=1)
        basis_pooled = B.median(axis=1, skipna=True)
        r = konkordanz(pooled.reindex(RI.index).dropna(), RI,
                        hintergrund=pooled, nziehungen=8000)
        r["label"] = label
        r["n_datensaetze"] = len(punkte)
        print("  POOLED %-14s (n=%d datasets): concordance %.3f  z %+6.2f  p %8.4g  n_gene=%s"
              % (label, len(punkte), r.get("konkordanz", np.nan),
                 r.get("konkordanz_z", np.nan), r.get("konkordanz_p", np.nan),
                 r.get("n", "-")))
        return r

    r_fail = gepoolt(punkte_fail, "durchgefallen")
    pd.DataFrame([r_fail]).to_csv(
        AUS / "ws6_p2_gepoolt_bestanden_vs_durchgefallen.csv", index=False)
    r_pass = gepoolt(punkte_pass, "bestanden")
    pd.DataFrame([r_fail, r_pass]).to_csv(
        AUS / "ws6_p2_gepoolt_bestanden_vs_durchgefallen.csv", index=False)

    print("\n" + "=" * 78)
    print("P1 -- cross-concordance between the arms (osteogenic, chondrogenic);")
    print("      ADIPOGEN/MYOGEN are NOT represented among the 18 datasets.")
    print("=" * 78)
    for a in ["osteogen", "chondrogen"]:
        punkte_a = [p for p in dwt_je_punkt if arm_map.get(p) == a]
        gepoolt(punkte_a, "arm_" + a)

    # cross-correlation of the two arm-pooled dWT vectors on the module
    # genes, against a null model of equal-sized random gene sets from the
    # joint background (Spearman correlation of the two arm vectors on
    # module genes vs on random genes).
    osteo_p = [p for p in dwt_je_punkt if arm_map.get(p) == "osteogen"]
    chond_p = [p for p in dwt_je_punkt if arm_map.get(p) == "chondrogen"]
    O = pd.concat([dwt_je_punkt[p].rename(p) for p in osteo_p], axis=1).median(axis=1, skipna=True)
    C = pd.concat([dwt_je_punkt[p].rename(p) for p in chond_p], axis=1).median(axis=1, skipna=True)
    gem = O.index.intersection(C.index)
    O, C = O.reindex(gem), C.reindex(gem)
    modulgene = [g for g in RI.index if g in gem]
    from scipy import stats
    rho_module, _ = stats.spearmanr(O[modulgene], C[modulgene])
    rng = np.random.default_rng(20260822)
    k = len(modulgene)
    null_rho = np.empty(5000)
    idx_all = np.arange(len(gem))
    Ov, Cv = O.values, C.values
    for i in range(5000):
        idx = rng.choice(idx_all, size=k, replace=False)
        null_rho[i] = stats.spearmanr(Ov[idx], Cv[idx])[0]
    z = (rho_module - null_rho.mean()) / null_rho.std(ddof=1)
    p_ein = (1 + int((null_rho >= rho_module).sum())) / (1 + len(null_rho))
    p = min(1.0, 2 * p_ein)
    print("  arm-arm Spearman rho of the module genes: %.3f (null mean %.3f +- %.3f, z %+.2f, p %.4g)"
          % (rho_module, null_rho.mean(), null_rho.std(ddof=1), z, p))
    pd.DataFrame([{"rho_modulegene": rho_module, "null_mittel": null_rho.mean(),
                   "null_sd": null_rho.std(ddof=1), "z": z, "p": p,
                   "n_modulegene_gemeinsam": k, "n_hintergrund": len(gem),
                   "n_datensaetze_osteo": len(osteo_p),
                   "n_datensaetze_chondro": len(chond_p)}]).to_csv(
        AUS / "ws6_p1_arm_kreuzkonkordanz.csv", index=False)


if __name__ == "__main__":
    main()

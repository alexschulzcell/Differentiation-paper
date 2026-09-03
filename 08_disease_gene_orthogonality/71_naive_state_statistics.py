# -*- coding: utf-8 -*-
"""
71_naive_state_statistics.py -- P1/P2/P2b/P3/P4 based on the already-built
derived_data/followup/ws3_zellen_tabelle.csv (from 70_naive_state_predicts_response.py) PLUS
the full gene vectors, which are built here once more from the same
54b_zellen loaders but immediately cast into fixed numpy matrices
(module genes x cells) -- this makes the permutation null vectorized and
fast (no pandas reindex per draw).
"""
from __future__ import annotations

import importlib
import pathlib
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

WURZEL = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WURZEL / "00_shared"))
AUS = WURZEL / "derived_data" / "followup"
AUS.mkdir(parents=True, exist_ok=True)

t0 = time.time()
Z54 = importlib.import_module("54b_zellen")
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import MODUL  # noqa: E402
print(f"[t={time.time()-t0:.1f}s] modules imported", flush=True)

MARKERSAETZE = {"OSTEOGEN": OSTEOGEN, "ADIPOGEN": ADIPOGEN,
                "MYOGEN": MYOGEN, "CHONDROGEN": CHONDROGEN, "NAIV": NAIV}
MODULGENE = list(MODUL.symbol)
RI = MODUL.set_index("symbol")["ri"].reindex(MODULGENE).to_numpy(dtype=float)
SEED = 20260822
NZIEH = 8000  # reduced from 20000 -- defensible at n=7/14 donors


def _identitaet(X):
    return X


MARKERGENE = sorted(set(OSTEOGEN) | set(CHONDROGEN))


def baue_matrizen():
    """Per cell: naiv_z (vector over MODULGENE), dwt (same),
    naiv_roh (raw-abundance mean of the naive samples, unweighted), plus
    naiv_marker (vector over OSTEOGEN+CHONDROGEN markers, z scale -- these
    genes generally do NOT lie in the 173-gene module and need their own
    column)."""
    zmat_original = Z54.zmat
    zellen_meta = []
    naiv_cols, dwt_cols, roh_cols, marker_cols = [], [], [], []
    for gse, (lader, kurz, herkunft) in Z54.STUDIEN.items():
        Z, zellen, proben = lader()
        Z54.zmat = _identitaet
        try:
            Zroh, _, _ = lader()
        finally:
            Z54.zmat = zmat_original
        for z in zellen:
            fehlt = [k for k in ("k_naiv", "k_diff", "l_naiv", "l_diff") if not z[k]]
            if fehlt:
                continue
            d, i, b = Z54.zelle(Z, z["k_naiv"], z["k_diff"], z["l_naiv"], z["l_diff"])
            nmod = len(set(MODUL.symbol) & set(d.dropna().index))
            if nmod < 60:
                continue
            name = f"{z['spender']}_{z['achse'][:6]}"
            b_roh = (Zroh[list(z["k_naiv"])].mean(axis=1)
                     if all(c in Zroh.columns for c in z["k_naiv"])
                     else pd.Series(dtype=float))
            naiv_cols.append(b.reindex(MODULGENE).to_numpy(dtype=float))
            dwt_cols.append(d.reindex(MODULGENE).to_numpy(dtype=float))
            roh_cols.append(b_roh.reindex(MODULGENE).to_numpy(dtype=float))
            marker_cols.append(b.reindex(MARKERGENE).to_numpy(dtype=float))
            zellen_meta.append(dict(studie=gse, zelle=name, spender=z["spender"],
                                    achse=z["achse"], e2=z["e2"],
                                    laesion=z["laesion"], herkunft=herkunft))
        print(f"[t={time.time()-t0:.1f}s] {gse} loaded, {len(zellen_meta)} cells cumulative", flush=True)
    naiv_mat = np.column_stack(naiv_cols)      # genes x cells, z scale
    dwt_mat = np.column_stack(dwt_cols)
    roh_mat = np.column_stack(roh_cols)
    marker_mat = np.column_stack(marker_cols)  # marker genes x cells, z scale
    meta = pd.DataFrame(zellen_meta)
    return naiv_mat, dwt_mat, roh_mat, marker_mat, meta


def modulwert_vec(mat, ri=RI):
    """ri-weighted column mean, NaN-robust."""
    m = np.where(np.isnan(mat), np.nan, mat * ri[:, None])
    return np.nanmean(m, axis=0)


def amplitude_vec(mat):
    return np.nanmean(np.abs(mat), axis=0)


def marker_score(mat_z, genliste, alle_gene):
    idx = [i for i, g in enumerate(alle_gene) if g in genliste]
    if len(idx) < 3:
        return np.full(mat_z.shape[1], np.nan)
    return np.nanmean(mat_z[idx, :], axis=0)


def spearman_col(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return np.nan
    return float(stats.spearmanr(x[ok], y[ok]).correlation)


def p1_p2(naiv_mat, dwt_mat, meta, idx, label, seed=SEED, nziehungen=NZIEH):
    idx = np.asarray(idx)
    n = len(idx)
    if n < 4:
        return {"variante": label, "n": n, "status": "too few donors"}
    naiv_score = modulwert_vec(naiv_mat[:, idx])
    dwt_amp = amplitude_vec(dwt_mat[:, idx])
    dwt_dir = modulwert_vec(dwt_mat[:, idx])
    rho_amp = spearman_col(naiv_score, dwt_amp)
    rho_dir = spearman_col(naiv_score, dwt_dir)

    diff_mat = dwt_mat[:, idx] + naiv_mat[:, idx]     # diff = dWT + naive
    naiv_sub = naiv_mat[:, idx]
    studien = meta["studie"].to_numpy()[idx]
    gruppen = {s: np.flatnonzero(studien == s) for s in np.unique(studien)}

    rng = np.random.default_rng(seed)
    null_amp = np.empty(nziehungen)
    null_dir = np.empty(nziehungen)
    perm_idx = np.tile(np.arange(n), (nziehungen, 1))
    for s, pos in gruppen.items():
        if len(pos) < 2:
            continue
        for it in range(nziehungen):
            perm_idx[it, pos] = rng.permutation(pos)

    for it in range(nziehungen):
        diff_perm = diff_mat[:, perm_idx[it]]
        dwt_null = diff_perm - naiv_sub
        null_amp[it] = spearman_col(naiv_score, amplitude_vec(dwt_null))
        null_dir[it] = spearman_col(naiv_score, modulwert_vec(dwt_null))

    def fasse(beob, null):
        null = null[np.isfinite(null)]
        m, s = float(np.mean(null)), float(np.std(null, ddof=1))
        z = (beob - m) / s if s > 0 else np.nan
        p = (1 + int((np.abs(null - m) >= abs(beob - m)).sum())) / (1 + len(null))
        return dict(beobachtet=beob, null_mittel=m, null_sd=s, z=float(z),
                    p=float(min(1.0, p)))

    aus = {"variante": label, "n": n, "status": "ok"}
    for k, v in fasse(rho_amp, null_amp).items():
        aus[f"amplitude_{k}"] = v
    for k, v in fasse(rho_dir, null_dir).items():
        aus[f"richtung_{k}"] = v
    return aus


def p2b_roh(roh_mat, dwt_mat, meta, idx, label, seed=SEED, nziehungen=NZIEH):
    idx = np.asarray(idx)
    n = len(idx)
    if n < 4:
        return {"variante": label, "n": n, "status": "too few donors"}
    naiv_roh = np.nanmean(roh_mat[:, idx], axis=0)   # unweighted mean
    dwt_amp = amplitude_vec(dwt_mat[:, idx])
    rho = spearman_col(naiv_roh, dwt_amp)
    studien = meta["studie"].to_numpy()[idx]
    gruppen = {s: np.flatnonzero(studien == s) for s in np.unique(studien)}
    rng = np.random.default_rng(seed)
    null = np.empty(nziehungen)
    for it in range(nziehungen):
        dwt_perm = dwt_amp.copy()
        for s, pos in gruppen.items():
            if len(pos) < 2:
                continue
            dwt_perm[pos] = dwt_amp[rng.permutation(pos)]
        null[it] = spearman_col(naiv_roh, dwt_perm)
    nn = null[np.isfinite(null)]
    m, s = float(np.mean(nn)), float(np.std(nn, ddof=1))
    z = (rho - m) / s if s > 0 else np.nan
    p = (1 + int((np.abs(nn - m) >= abs(rho - m)).sum())) / (1 + len(nn))
    return {"variante": label, "n": n, "status": "ok", "beobachtet": rho,
            "null_mittel": m, "null_sd": s, "z": float(z), "p": float(min(1.0, p))}


def main():
    naiv_mat, dwt_mat, roh_mat, marker_mat, meta = baue_matrizen()
    print(f"[t={time.time()-t0:.1f}s] matrices built: {naiv_mat.shape}", flush=True)

    eich = pd.read_csv(WURZEL / "derived_data" / "M_donoren" / "eichung.csv").set_index("zelle")["bestanden"]
    meta["bestanden"] = meta["zelle"].map(eich).fillna(False).astype(bool)
    meta.to_csv(AUS / "ws3_zellen_meta.csv", index=False)

    idx_geeicht = np.flatnonzero(meta["bestanden"].to_numpy())
    idx_alle = np.arange(len(meta))
    print(f"[t={time.time()-t0:.1f}s] calibrated={len(idx_geeicht)} all={len(idx_alle)}", flush=True)

    # ---------------------------------------------------------------- P1/P2
    erg = [p1_p2(naiv_mat, dwt_mat, meta, idx_geeicht, "primaer_n7_geeicht"),
           p1_p2(naiv_mat, dwt_mat, meta, idx_alle, "sensitivitaet_n14_alle")]
    pd.DataFrame(erg).to_csv(AUS / "ws3_p1_p2_korrelation_z.csv", index=False)
    print(f"[t={time.time()-t0:.1f}s] P1/P2 done", flush=True)
    for r in erg:
        print("  ", r)

    # ---------------------------------------------------------------- P2b
    erg2b = [p2b_roh(roh_mat, dwt_mat, meta, idx_geeicht, "primaer_n7_geeicht"),
             p2b_roh(roh_mat, dwt_mat, meta, idx_alle, "sensitivitaet_n14_alle")]
    pd.DataFrame(erg2b).to_csv(AUS / "ws3_p2b_korrelation_roh.csv", index=False)
    print(f"[t={time.time()-t0:.1f}s] P2b done", flush=True)
    for r in erg2b:
        print("  ", r)

    # ---------------------------------------------------------------- P3
    achse = meta["achse"].to_numpy()
    sc_osteogen = marker_score(marker_mat, OSTEOGEN, MARKERGENE)
    sc_chondrogen = marker_score(marker_mat, CHONDROGEN, MARKERGENE)
    eigen = np.where(achse == "osteogen", sc_osteogen, sc_chondrogen)
    fremd = np.where(achse == "osteogen", sc_chondrogen, sc_osteogen)
    achsendiff = eigen - fremd
    p3_rows = []
    for label, idx in (("primaer_n7_geeicht", idx_geeicht), ("sensitivitaet_n14_alle", idx_alle)):
        ach = achse[idx]
        ad = achsendiff[idx]
        ok = np.isfinite(ad)
        ach, ad = ach[ok], ad[ok]
        if len(set(ach)) < 2 or len(ach) < 4:
            p3_rows.append({"variante": label, "status": "one axis or too few cells", "n": len(ach)})
            continue
        y = (ach == "osteogen").astype(int)
        u = stats.mannwhitneyu(ad[y == 1], ad[y == 0], alternative="two-sided")
        auc = float(u.statistic / (max(1, (y == 1).sum()) * max(1, (y == 0).sum())))
        p3_rows.append({"variante": label, "status": "ok", "n": len(ach),
                        "n_osteogen": int((y == 1).sum()), "n_chondrogen": int((y == 0).sum()),
                        "auc": auc, "mannwhitney_p": float(u.pvalue)})
    pd.DataFrame(p3_rows).to_csv(AUS / "ws3_p3_achsenvorhersage.csv", index=False)
    print(f"[t={time.time()-t0:.1f}s] P3 done:", p3_rows, flush=True)

    # ---------------------------------------------------------------- P4
    naiv_score_all = modulwert_vec(naiv_mat)
    naiv_roh_all = np.nanmean(roh_mat, axis=0)
    p4_rows = []
    for label, vals in (("naiv_modulewert_z", naiv_score_all),
                        ("naiv_modulewert_roh", naiv_roh_all),
                        ("naiv_marker_achsendiff", achsendiff)):
        a = vals[meta["bestanden"].to_numpy()]
        b = vals[~meta["bestanden"].to_numpy()]
        a, b = a[np.isfinite(a)], b[np.isfinite(b)]
        if len(a) < 3 or len(b) < 3:
            p4_rows.append({"variable": label, "status": "too few cells",
                            "n_bestanden": len(a), "n_durchgefallen": len(b)})
            continue
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        p4_rows.append({"variable": label, "status": "ok", "n_bestanden": len(a),
                        "n_durchgefallen": len(b), "median_bestanden": float(np.median(a)),
                        "median_durchgefallen": float(np.median(b)),
                        "U": float(u.statistic), "p": float(u.pvalue)})
    pd.DataFrame(p4_rows).to_csv(AUS / "ws3_p4_eichung_naivunterschied.csv", index=False)
    print(f"[t={time.time()-t0:.1f}s] P4 done:", p4_rows, flush=True)
    print(f"[t={time.time()-t0:.1f}s] ALL DONE", flush=True)


if __name__ == "__main__":
    main()

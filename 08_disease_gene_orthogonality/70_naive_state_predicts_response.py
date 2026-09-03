# -*- coding: utf-8 -*-
"""
70_naive_state_predicts_response.py -- work package WS3: does the naive state (day 0) carry
information about the later differentiation response of the same donor?

Uses the six loaders from 08_disease_gene_orthogonality/31_donor_cells_build_calibrate.py (not rewritten,
imported) to form, per donor cell, the module value IN THE NAIVE ARM
(`basis`, z scale) -- 31_donor_cells_build_calibrate.py stores this value only AVERAGED OVER
CELLS per study in zellen.pkl; here the same loop is repeated and the value
is recorded PER CELL (it exists in the code as `basis_sp[name]` but is
averaged away before saving -- 31_donor_cells_build_calibrate.py line ~372).

P1  correlation module value(naive) vs amplitude/direction of dWT, per donor.
P2  circularity trap: dWT is ALGEBRAICALLY `diff - naive`. Any correlation
    of naive with dWT thereby has a built-in negative bias, without any
    biology. Null model: `diff := dWT_obs + basis_obs` is PERMUTED between
    donors per study (the naive value stays with its real donor),
    dWT_null := diff_permuted - naive_real is recomputed, correlation recomputed.
    This reproduces exactly the subtraction mechanics but destroys the donor
    pairing -- the only quantity that could carry real prediction.
P2b alternative scale: naive on ABSOLUTE scale (log2 abundance before z
    standardization per study), response still on the contrast scale.
P3  composition: share of lineage markers of the own axis in the naive state
    as predictor instead of the global module value.
P4  do the 7 cells that fail their own calibration already differ in the
    naive state from the 7 passing ones?

Only the 7 calibrated cells carry numbers (rule of the assignment); the
other 7 are carried along and reported, but P1/P2 are computed primarily on
n=7, n=14 supplied as sensitivity.
"""
from __future__ import annotations

import importlib
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

WURZEL = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WURZEL / "00_shared"))
AUS = WURZEL / "derived_data" / "followup"
AUS.mkdir(parents=True, exist_ok=True)

Z54 = importlib.import_module("54b_zellen")
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import MODUL  # noqa: E402

SEED = 20260822
NZIEH = 20000
MARKERSAETZE = {"OSTEOGEN": OSTEOGEN, "ADIPOGEN": ADIPOGEN,
                "MYOGEN": MYOGEN, "CHONDROGEN": CHONDROGEN, "NAIV": NAIV}

RI = dict(zip(MODUL.symbol, MODUL.ri))
MODULGENE = list(RI.keys())

# ---------------------------------------------------------------------------
# 1. rebuild all 14 cells -- once on the z scale (original), once on absolute
#    log2 abundance scale (zmat replaced by identity)
# ---------------------------------------------------------------------------
def _identitaet(X):
    return X


def baue_zellen():
    """Reproduces 54b_zellen.main()'s core loop, but keeps basis_sp
    PER CELL (not averaged) and additionally the raw-abundance variant."""
    zeilen = []
    zmat_original = Z54.zmat
    for gse, (lader, kurz, herkunft) in Z54.STUDIEN.items():
        # -- z scale (as in the original) --
        Z, zellen, proben = lader()
        # -- raw abundance scale: replace zmat by identity, reload --
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
            # raw-abundance baseline: simple mean of the naive control samples
            # on the log2 scale (no z per gene) -- same samples, same matrix,
            # just without standardization.
            b_roh = Zroh[list(z["k_naiv"])].mean(axis=1) if all(
                c in Zroh.columns for c in z["k_naiv"]) else pd.Series(dtype=float)
            zeilen.append(dict(studie=gse, zelle=name, spender=z["spender"],
                               achse=z["achse"], e2=z["e2"], laesion=z["laesion"],
                               herkunft=herkunft, dwt=d, iv=i, basis_z=b,
                               basis_roh=b_roh))
    return zeilen


def modulwert(serie: pd.Series, gene=None, gewichtet=True) -> float:
    """Mean (ri-weighted) value of a gene series over the module."""
    gene = gene or MODULGENE
    s = serie.reindex(gene).dropna()
    if len(s) < 8:
        return np.nan
    if gewichtet:
        ri = pd.Series(RI).reindex(s.index)
        return float((s * ri).mean())
    return float(s.mean())


def amplitude(serie: pd.Series, gene=None) -> float:
    gene = gene or MODULGENE
    s = serie.reindex(gene).dropna()
    if len(s) < 8:
        return np.nan
    return float(s.abs().mean())


def marker_anteil(serie: pd.Series, satz: list) -> float:
    """Mean of a marker set on a gene series (z scale)."""
    s = serie.reindex(satz).dropna()
    if len(s) < 3:
        return np.nan
    return float(s.mean())


# ---------------------------------------------------------------------------
def haupt():
    zellen = baue_zellen()
    eich = pd.read_csv(WURZEL / "derived_data" / "M_donoren" / "eichung.csv")
    eich = eich.set_index("zelle")["bestanden"]

    zeilen = []
    for z in zellen:
        row = dict(studie=z["studie"], zelle=z["zelle"], spender=z["spender"],
                   achse=z["achse"], e2=z["e2"], laesion=z["laesion"],
                   herkunft=z["herkunft"],
                   bestanden=bool(eich.get(z["zelle"], False)))
        row["naiv_modulewert_z"] = modulwert(z["basis_z"])
        row["naiv_modulewert_roh"] = modulwert(z["basis_roh"], gewichtet=False)
        row["dwt_modulewert"] = modulwert(z["dwt"])
        row["dwt_amplitude"] = amplitude(z["dwt"])
        row["n_modulegene"] = int(z["dwt"].reindex(MODULGENE).notna().sum())
        satz_eigen = "OSTEOGEN" if z["achse"] == "osteogen" else "CHONDROGEN"
        satz_fremd = "CHONDROGEN" if z["achse"] == "osteogen" else "OSTEOGEN"
        row["naiv_marker_eigen"] = marker_anteil(z["basis_z"], MARKERSAETZE[satz_eigen])
        row["naiv_marker_fremd"] = marker_anteil(z["basis_z"], MARKERSAETZE[satz_fremd])
        row["naiv_marker_achsendiff"] = (row["naiv_marker_eigen"] - row["naiv_marker_fremd"]
                                          if pd.notna(row["naiv_marker_eigen"]) and
                                          pd.notna(row["naiv_marker_fremd"]) else np.nan)
        zeilen.append(row)
    T = pd.DataFrame(zeilen)
    T.to_csv(AUS / "ws3_zellen_tabelle.csv", index=False)

    # ---------------------------------------------------------------- P1/P2
    rng = np.random.default_rng(SEED)
    daten_je_studie = {z["studie"]: z for z in zellen}
    studien_index = {}
    for i, z in enumerate(zellen):
        studien_index.setdefault(z["studie"], []).append(i)

    def diff_serie(z):
        # diff := dWT + naive (algebraic inverse of dWT = diff - naive)
        return z["dwt"].add(z["basis_z"], fill_value=np.nan)

    diffe = [diff_serie(z) for z in zellen]

    def korrelationslauf(teilmenge_idx, nziehungen=NZIEH, seed=SEED):
        idx = list(teilmenge_idx)
        if len(idx) < 4:
            return {"n": len(idx), "status": "too few donors"}
        naiv = np.array([modulwert(zellen[i]["basis_z"]) for i in idx])
        dwt_amp = np.array([amplitude(zellen[i]["dwt"]) for i in idx])
        dwt_dir = np.array([modulwert(zellen[i]["dwt"]) for i in idx])
        ok = np.isfinite(naiv) & np.isfinite(dwt_amp)
        naiv, dwt_amp, dwt_dir = naiv[ok], dwt_amp[ok], dwt_dir[ok]
        idx_ok = [i for i, o in zip(idx, ok) if o]
        n = len(idx_ok)
        if n < 4:
            return {"n": n, "status": "too few donors after NA filter"}
        rho_amp = float(stats.spearmanr(naiv, dwt_amp).correlation)
        rho_dir = float(stats.spearmanr(naiv, dwt_dir).correlation)

        rng_l = np.random.default_rng(seed)
        null_amp = np.empty(nziehungen)
        null_dir = np.empty(nziehungen)
        # per study, permute among its own donors (diff stays with the
        # study, naive stays with the real donor)
        by_study = {}
        for pos, i in enumerate(idx_ok):
            by_study.setdefault(zellen[i]["studie"], []).append(pos)

        naiv_arr = naiv  # fixed order like idx_ok
        for it in range(nziehungen):
            dwt_amp_null = np.empty(n)
            dwt_dir_null = np.empty(n)
            for studie, posliste in by_study.items():
                posliste_arr = np.array(posliste)
                perm = rng_l.permutation(posliste_arr)
                for p_orig, p_perm in zip(posliste_arr, perm):
                    i_orig = idx_ok[p_orig]
                    i_perm = idx_ok[p_perm]
                    diff_perm = diffe[i_perm]
                    naiv_echt = zellen[i_orig]["basis_z"]
                    dwt_null_serie = diff_perm.reindex(naiv_echt.index) - naiv_echt
                    dwt_amp_null[p_orig] = amplitude(dwt_null_serie)
                    dwt_dir_null[p_orig] = modulwert(dwt_null_serie)
            null_amp[it] = stats.spearmanr(naiv_arr, dwt_amp_null).correlation
            null_dir[it] = stats.spearmanr(naiv_arr, dwt_dir_null).correlation

        def fasse(beob, null):
            null = null[np.isfinite(null)]
            m, s = float(np.mean(null)), float(np.std(null, ddof=1))
            z = (beob - m) / s if s > 0 else np.nan
            p = (1 + int((np.abs(null - m) >= abs(beob - m)).sum())) / (1 + len(null))
            return {"beobachtet": beob, "null_mittel": m, "null_sd": s,
                    "z": float(z), "p": float(min(1.0, p)),
                    "mde80_pos": m + 2.8 * s, "mde80_neg": m - 2.8 * s}

        aus = {"n": n, "status": "ok"}
        for schl, v in fasse(rho_amp, null_amp).items():
            aus[f"amplitude_{schl}"] = v
        for schl, v in fasse(rho_dir, null_dir).items():
            aus[f"richtung_{schl}"] = v
        return aus

    idx_alle = list(range(len(zellen)))
    idx_geeicht = [i for i, z in enumerate(zellen)
                   if bool(eich.get(z["zelle"], False))]

    erg_p1 = []
    for label, idx in (("primaer_n7_geeicht", idx_geeicht),
                        ("sensitivitaet_n14_alle", idx_alle)):
        r = korrelationslauf(idx)
        r["variante"] = label
        r["skala_naiv"] = "z (studienintern)"
        erg_p1.append(r)
    pd.DataFrame(erg_p1).to_csv(AUS / "ws3_p1_p2_korrelation_z.csv", index=False)

    # ---------------------------------------------------------------- P2b
    # the same computation with naive on raw abundance scale instead of z
    def korrelationslauf_roh(teilmenge_idx, nziehungen=NZIEH, seed=SEED):
        idx = list(teilmenge_idx)
        naiv = np.array([modulwert(zellen[i]["basis_roh"], gewichtet=False) for i in idx])
        dwt_amp = np.array([amplitude(zellen[i]["dwt"]) for i in idx])
        ok = np.isfinite(naiv) & np.isfinite(dwt_amp)
        naiv, dwt_amp = naiv[ok], dwt_amp[ok]
        idx_ok = [i for i, o in zip(idx, ok) if o]
        n = len(idx_ok)
        if n < 4:
            return {"n": n, "status": "too few donors"}
        rho = float(stats.spearmanr(naiv, dwt_amp).correlation)
        # null: permute the donor label of the dWT amplitude per study
        # (no algebraic coupling on the raw scale -- the null here is a
        # simple permutation null, no circularity artifact to reproduce,
        # since naive (raw) does not enter dWT (z))
        rng_l = np.random.default_rng(seed)
        by_study = {}
        for pos, i in enumerate(idx_ok):
            by_study.setdefault(zellen[i]["studie"], []).append(pos)
        null = np.empty(nziehungen)
        for it in range(nziehungen):
            dwt_perm = dwt_amp.copy()
            for studie, posliste in by_study.items():
                posliste_arr = np.array(posliste)
                perm = rng_l.permutation(posliste_arr)
                dwt_perm[posliste_arr] = dwt_amp[perm]
            null[it] = stats.spearmanr(naiv, dwt_perm).correlation
        m, s = float(np.mean(null)), float(np.std(null, ddof=1))
        z = (rho - m) / s if s > 0 else np.nan
        p = (1 + int((np.abs(null - m) >= abs(rho - m)).sum())) / (1 + len(null))
        return {"n": n, "status": "ok", "beobachtet": rho, "null_mittel": m,
                "null_sd": s, "z": float(z), "p": float(min(1.0, p))}

    erg_p2b = []
    for label, idx in (("primaer_n7_geeicht", idx_geeicht),
                        ("sensitivitaet_n14_alle", idx_alle)):
        r = korrelationslauf_roh(idx)
        r["variante"] = label
        r["skala_naiv"] = "roh (log2-Abundanz, unstandardisiert)"
        erg_p2b.append(r)
    pd.DataFrame(erg_p2b).to_csv(AUS / "ws3_p2b_korrelation_roh.csv", index=False)

    # ---------------------------------------------------------------- P3
    def achsen_test(idx):
        rows = []
        for i in idx:
            z = zellen[i]
            rows.append(dict(zelle=z["zelle"], achse=z["achse"],
                             achsendiff=T.loc[T.zelle == z["zelle"],
                                              "naiv_marker_achsendiff"].values[0]))
        d = pd.DataFrame(rows).dropna()
        return d

    d_achse = achsen_test(idx_geeicht)
    d_achse_alle = achsen_test(idx_alle)
    if d_achse.achse.nunique() > 1:
        auc_rows = []
        for label, dd in (("primaer_n7_geeicht", d_achse),
                           ("sensitivitaet_n14_alle", d_achse_alle)):
            if dd.achse.nunique() < 2 or len(dd) < 4:
                auc_rows.append({"variante": label, "status": "one axis or too few cells"})
                continue
            y = (dd.achse == "osteogen").astype(int).values
            x = dd.achsendiff.values
            auc = float(stats.mannwhitneyu(x[y == 1], x[y == 0],
                                            alternative="two-sided").statistic
                        / (max(1, (y == 1).sum()) * max(1, (y == 0).sum())))
            auc_rows.append({"variante": label, "status": "ok", "n": len(dd),
                            "n_osteogen": int((y == 1).sum()),
                            "n_chondrogen": int((y == 0).sum()), "auc": auc})
        pd.DataFrame(auc_rows).to_csv(AUS / "ws3_p3_achsenvorhersage.csv", index=False)
    else:
        pd.DataFrame([{"status": "only one axis among the calibrated cells -- "
                                  "P3 not testable for axis prediction"}]
                     ).to_csv(AUS / "ws3_p3_achsenvorhersage.csv", index=False)

    # ---------------------------------------------------------------- P4
    p4rows = []
    for spalte in ("naiv_modulewert_z", "naiv_modulewert_roh", "naiv_marker_achsendiff"):
        a = T.loc[T.bestanden, spalte].dropna()
        b = T.loc[~T.bestanden, spalte].dropna()
        if len(a) < 3 or len(b) < 3:
            p4rows.append({"variable": spalte, "status": "too few cells",
                          "n_bestanden": len(a), "n_durchgefallen": len(b)})
            continue
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        p4rows.append({"variable": spalte, "status": "ok",
                       "n_bestanden": len(a), "n_durchgefallen": len(b),
                       "median_bestanden": float(a.median()),
                       "median_durchgefallen": float(b.median()),
                       "U": float(u.statistic), "p": float(u.pvalue)})
    pd.DataFrame(p4rows).to_csv(AUS / "ws3_p4_eichung_naivunterschied.csv", index=False)

    print("done. rows T:", len(T), "| calibrated:", len(idx_geeicht))
    print(T[["zelle", "bestanden", "naiv_modulewert_z", "dwt_amplitude"]])


if __name__ == "__main__":
    haupt()

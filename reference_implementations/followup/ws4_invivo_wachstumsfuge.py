# -*- coding: utf-8 -*-
"""
ws4_invivo_wachstumsfuge.py -- WS4: anchoring the fixed 173-gene program in
real (in-vivo) skeletal development.

Data source: human limb single-cell atlas (Nature 2023,
s41586-023-06806-x), Pcw5.1-9.3, `celltype` carries the chondrogenesis axis
MesCond -> ChondroProg -> RestingChon -> ProlifChon -> PrehyperChon ->
HyperChon (growth-plate zoning in the condensation/chondrogenesis window).

Unit: `adj_sample` (stage x region x specimen) -- the finest biological
sample unit available in the dataset, NOT the cell.

Statistic: `reference_implementations/_module.py: kontrast` (contrast of two gene sets against
the same background-drawn null), no second implementation.

Positive control (mandatory, see BRIEFING/WS4): chondrogenic markers
(`_marker.py`) must rise along the axis, naive/proliferation markers fall --
otherwise the layer carries no finding about the module.
"""
from __future__ import annotations

import sys
import pathlib

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats

WURZEL = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(WURZEL / "reference_implementations"))
from _module import MODUL, kontrast          # noqa: E402
from _marker import OSTEOGEN, ADIPOGEN, MYOGEN, CHONDROGEN, NAIV  # noqa: E402

REF = WURZEL.parent / "Referenzdaten" / "Limb_Nature2023_s41586-023-06806-x"
H5AD = REF / "221114LimbCellranger3annotated.minimal.h5ad"
AUS_A = WURZEL / "reference_implementations" / "followup"
AUS_E = WURZEL / "derived_data" / "followup"
AUS_B = WURZEL / "Neu" / "Befunde"
for d in (AUS_A, AUS_E, AUS_B):
    d.mkdir(parents=True, exist_ok=True)

MIN_ZELLEN = 5          # per (zone, sample) -- otherwise no pseudobulk point
NZIEHUNGEN = 20000
SEED = 20260822

# The growth-plate/chondrogenesis axis, as annotated in the atlas.
ZONEN = ["MesCond", "ChondroProg", "RestingChon", "ProlifChon",
         "PrehyperChon", "HyperChon"]
RANG = {z: i for i, z in enumerate(ZONEN)}


def lade_obs(f):
    def cat(name):
        g = f["obs"][name]
        cats = np.array([c.decode() if isinstance(c, bytes) else c
                          for c in g["categories"][:]])
        codes = g["codes"][:]
        return pd.Series(np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA"))
    obs = pd.DataFrame({
        "celltype": cat("celltype"),
        "adj_sample": cat("adj_sample"),
        "adj_stage": cat("adj_stage"),
    })
    return obs


def main():
    print("loading h5ad (backed) ...", flush=True)
    f = h5py.File(H5AD, "r")
    obs = lade_obs(f)
    var_index = np.array([g.decode() if isinstance(g, bytes) else g
                           for g in f["var"]["_index"][:]])
    symbol_zu_idx = {}
    for i, s in enumerate(var_index):
        symbol_zu_idx.setdefault(s, i)  # first occurrence

    # -- gene lists -----------------------------------------------------
    modul_up = MODUL.loc[MODUL.ri > 0, "symbol"].tolist()
    modul_dn = MODUL.loc[MODUL.ri < 0, "symbol"].tolist()
    modul_alle = MODUL["symbol"].tolist()

    panels = pd.read_csv(WURZEL / "derived_data" / "M_humangenetik" / "panels.csv")
    pa309 = panels.loc[panels.panel == "PA309", "symbol"].dropna().unique().tolist()
    noso = panels.loc[panels.panel == "NOSO", "symbol"].dropna().unique().tolist()

    benoetigt = sorted(set(modul_alle) | set(OSTEOGEN) | set(ADIPOGEN)
                        | set(MYOGEN) | set(CHONDROGEN) | set(NAIV)
                        | set(pa309) | set(noso))
    benoetigt = [g for g in benoetigt if g in symbol_zu_idx]
    print(f"required genes found in the atlas: {len(benoetigt)}", flush=True)

    # background for the null: 4000 random genes expressed in the atlas
    # (mean raw count > 0 over all cells), plus all required ones.
    rng = np.random.default_rng(SEED)
    n_genes = len(var_index)
    hg_idx = rng.choice(n_genes, size=min(4000, n_genes), replace=False)
    hg_symbole = var_index[hg_idx]
    alle_symbole = sorted(set(benoetigt) | set(hg_symbole))
    alle_symbole = [g for g in alle_symbole if g in symbol_zu_idx]
    spalten_idx = np.array([symbol_zu_idx[g] for g in alle_symbole])
    print(f"loading column block: {len(spalten_idx)} genes x {len(obs)} cells ...",
          flush=True)

    # -- load the CSR completely, then subset columns (faster than H5 fancy indexing)
    indptr = f["X"]["indptr"][:]
    indices = f["X"]["indices"][:]
    data = f["X"]["data"][:]
    X = sparse.csr_matrix((data, indices, indptr),
                          shape=(len(obs), n_genes))
    del data, indices, indptr
    Xs = X[:, spalten_idx].tocsc()
    del X
    print("column block loaded.", flush=True)

    # -- pseudobulk per (zone, sample) ------------------------------------
    zeilen = []
    matrizen = {}   # (zone) -> DataFrame Probe x Gen
    for zone in ZONEN:
        maske_zone = (obs["celltype"].values == zone)
        proben = obs.loc[maske_zone, "adj_sample"].unique()
        werte = {}
        for p in proben:
            m = maske_zone & (obs["adj_sample"].values == p)
            n = int(m.sum())
            if n < MIN_ZELLEN:
                continue
            mw = np.asarray(Xs[m, :].mean(axis=0)).ravel()
            werte[p] = mw
            zeilen.append({"zone": zone, "probe": p, "n_zellen": n})
        if werte:
            df = pd.DataFrame.from_dict(werte, orient="index",
                                         columns=alle_symbole)
            matrizen[zone] = df
    proben_uebersicht = pd.DataFrame(zeilen)
    proben_uebersicht.to_csv(AUS_E / "ws4_proben_je_zone.csv", index=False)
    print(proben_uebersicht.groupby("zone")["probe"].count())

    f.close()

    # -- positive control: chondrogenic markers up, naive markers down ----
    kontrolle_zeilen = []
    for zone in ZONEN:
        if zone not in matrizen:
            continue
        df = matrizen[zone]
        for probe, row in df.iterrows():
            r = kontrast(row, CHONDROGEN, NAIV,
                         nziehungen=NZIEHUNGEN, seed=SEED)
            r.update({"zone": zone, "probe": probe, "vergleich": "chondrogen_vs_naiv"})
            kontrolle_zeilen.append(r)
            r2 = kontrast(row, OSTEOGEN, ADIPOGEN,
                          nziehungen=NZIEHUNGEN, seed=SEED)
            r2.update({"zone": zone, "probe": probe, "vergleich": "osteogen_vs_adipogen"})
            kontrolle_zeilen.append(r2)
    kontrolle = pd.DataFrame(kontrolle_zeilen)
    kontrolle.to_csv(AUS_E / "ws4_positivkontrolle_je_probe.csv", index=False)

    def trendtest(df, wertspalte="z"):
        """Jonckheere-like trend: Spearman of zone rank against the sample
        value, against a null that permutes the sample values (zone
        assignment redrawn at random, sample count per zone preserved)."""
        sub = df.dropna(subset=[wertspalte]).copy()
        sub["rang"] = sub["zone"].map(RANG)
        if sub["rang"].nunique() < 3 or len(sub) < 6:
            return {"status": "too few zones/samples", "n": len(sub)}
        beob = stats.spearmanr(sub["rang"], sub[wertspalte]).statistic
        rng2 = np.random.default_rng(SEED)
        werte = sub[wertspalte].to_numpy()
        raenge = sub["rang"].to_numpy()
        null = np.empty(NZIEHUNGEN)
        for i in range(NZIEHUNGEN):
            perm = rng2.permutation(werte)
            null[i] = stats.spearmanr(raenge, perm).statistic
        m, s = float(np.nanmean(null)), float(np.nanstd(null, ddof=1))
        z = (beob - m) / s if s > 0 else np.nan
        p_ein = (1 + int((null >= beob).sum())) / (1 + len(null))
        return {"status": "ok", "n": len(sub), "spearman_rho": float(beob),
                "null_mittel": m, "null_sd": s, "z": float(z),
                "p": float(min(1.0, 2 * p_ein)), "mde80": m + 2.8 * s}

    kontrolle_chondro = kontrolle[kontrolle.vergleich == "chondrogen_vs_naiv"]
    trend_kontrolle = trendtest(kontrolle_chondro, "kontrast")
    print("positive control (chondrogenic-naive, trend across zones):", trend_kontrolle)

    besteht_eichung = (trend_kontrolle.get("status") == "ok"
                        and trend_kontrolle.get("z", np.nan) > 2.8)

    ausgabe = {"positivkontrolle_trend": trend_kontrolle,
               "besteht_eichung": bool(besteht_eichung)}

    if not besteht_eichung:
        print("CALIBRATION NOT PASSED -- no further numbers on the module.")
        pd.Series(ausgabe).to_json(AUS_E / "ws4_status.json")
        return ausgabe

    # -- P1: module value along the axis ----------------------------------
    modul_zeilen = []
    for zone in ZONEN:
        if zone not in matrizen:
            continue
        df = matrizen[zone]
        for probe, row in df.iterrows():
            r = kontrast(row, modul_up, modul_dn, nziehungen=NZIEHUNGEN, seed=SEED)
            r.update({"zone": zone, "probe": probe})
            modul_zeilen.append(r)
    modul_df = pd.DataFrame(modul_zeilen)
    modul_df.to_csv(AUS_E / "ws4_modulwert_je_probe.csv", index=False)
    trend_module = trendtest(modul_df, "kontrast")
    print("P1 Modultrend:", trend_module)
    ausgabe["p1_moduletrend"] = trend_module

    # -- P2: which module genes are recoverable in vivo? -------------------
    # per module gene: concordance of sign(HyperChon pseudobulk - MesCond) with ri,
    # over all sample pairs (zone early vs late, same sample when possible,
    # otherwise zone means).
    frueh, spaet = "MesCond", "HyperChon"
    p2 = None
    if frueh in matrizen and spaet in matrizen:
        mittel_frueh = matrizen[frueh].mean(axis=0)
        mittel_spaet = matrizen[spaet].mean(axis=0)
        delta = (mittel_spaet - mittel_frueh)
        rows = []
        for _, g in MODUL.iterrows():
            sym = g["symbol"]
            if sym in delta.index:
                d = delta[sym]
                rows.append({"symbol": sym, "ensembl": g["ensembl"], "ri": g["ri"],
                             "delta_hyper_minus_mes": float(d),
                             "konkordant": bool(np.sign(d) == g["ri"]) if d != 0 else None})
        p2 = pd.DataFrame(rows)
        p2.to_csv(AUS_E / "ws4_p2_gen_konkordanz.csv", index=False)
        n_konk = int(p2["konkordant"].sum())
        n_tot = int(p2["konkordant"].notna().sum())
        print(f"P2: {n_konk}/{n_tot} module genes concordant in vivo (Mes->Hyper) "
              f"with the predicted direction.")
        ausgabe["p2_anteil_konkordant"] = n_konk / n_tot if n_tot else np.nan
        ausgabe["p2_n"] = n_tot

    # -- P3: skeletal dysplasia genes (PA309/NOSO) in the same zone? -------
    p3_zeilen = []
    pa309_vorh = [g for g in pa309 if g in alle_symbole]
    for zone in ZONEN:
        if zone not in matrizen:
            continue
        df = matrizen[zone]
        for probe, row in df.iterrows():
            r = kontrast(row, pa309_vorh, modul_alle,
                         nziehungen=NZIEHUNGEN, seed=SEED)
            r.update({"zone": zone, "probe": probe, "vergleich": "PA309_vs_Modul"})
            p3_zeilen.append(r)
    p3 = pd.DataFrame(p3_zeilen)
    p3.to_csv(AUS_E / "ws4_p3_panel_vs_modul.csv", index=False)
    trend_p3 = trendtest(p3, "kontrast")
    print("P3 trend PA309 vs module across zones:", trend_p3)
    ausgabe["p3_panel_vs_module_trend"] = trend_p3

    import json
    with open(AUS_E / "ws4_status.json", "w") as fh:
        json.dump(ausgabe, fh, indent=2, default=float)
    return ausgabe


if __name__ == "__main__":
    main()

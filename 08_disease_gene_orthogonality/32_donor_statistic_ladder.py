# -*- coding: utf-8 -*-
"""
32_donor_statistic_ladder.py -- phase M-D: the preregistered statistical ladder.

Preregistration `PRAEREG_M_D.md` §7-§9 including addendum 1. Computed
exclusively on cells with a **passed built-in calibration**
(`31_donor_cells_build_calibrate.py`); failed cells carry no number.

  S1 (primary)  mean pairwise Spearman correlation between different
                donors, per cell and pooled
  S2            share of the first principal component (uncentered)
  S3a           directed sign concordance against `ri`
  S3b           direction-free agreement

Three sets of vectors are computed separately:

  Program               `dWT` of all calibrated cells (E2 does not matter
                        here -- the differentiation arm is untouched)
  Lesion response       `iv` of the calibrated cells with **E2 fulfilled**
  Engineering response  `iv` of the calibrated cells with E2 not fulfilled --
                        reported separately, never mixed with the lesion
                        response into one number (addendum 1 b)

All with baseline-stratified null, 20 000 draws, seed 20260822;
synthesis exclusively against the donor-flip null, never as an average;
leave-one-out computation complete.

Output: derived_data/M_donoren/{statistik.csv, je_zelle.csv, synthese.csv,
         auslassung.csv, 54c_log.txt}
"""
from __future__ import annotations

import pathlib
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import (ERGEBNISSE, MODUL, NZIEHUNGEN, SEED_D,  # noqa: E402
                    leiter, synthese_flip)

AUS = ERGEBNISSE / "M_donoren"
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def laesionssatz() -> pd.DataFrame:
    """The equal-size 173-gene set of highest `iv` consistency.
    Taken unchanged from phase M-B, not rebuilt."""
    return pd.read_csv(ERGEBNISSE / "M_patienten" / "laesionssatz_173.csv")


def sammle(daten: dict, zellen: list[str], groesse: str
           ) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """Cells x gene over the genes common to all participating studies."""
    teile, basen = [], []
    for gse, d in daten.items():
        sp = [c for c in d[groesse].columns if c in zellen]
        if not sp:
            continue
        teile.append(d[groesse][sp].T)
        basen.append(d["basis"])
    X = pd.concat(teile, axis=0, join="inner")
    basis = pd.concat(basen, axis=1).mean(axis=1).reindex(X.columns).dropna()
    X = X[basis.index]
    X = X.loc[:, X.notna().all(axis=0)]
    return X, basis.reindex(X.columns), X


def rechne(name: str, X: pd.DataFrame, satz: pd.DataFrame, basis: pd.Series,
           spender: list, nziehungen: int) -> tuple[dict, pd.DataFrame]:
    ri = satz.set_index("symbol").ri
    gene = [g for g in satz.symbol if g in X.columns]
    r = leiter(X, gene, ri, X, basis, spender=spender,
               nziehungen=nziehungen, seed=SEED_D)
    if r.get("status") != "ok":
        log("   %-22s %s" % (name, r.get("status")))
        return r, pd.DataFrame()
    log("   %-22s %d cells, %d donor pairs, %d genes" %
        (name, r["n_zellen"], r["n_paare"], r["n_gene"]))
    for s in ("S1", "S2", "S3a", "S3b"):
        log("      %-4s %+7.4f | null %+7.4f +- %6.4f | z %+6.2f | "
            "p %8.4g | MDE80 %+7.4f %s"
            % (s, r[f"{s}_beobachtet"], r[f"{s}_null_mittel"],
               r[f"{s}_null_sd"], r[f"{s}_z"], r[f"{s}_p"], r[f"{s}_mde80"],
               "ABOVE MDE80" if r[f"{s}_beobachtet"] > r[f"{s}_mde80"] else ""))
    j = r.pop("je_zelle")
    j.insert(0, "groesse", name)
    return r, j


def main() -> None:
    log("=" * 78)
    log("Phase M-D  --  the preregistered statistical ladder")
    log("PRAEREG_M_D.md §7-§9, seed %d, %d draws" % (SEED_D, NZIEHUNGEN))
    log("=" * 78)

    with open(AUS / "zellen.pkl", "rb") as f:
        daten = pickle.load(f)
    E = pd.read_csv(AUS / "eichung.csv")
    geeicht = E[E.bestanden]
    log("\nCalibrated cells: %d of %d (%d studies)"
        % (len(geeicht), len(E), geeicht.studie.nunique()))
    for _, r in geeicht.iterrows():
        log("   %-24s %-11s donor %-16s E2 %s  (%s)"
            % (r.zelle, r.achse, r.spender, "yes" if r.e2 else "no", r.herkunft))

    programm = MODUL[["symbol", "ri"]].copy()
    laesion = laesionssatz()[["symbol", "ri"]]

    saetze = [
        ("Programm", "dwt", list(geeicht.zelle), programm),
        ("Laesionsantwort", "iv", list(geeicht[geeicht.e2].zelle), laesion),
        ("Engineering-Antwort", "iv", list(geeicht[~geeicht.e2].zelle), laesion),
    ]

    zeilen, je_zelle = [], []
    matrizen = {}
    for name, groesse, zellen, satz in saetze:
        log("\n" + "-" * 78)
        log("%s  --  %d calibrated cell(s)" % (name, len(zellen)))
        if len(zellen) < 2:
            log("   fewer than two cells -- no between-donor statistic "
                "possible. NO number is reported (preregistration §9).")
            zeilen.append(dict(groesse=name, n_zellen=len(zellen),
                               status="too few calibrated cells"))
            continue
        X, basis, _ = sammle(daten, zellen, groesse)
        sp = list(geeicht.set_index("zelle").spender.reindex(X.index))
        r, j = rechne(name, X, satz, basis, sp, NZIEHUNGEN)
        r.update(groesse=name, n_spender=len(set(sp)),
                 studien="|".join(sorted(set(
                     geeicht.set_index("zelle").studie.reindex(X.index)))))
        zeilen.append(r)
        if len(j):
            je_zelle.append(j)
            # per-study matrices for the donor-flip synthesis
            ri = satz.set_index("symbol").ri
            gene = [g for g in satz.symbol if g in X.columns]
            for gse in sorted(set(geeicht.set_index("zelle").studie.reindex(X.index))):
                zl = [c for c in X.index
                      if geeicht.set_index("zelle").studie.get(c) == gse]
                matrizen.setdefault(name, {})[gse] = (
                    np.nan_to_num(X.loc[zl, gene].to_numpy(float)),
                    ri.reindex(gene).to_numpy(float),
                    list(geeicht.set_index("zelle").spender.reindex(zl)))

    S = pd.DataFrame(zeilen)
    S.to_csv(AUS / "statistik.csv", index=False)
    if je_zelle:
        J = pd.concat(je_zelle, ignore_index=True)
        J = J.merge(geeicht[["zelle", "studie", "spender", "achse", "e2",
                             "herkunft", "n_proben"]], on="zelle", how="left")
        J.to_csv(AUS / "je_zelle.csv", index=False)
        log("\n" + "-" * 78)
        log("Per cell -- S1 against its own MDE80 (primary statistic):")
        for _, r in J[J.groesse == "Programm"].iterrows():
            log("   %-24s S1 %+6.3f | null %+6.3f +- %5.3f | z %+6.2f | "
                "MDE80 %+6.3f -> %s"
                % (r.zelle, r.S1, r.S1_null_mittel, r.S1_null_sd, r.S1_z,
                   r.S1_mde80, "ABOVE" if r.S1_ueber_mde80 else "below"))

    # ---- study synthesis against the donor-flip null ----------------------
    syn_zeilen, syn_studien = [], []
    for name, m in matrizen.items():
        if len(m) < 2:
            log("\nSynthesis %s: only one study -- no synthesis." % name)
            continue
        log("\n" + "-" * 78)
        log("Study synthesis %s (donor-flip null, never an average)" % name)
        for kz in ("S1", "S2", "S3a", "S3b"):
            r = synthese_flip(m, kennzahl=kz, nziehungen=10000, seed=SEED_D)
            js = r.pop("je_studie")
            js.insert(0, "groesse", name)
            syn_studien.append(js)
            r["groesse"] = name
            syn_zeilen.append(r)
            wirksam = int(js.beobachtet.notna().sum())
            if r.get("entartet"):
                log("   %-4s observed %+7.4f | null DEGENERATE (spread 0) -- "
                    "the statistic is invariant under donor flips. "
                    "NO number reported." % (kz, r["beobachtet"]))
            else:
                log("   %-4s observed %+7.4f | null %+7.4f +- %6.4f | z %+6.2f "
                    "| p %8.4g | MDE80 %+7.4f | studies above MDE80: %d/%d "
                    "| studies with a defined statistic: %d"
                    % (kz, r["beobachtet"], r["null_mittel"], r["null_sd"],
                       r["z"], r["p"], r["mde80"], int(js.ueber_mde80.sum()),
                       len(js), wirksam))
            r["n_studien_wirksam"] = wirksam
    if syn_zeilen:
        pd.DataFrame(syn_zeilen).to_csv(AUS / "synthese.csv", index=False)
        pd.concat(syn_studien, ignore_index=True).to_csv(
            AUS / "synthese_je_studie.csv", index=False)

    # ---- leave-one-out: each cell dropped once ----------------------------
    log("\n" + "-" * 78)
    log("Leave-one-out -- each cell dropped once (S1, program)")
    aus_zeilen = []
    zellen = list(geeicht.zelle)
    X, basis, _ = sammle(daten, zellen, "dwt")
    ri = programm.set_index("symbol").ri
    gene = [g for g in programm.symbol if g in X.columns]
    for weg in [None] + zellen:
        bleib = [c for c in X.index if c != weg]
        if len(bleib) < 2:
            continue
        sp = list(geeicht.set_index("zelle").spender.reindex(bleib))
        r = leiter(X.loc[bleib], gene, ri, X.loc[bleib], basis, spender=sp,
                   nziehungen=2000, seed=SEED_D)
        r.pop("je_zelle", None)
        aus_zeilen.append(dict(ohne=weg or "-- (all)", n_zellen=len(bleib),
                               S1=r["S1_beobachtet"], S1_z=r["S1_z"],
                               S1_p=r["S1_p"], S1_mde80=r["S1_mde80"],
                               S2_z=r["S2_z"], S3a_z=r["S3a_z"],
                               S3b_z=r["S3b_z"],
                               nziehungen=2000))
        log("   without %-20s n %2d | S1 %+6.3f | z %+6.2f | p %8.4g"
            % (aus_zeilen[-1]["ohne"], len(bleib), r["S1_beobachtet"],
               r["S1_z"], r["S1_p"]))
    pd.DataFrame(aus_zeilen).to_csv(AUS / "auslassung.csv", index=False)

    # ---- separate analysis: cells with few control samples ----------------
    log("\n" + "-" * 78)
    log("Separate analysis: cells with fewer than 3 control samples")
    wenig = list(E.set_index("zelle").n_proben[
        E.set_index("zelle").n_proben < 6].index)
    viele = [c for c in zellen if c not in wenig]
    log("   cells with < 3 control samples (< 6 samples in the 2x2): %s"
        % (", ".join(w for w in wenig if w in zellen) or "none"))
    if len(viele) >= 2 and len(viele) < len(zellen):
        sp = list(geeicht.set_index("zelle").spender.reindex(viele))
        r = leiter(X.loc[viele], gene, ri, X.loc[viele], basis, spender=sp,
                   nziehungen=2000, seed=SEED_D)
        log("   only cells with >= 3 control samples: S1 %+6.3f | z %+6.2f | "
            "p %8.4g" % (r["S1_beobachtet"], r["S1_z"], r["S1_p"]))

    log("\n" + "=" * 78)
    (AUS / "54c_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()

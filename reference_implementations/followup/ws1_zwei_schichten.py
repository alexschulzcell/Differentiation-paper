# -*- coding: utf-8 -*-
"""
ws1_zwei_schichten.py -- work package WS1: "two layers".

Hypothesis: skeletal-dysplasia/short-stature genes are constitutively highly
expressed, dose-sensitive secretion/matrix infrastructure -- not
transcriptionally dynamic differentiation genes.

Design word-for-word modeled on
`reference_implementations/52_human_genetics_anchor.py` (grid from
expression x length deciles, matched draw null, MDE80 = mu + 2.8*sd),
here generalized to CONTINUOUS quantities per gene (not only count-in-panel):
mean absolute expression (rank), |dWT_med|, LOEUF.

Status: EXPLORATORY. Nothing here is preregistered.

Output: derived_data/followup/ws1_*.csv
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

WURZEL = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(WURZEL / "reference_implementations"))
from _marker import CHONDROGEN, OSTEOGEN  # noqa: E402
from _module import ERGEBNISSE, SEED, NZIEHUNGEN  # noqa: E402

NEU = WURZEL / "Neu"
AUS = NEU / "derived_data"
AUS.mkdir(parents=True, exist_ok=True)
NDEZ = 10
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------------------
def lade_genkarte() -> pd.DataFrame:
    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    return K


def lade_exonlaengen() -> pd.Series:
    d = pd.read_csv(ERGEBNISSE / "M_humangenetik" / "_exonlaengen.csv")
    return pd.Series(d.laenge.values, index=d.ensembl)


def lade_expr18() -> pd.DataFrame:
    """Per gene: median expression rank over the 18 datasets (rank WITHIN
    each dataset, then median over datasets -- no arithmetic mean over
    cohorts, rule 4), plus the number of datasets in which the gene was
    measurable at all (rowVar>0 filter resp. rlog filter already applied in
    the loader)."""
    frames = []
    for i in range(1, 19):
        f = AUS / "ws1_expr" / f"punkt_{i:02d}.csv"
        d = pd.read_csv(f)
        frames.append(d[["ensembl", "expr_rank"]])
    all_d = pd.concat(frames, ignore_index=True)
    g = all_d.groupby("ensembl")["expr_rank"]
    out = g.median().rename("expr_rank_med").to_frame()
    out["n_gemessen"] = g.size()
    return out.reset_index()


def lade_constraint() -> pd.DataFrame:
    d = pd.read_csv(WURZEL / "data_raw" / "_referenz" / "gnomad_constraint" /
                     "constraint_reduziert.csv")
    d = d.rename(columns={"ensembl": "ensembl_gnomad"})
    # gnomAD v2.1.1 gene_id is partly an outdated Ensembl ID -- merging via
    # symbol is more robust (symbols are known from the gene map).
    d = d[["symbol", "loeuf", "pli", "mis_z"]].drop_duplicates(subset="symbol")
    return d


def sym2ens(symbole, K: pd.DataFrame) -> set:
    m = dict(zip(K.symbol.astype(str), K.ensembl))
    return {m[s] for s in symbole if s in m}


def raster(G: pd.DataFrame) -> pd.DataFrame:
    """Background grid: decile(basis_med) x decile(length) -- AS IN
    52_human_genetics_anchor.py, so that the same positive control applies."""
    d = G[["ensembl", "basis_med", "laenge"]].dropna().copy()
    d["e_dez"] = pd.qcut(d.basis_med.rank(method="first"), NDEZ, labels=False)
    d["l_dez"] = pd.qcut(d.laenge.rank(method="first"), NDEZ, labels=False)
    d["zelle"] = d.e_dez * NDEZ + d.l_dez
    return d.set_index("ensembl")


def gezogene_mittelwertnull(zellen_map: dict, z_satz: np.ndarray,
                             werte: pd.Series, nzieh: int, seed: int) -> np.ndarray:
    """Draws for each set gene a background gene from the same grid cell
    (with replacement across draws), computes the mean of the target
    quantity `werte` for the draw. Returns the null distribution of the
    means (length nzieh)."""
    rng = np.random.default_rng(seed)
    null = np.empty(nzieh)
    wv = {i: v for i, v in enumerate(werte.values)}
    for i in range(nzieh):
        idx = np.array([rng.choice(zellen_map[z]) for z in z_satz])
        null[i] = werte.values[idx].mean()
    return null


def stetiger_test(satz: set, HG: pd.DataFrame, spalte: str,
                   nzieh: int = NZIEHUNGEN, seed: int = SEED,
                   min_n: int = 8) -> dict:
    """Tests whether the mean of `spalte` for `satz` deviates from the matched
    background null. `HG` must contain `spalte` and `zelle` and must have no
    missing values in `spalte` (filter beforehand)."""
    d = HG.dropna(subset=[spalte])
    s = [g for g in satz if g in d.index]
    n = len(s)
    if n < min_n:
        return {"status": "too small", "n_satz": n}

    zellen = {}
    for z, idx in d.groupby("zelle").groups.items():
        zellen[z] = np.array([d.index.get_loc(g) for g in idx])
    z_satz = d.loc[s, "zelle"].values
    werte = d[spalte]

    beob = float(werte.loc[s].mean())
    null = gezogene_mittelwertnull(zellen, z_satz, werte, nzieh, seed)
    mu, sd = float(null.mean()), float(null.std(ddof=1))
    z = (beob - mu) / sd if sd > 0 else np.nan
    p_ein = (1 + int((null >= beob).sum())) / (1 + nzieh) if beob >= mu else \
            (1 + int((null <= beob).sum())) / (1 + nzieh)
    p_emp = float(min(1.0, 2 * p_ein))
    return {"status": "ok", "n_satz": n, "spalte": spalte,
            "beobachtet": beob, "null_mittel": mu, "null_sd": sd,
            "z": float(z), "p": p_emp,
            "mde80_diff": 2.8 * sd,
            "diff": beob - mu}


def main() -> None:
    log("=" * 78)
    log("WS1 -- two layers: constitutive infrastructure vs induced program")
    log("EXPLORATORY -- nothing here is preregistered. Seed %d, %d draws."
        % (SEED, NZIEHUNGEN))
    log("=" * 78)

    K = lade_genkarte()
    L = lade_exonlaengen()
    E = lade_expr18()
    C = lade_constraint()

    G = K.copy()
    G["laenge"] = G.ensembl.map(L)
    G = G.merge(E, on="ensembl", how="left")
    G = G.merge(C, on="symbol", how="left")
    G["dWT_abs"] = G.dWT_med.abs()
    G.to_csv(AUS / "ws1_genkarte_erweitert.csv", index=False)
    log("\nGene map extended: %d genes, of which %d with length, %d with "
        "expression rank (>=1 of the 18 points), %d with gnomAD constraint."
        % (len(G), G.laenge.notna().sum(), G.expr_rank_med.notna().sum(),
           G.loeuf.notna().sum()))

    HG = raster(G)
    HG = HG.join(G.set_index("ensembl")[["expr_rank_med", "dWT_abs", "loeuf",
                                          "n_gemessen"]])
    log("Background grid: %d genes (expression AND length known), %dx%d cells"
        % (len(HG), NDEZ, NDEZ))

    P = pd.read_csv(ERGEBNISSE / "M_humangenetik" / "panels.csv")
    P = P[P.ensembl.notna()]
    panels = {k: set(v) for k, v in P.groupby("panel").ensembl}

    programm = set(pd.read_csv(WURZEL / "derived_data" / "reference_tables" /
                                "S5_konvergente_gene.csv").rename(
                                columns={"gen": "ensembl"}).ensembl)
    laesion = set(pd.read_csv(ERGEBNISSE / "M_patienten" /
                              "laesionssatz_173.csv").ensembl)
    zellzyklus = set(pd.read_csv(AUS / "ws1_zellzyklus.csv").ensembl)
    marker = sym2ens(OSTEOGEN + CHONDROGEN, K)

    saetze = {
        "KLEIN": panels.get("KLEIN", set()),
        "KLEIN_BREIT": panels.get("KLEIN_BREIT", set()),
        "NOSO": panels.get("NOSO", set()),
        "NOSO_BREIT": panels.get("NOSO_BREIT", set()),
        "PA309": panels.get("PA309", set()),
        "PA1471": panels.get("PA1471", set()),
        "GWAS": panels.get("GWAS", set()),
        "Programm(173)": programm,
        "Laesionsantwort(173)": laesion,
        "Zellzyklus(Negativkontrolle)": zellzyklus,
        "Linienmarker(Positivkontrolle)": marker,
    }
    log("\nSet sizes in the background: " + ", ".join(
        "%s %d" % (k, len(v & set(HG.index))) for k, v in saetze.items()))

    # -------------------------------------------------------------- P1 + P2
    log("\n" + "-" * 78)
    log("P1 -- absolute expression (rank over 18 datasets, median) per set "
        "against matched null")
    log("P2 -- |dWT_med| (dynamics, z scale from the gene map) per set against "
        "matched null")
    zeilen = []
    for sn, satz in saetze.items():
        for spalte, tag in (("expr_rank_med", "P1_expr"), ("dWT_abs", "P2_dyn"),
                             ("loeuf", "P3_loeuf")):
            r = stetiger_test(satz, HG, spalte)
            r.update(satz=sn, tag=tag)
            zeilen.append(r)
            if r["status"] != "ok":
                log("%-32s %-9s %s" % (sn, tag, r["status"]))
                continue
            log("%-32s %-9s n=%4d beob=%7.4f null=%7.4f diff=%+7.4f z=%+6.2f "
                "p=%8.4g MDE80(diff)=%.4f"
                % (sn, tag, r["n_satz"], r["beobachtet"], r["null_mittel"],
                   r["diff"], r["z"], r["p"], r["mde80_diff"]))
    R = pd.DataFrame(zeilen)
    R.to_csv(AUS / "ws1_stetige_tests.csv", index=False)

    # -------------------------------------------------------------- P4
    log("\n" + "-" * 78)
    log("P4 -- share of panel genes measured at all in the 18 datasets "
        "(trivial alternative explanation: not expressed)")
    zeilen4 = []
    for sn, satz in saetze.items():
        s = [g for g in satz if g in set(G.ensembl)]
        n = len(s)
        if n == 0:
            continue
        sub = G[G.ensembl.isin(s)]
        gemessen = sub.n_gemessen.fillna(0)
        anteil_irgendwo = float((gemessen > 0).mean())
        anteil_mehrheit = float((gemessen >= 9).mean())  # in >=9 of 18 points
        zeilen4.append(dict(satz=sn, n=n,
                             anteil_mind1_datensatz=anteil_irgendwo,
                             anteil_mind_halbe_datensaetze=anteil_mehrheit,
                             median_n_gemessen=float(gemessen.median())))
        log("%-32s n=%4d in>=1 point: %5.1f%%  in>=9 points: %5.1f%%  "
            "median n_gemessen=%.0f/18"
            % (sn, n, 100 * anteil_irgendwo, 100 * anteil_mehrheit,
               gemessen.median()))
    pd.DataFrame(zeilen4).to_csv(AUS / "ws1_p4_messbarkeit.csv", index=False)

    log("=" * 78)
    (AUS / "ws1_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()

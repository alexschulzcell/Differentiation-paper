# -*- coding: utf-8 -*-
"""
52_human_genetics_anchor.py -- phase M-A: the human-genetic anchor.

Preregistration `PRAEREG_M_A.md` (2026-08-21) including addendum 1, both
dated before the first statistic of this phase.

The question is whether the 173 convergent program genes are enriched for
skeletal dysplasia, short-stature or body-height genes -- with considerably
better power than in the old version (357 and 50 genes, detection limit
~OR 1.6). The counter-check is the equal-size lesion response set.

The background is the gene pool of the computation, **expression- and
length-matched**: the null draws for each set gene a background gene from
the same cell of the grid (decile of baseline expression x decile of union
exon length). Without this matching the OR is not interpretable -- disease
genes are longer and more highly expressed than average.

The matching is **not an adjustment of the target quantity**: panel
membership is an external annotation, not a property of these datasets.
The guard "matching also matches the target quantity" concerns covariates
of the baseline on the z scale and does not apply here.

Output: derived_data/M_humangenetik/{anker,anker_power,eichung_A}.csv
"""
from __future__ import annotations

import gzip
import pathlib
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _marker import CHONDROGEN, OSTEOGEN  # noqa: E402
from _module import DATEN, ERGEBNISSE, MODUL, NZIEHUNGEN, SEED  # noqa: E402

AUS = ERGEBNISSE / "M_humangenetik"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []
NDEZ = 10
BONFERRONI = 14          # 7 panels x 2 gene sets, addendum 1 (c)


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------------------
def exonlaengen() -> pd.Series:
    """Union exon length per gene from the project's Gencode reference."""
    zwi = AUS / "_exonlaengen.csv"
    if zwi.exists():
        d = pd.read_csv(zwi)
        return pd.Series(d.laenge.values, index=d.ensembl)
    gtf = DATEN / "_referenz" / "gencode.v46.annotation.hg38.gtf.gz"
    stuecke: dict[str, list] = defaultdict(list)
    with gzip.open(gtf, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln[0] == "#":
                continue
            t = ln.split("\t")
            if len(t) < 9 or t[2] != "exon":
                continue
            g = re.search(r'gene_id "([^".]+)', t[8])
            if g:
                stuecke[g.group(1)].append((int(t[3]), int(t[4])))
    aus = {}
    for g, iv in stuecke.items():
        iv.sort()
        ges, ende = 0, -1
        for a, e in iv:
            a = max(a, ende + 1)
            if e >= a:
                ges += e - a + 1
                ende = max(ende, e)
        aus[g] = ges
    s = pd.Series(aus, name="laenge")
    s.rename_axis("ensembl").reset_index().to_csv(zwi, index=False)
    return s


def raster(K: pd.DataFrame) -> pd.DataFrame:
    """Background pool with expression and length deciles."""
    L = exonlaengen()
    d = K[["ensembl", "basis_med"]].dropna().copy()
    d["laenge"] = d.ensembl.map(L)
    d = d.dropna()
    d["e_dez"] = pd.qcut(d.basis_med.rank(method="first"), NDEZ, labels=False)
    d["l_dez"] = pd.qcut(d.laenge.rank(method="first"), NDEZ, labels=False)
    d["zelle"] = d.e_dez * NDEZ + d.l_dez
    return d.set_index("ensembl")


def gematcht(satz: set, panel: set, HG: pd.DataFrame,
             nzieh: int = NZIEHUNGEN, seed: int = SEED) -> dict:
    """Overlap set x panel against the matched background null."""
    s = [g for g in satz if g in HG.index]
    p = {g for g in panel if g in HG.index}
    n, m, N = len(s), len(p), len(HG)
    # Minimum size as in `_module.konkordanz`: below 8 measurable genes no
    # computation. The canonical marker sets are small (30 symbols, some of
    # which are not expressed in the computation pool) -- a higher bound
    # would let the positive control fail for purely technical reasons.
    if n < 8 or m < 10:
        return {"status": "too small", "n_satz": n, "n_panel": m}
    beob = len(set(s) & p)

    zellen = defaultdict(list)
    ist_panel = HG.index.isin(p)
    for i, z in enumerate(HG.zelle.values):
        zellen[z].append(i)
    zellen = {z: np.array(v) for z, v in zellen.items()}
    z_satz = HG.zelle.reindex(s).values

    rng = np.random.default_rng(seed)
    null = np.empty(nzieh)
    for i in range(nzieh):
        idx = np.array([rng.choice(zellen[z]) for z in z_satz])
        null[i] = ist_panel[idx].sum()
    mu, sd = float(null.mean()), float(null.std(ddof=1))
    z = (beob - mu) / sd if sd > 0 else np.nan
    p_ein = (1 + int((null >= beob).sum())) / (1 + nzieh)
    p_emp = float(min(1.0, 2 * min(p_ein, 1 - p_ein + 1 / (1 + nzieh))))

    def als_or(k: float) -> float:
        """Count -> OR against the expectation of the matched null."""
        if k <= 0 or k >= n or mu <= 0 or mu >= n:
            return np.nan
        return (k / (n - k)) / (mu / (n - mu))

    # raw 2x2 OR, as in f4_krankheitsanreicherung.csv
    a, b = beob, n - beob
    c, d_ = m - beob, N - n - (m - beob)
    roh_or, roh_p = stats.fisher_exact([[a, b], [c, d_]])

    return {"status": "ok", "n_satz": n, "n_panel": m, "n_hg": N,
            "beobachtet": beob, "null_mittel": mu, "null_sd": sd,
            "z": float(z), "p": p_emp,
            "OR_gematcht": als_or(beob),
            "mde80_anzahl": mu + 2.8 * sd,
            "OR_mde80": als_or(mu + 2.8 * sd),
            "OR_roh": float(roh_or), "p_roh_fisher": float(roh_p)}


def sym2ens(symbole) -> set:
    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    m = dict(zip(K.symbol.astype(str), K.ensembl))
    return {m[s] for s in symbole if s in m}


def main() -> None:
    log("=" * 78)
    log("Phase M-A  --  human-genetic anchor")
    log("PRAEREG_M_A.md (2026-08-21) + addendum 1 | seed %d | %d draws"
        % (SEED, NZIEHUNGEN))
    log("=" * 78)

    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    HG = raster(K)
    log("\nBackground: %d genes with expression AND length, %d x %d grid"
        % (len(HG), NDEZ, NDEZ))

    P = pd.read_csv(AUS / "panels.csv")
    P = P[P.ensembl.notna()]
    panels = {k: set(v) for k, v in P.groupby("panel").ensembl}
    log("Panels: " + ", ".join("%s %d" % (k, len(v & set(HG.index)))
                               for k, v in sorted(panels.items())))

    programm = set(MODUL.ensembl)
    laesion = set(pd.read_csv(
        ERGEBNISSE / "M_patienten" / "laesionssatz_173.csv").ensembl)
    saetze = {"Programm": programm, "Laesionsantwort": laesion}

    # ---------------------------------------------------------- gate A (a)
    log("\n" + "-" * 78)
    log("POSITIVE CONTROL (a) -- lineage markers in the dysplasia panel")
    eich = []
    marker = sym2ens(OSTEOGEN + CHONDROGEN)
    for pn in ("NOSO", "NOSO_BREIT", "PA309"):
        r = gematcht(marker, panels[pn], HG)
        r.update(teil="a_linienmarker", panel=pn,
                 satz="OSTEOGEN+CHONDROGEN (_marker.py)")
        eich.append(r)
        if r["status"] != "ok":
            log("   %-11s %s" % (pn, r["status"]))
            continue
        log("   %-11s %2d of %2d marker genes in the panel (null %.2f) | "
            "OR %.2f | z %+.2f | p %.4g"
            % (pn, r["beobachtet"], r["n_satz"], r["null_mittel"],
               r["OR_gematcht"], r["z"], r["p"]))
    a_ok = any(e["status"] == "ok" and e["panel"] == "NOSO"
               and e["p"] < 0.05 and e["z"] > 0 for e in eich)
    log("   -> part (a) %s" % ("PASSED" if a_ok else "FAILED"))

    # ---------------------------------------------------------- gate A (b)
    log("\nPOSITIVE CONTROL (b) -- the anchor: distal vs biosynthetic secretion")
    G = pd.read_csv(AUS / "go_saetze.csv")
    S_DISTAL = set(G.ensembl[G.satz == "S_DISTAL"]) & set(HG.index)
    S_BIOSYN = set(G.ensembl[G.satz == "S_BIOSYN"]) & set(HG.index)
    log("   S_DISTAL %d | S_BIOSYN %d genes in the background"
        % (len(S_DISTAL), len(S_BIOSYN)))
    b_ok = False
    for pn in ("NOSO", "NOSO_BREIT", "PA309"):
        pa = panels[pn] & set(HG.index)
        a1, a0 = len(S_DISTAL & pa), len(S_DISTAL - pa)
        b1, b0 = len(S_BIOSYN & pa), len(S_BIOSYN - pa)
        odds, pw = stats.fisher_exact([[a1, a0], [b1, b0]])
        log("   %-11s distal %3d/%4d vs biosyn %3d/%4d | OR %5.2f | p %.3g"
            % (pn, a1, a1 + a0, b1, b1 + b0, odds, pw))
        eich.append(dict(teil="b_anker", panel=pn, satz="S_DISTAL vs S_BIOSYN",
                         status="ok", beobachtet=a1, n_satz=a1 + a0,
                         n_panel=b1 + b0, OR_roh=float(odds),
                         p_roh_fisher=float(pw)))
        if pn == "NOSO":
            b_ok = bool(odds > 2 and pw < 1e-3)
    log("   -> part (b) %s" % ("PASSED" if b_ok else "FAILED"))

    tor_a = a_ok and b_ok
    log("\nGATE A: %s" % ("PASSED -- the numbers of this phase carry"
                          if tor_a else
                          "FAILED -- no number of this phase is reported"))
    pd.DataFrame(eich).to_csv(AUS / "eichung_A.csv", index=False)

    # ------------------------------------------------------- main computation
    log("\n" + "-" * 78)
    log("MAIN COMPUTATION -- enrichment against the matched null")
    log("Bonferroni over %d comparisons: threshold p < %.4f"
        % (BONFERRONI, 0.05 / BONFERRONI))
    log("\n%-16s %-12s %5s %5s %6s %7s %6s %8s %8s %8s"
        % ("set", "panel", "n", "obs", "null", "OR", "z", "p", "OR_MDE80",
           "OR_raw"))
    zeilen = []
    for sn, satz in saetze.items():
        for pn in sorted(panels):
            r = gematcht(satz, panels[pn], HG)
            r.update(satz=sn, panel=pn, tor_a=tor_a)
            zeilen.append(r)
            if r["status"] != "ok":
                log("%-16s %-12s %s" % (sn, pn, r["status"]))
                continue
            log("%-16s %-12s %5d %5d %6.1f %7.2f %+6.2f %8.4g %8.2f %8.2f"
                % (sn, pn, r["n_satz"], r["beobachtet"], r["null_mittel"],
                   r["OR_gematcht"], r["z"], r["p"], r["OR_mde80"],
                   r["OR_roh"]))
    R = pd.DataFrame(zeilen)
    R.to_csv(AUS / "anker.csv", index=False)

    # detection limit separately, as the plan demands
    W = R[R.status == "ok"][["satz", "panel", "n_satz", "n_panel",
                             "null_mittel", "null_sd", "mde80_anzahl",
                             "OR_mde80"]]
    W.to_csv(AUS / "anker_power.csv", index=False)
    log("\nDetection limit (OR at 80 %% power, matched null):")
    log("   program:         %s" % ", ".join(
        "%s %.2f" % (r.panel, r.OR_mde80)
        for _, r in W[W.satz == "Programm"].iterrows()))
    log("   old version for comparison: roughly OR 1.6 at 357 and 50 panel genes.")

    sig = R[(R.status == "ok") & (R.p < 0.05 / BONFERRONI)]
    log("\nAbove the Bonferroni threshold: %d of %d comparisons." % (len(sig),
                                                                     len(R)))
    for _, r in sig.iterrows():
        log("   %s x %s: OR %.2f, z %+.2f, p %.3g"
            % (r.satz, r.panel, r.OR_gematcht, r.z, r.p))
    log("=" * 78)
    (AUS / "52_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()

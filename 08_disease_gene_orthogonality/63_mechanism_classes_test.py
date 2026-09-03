# -*- coding: utf-8 -*-
"""
63_mechanism_classes_test.py -- work package WS2: tests mechanism classes of the
disease genes against the fixed 173/147-gene differentiation program.

Guard compliance: this is NOT a new convergence-axis search (that is spent,
per `Vorwissen/gefallene-hypothesen-guards.md`). The program (`MODUL`,
173 genes resp. 147 in the background grid) and its definition are NOT
changed. What is classified is exclusively the DISEASE-GENE side (the
panels) by external, versioned GO sets from
`08_disease_gene_orthogonality/61_mechanism_classes_go_sets.R`. The draw mechanism (expression- and
length-matched null) is imported verbatim from
`08_disease_gene_orthogonality/11_human_genetics_anchor.py`, not reimplemented.

Output:
  derived_data/followup/ws2_klassentest.csv        -- class x panel, fixed computation
  derived_data/followup/ws2_dwt_test.csv           -- class x |dWT| / dWT direction
  derived_data/followup/ws2_positivkontrolle.csv   -- gate A (a)+(b), reproduced
  derived_data/followup/ws2_power.csv              -- detection limit per class x panel
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

WURZEL = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WURZEL / "00_shared"))

# -- reference implementation, not rewritten -------------------------------
from _module import ERGEBNISSE, MODUL, NZIEHUNGEN, SEED  # noqa: E402
import importlib.util

spec = importlib.util.spec_from_file_location(
    "anker52", WURZEL / "08_disease_gene_orthogonality" / "11_human_genetics_anchor.py")
anker52 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(anker52)  # provides raster(), gematcht()

AUS = WURZEL / "derived_data" / "followup"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


PANELE = ["PA309", "NOSO", "NOSO_BREIT", "KLEIN", "KLEIN_BREIT", "PA1471"]


def main() -> None:
    log("=" * 78)
    log("WS2 -- mechanism classes of the disease genes x differentiation program")
    log("Draw mechanism verbatim from 08_disease_gene_orthogonality/11_human_genetics_anchor.py")
    log("Seed %d | %d draws" % (SEED, NZIEHUNGEN))
    log("=" * 78)

    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    HG = anker52.raster(K)
    log("Background: %d genes (expression AND length known)" % len(HG))

    P = pd.read_csv(ERGEBNISSE / "M_humangenetik" / "panels.csv")
    P = P[P.ensembl.notna()]
    panels = {k: set(v) for k, v in P.groupby("panel").ensembl}

    KL = pd.read_csv(AUS / "ws2_mechanismusklassen_go.csv")
    klassen = {k: set(v) for k, v in KL.groupby("klasse").ensembl}
    log("Mechanism classes: " + ", ".join(
        "%s %d" % (k, len(v & set(HG.index))) for k, v in sorted(klassen.items())))

    programm = set(MODUL.ensembl)  # 173 program genes (147 in the background grid)

    # ---------------------------------------------------- positive control
    log("\n" + "-" * 78)
    log("POSITIVE CONTROL (reproduced from 08_disease_gene_orthogonality/11_human_genetics_anchor.py)")
    from _marker import CHONDROGEN, OSTEOGEN  # noqa: E402

    def sym2ens(symbole):
        m = dict(zip(K.symbol.astype(str), K.ensembl))
        return {m[s] for s in symbole if s in m}

    eich = []
    marker = sym2ens(OSTEOGEN + CHONDROGEN)
    for pn in ("NOSO", "NOSO_BREIT", "PA309"):
        r = anker52.gematcht(marker, panels[pn], HG)
        r.update(teil="a_linienmarker", panel=pn)
        eich.append(r)
        if r["status"] == "ok":
            log("   (a) %-11s OR %.2f | z %+.2f | p %.4g"
                % (pn, r["OR_gematcht"], r["z"], r["p"]))
    a_ok = any(e.get("status") == "ok" and e["panel"] == "NOSO"
               and e["p"] < 0.05 and e["z"] > 0 for e in eich)

    G = pd.read_csv(ERGEBNISSE / "M_humangenetik" / "go_saetze.csv")
    S_DISTAL = set(G.ensembl[G.satz == "S_DISTAL"]) & set(HG.index)
    S_BIOSYN = set(G.ensembl[G.satz == "S_BIOSYN"]) & set(HG.index)
    b_ok = False
    for pn in ("NOSO", "NOSO_BREIT", "PA309"):
        pa = panels[pn] & set(HG.index)
        a1, a0 = len(S_DISTAL & pa), len(S_DISTAL - pa)
        b1, b0 = len(S_BIOSYN & pa), len(S_BIOSYN - pa)
        odds, pw = stats.fisher_exact([[a1, a0], [b1, b0]])
        log("   (b) %-11s distal %d/%d vs biosyn %d/%d | OR %.2f | p %.3g"
            % (pn, a1, a1 + a0, b1, b1 + b0, odds, pw))
        eich.append(dict(teil="b_anker", panel=pn, OR_roh=float(odds), p_roh_fisher=float(pw)))
        if pn == "NOSO":
            b_ok = bool(odds > 2 and pw < 1e-3)
    tor_a = a_ok and b_ok
    log("   -> positive control a %s | b %s | GATE A overall %s"
        % ("PASSED" if a_ok else "FAILED",
           "PASSED" if b_ok else "FAILED",
           "PASSED" if tor_a else "FAILED"))
    pd.DataFrame(eich).to_csv(AUS / "ws2_positivkontrolle.csv", index=False)

    if not tor_a:
        log("\nABORT: positive control failed -- no number of this "
            "phase is reported.")
        (AUS / "ws2_log.txt").write_text("\n".join(LOG), encoding="utf-8")
        return

    # ---------------------------------------------------- main computation: fixed
    log("\n" + "-" * 78)
    n_tests = len(klassen) * len(PANELE)
    log("MAIN COMPUTATION -- mechanism class x panel against the fixed program")
    log("Number of tests: %d classes x %d panels = %d" % (len(klassen), len(PANELE), n_tests))
    schwelle = 0.05 / n_tests
    log("Bonferroni threshold: p < %.5f" % schwelle)

    zeilen = []
    for kn, satz in sorted(klassen.items()):
        for pn in PANELE:
            r = anker52.gematcht(satz, panels[pn], HG)
            r.update(klasse=kn, panel=pn)
            zeilen.append(r)
            if r["status"] != "ok":
                log("%-16s %-12s %s" % (kn, pn, r["status"]))
                continue
            log("%-16s %-12s n=%3d obs=%3d null=%6.2f OR=%6.2f z=%+6.2f p=%8.4g OR_mde80=%6.2f"
                % (kn, pn, r["n_satz"], r["beobachtet"], r["null_mittel"],
                   r["OR_gematcht"], r["z"], r["p"], r["OR_mde80"]))
    R = pd.DataFrame(zeilen)
    R["p_bonferroni_signifikant"] = R.status.eq("ok") & (R.p < schwelle)
    R["n_tests_bonferroni"] = n_tests
    R.to_csv(AUS / "ws2_klassentest.csv", index=False)

    sig = R[R.p_bonferroni_signifikant]
    log("\nAbove the Bonferroni threshold: %d of %d comparisons (of which %d 'ok')."
        % (len(sig), len(R), (R.status == "ok").sum()))
    for _, r in sig.iterrows():
        log("   %s x %s: OR %.2f, z %+.2f, p %.3g" % (r.klasse, r.panel, r.OR_gematcht, r.z, r.p))

    # FDR (Benjamini-Hochberg) as a second, less conservative view
    ok = R[R.status == "ok"].copy()
    ok = ok.sort_values("p").reset_index(drop=True)
    m = len(ok)
    ok["rang"] = np.arange(1, m + 1)
    ok["fdr_schwelle"] = ok["rang"] / m * 0.05
    ok["fdr_signifikant"] = ok["p"] <= ok["fdr_schwelle"]
    # BH: largest rank with p <= threshold; all below count as well
    if ok["fdr_signifikant"].any():
        kk = ok.index[ok["fdr_signifikant"]].max()
        ok.loc[:kk, "fdr_signifikant"] = True
    ok.to_csv(AUS / "ws2_klassentest_fdr.csv", index=False)
    log("FDR (BH, 5%%) significant: %d of %d 'ok' tests" % (ok.fdr_signifikant.sum(), m))

    # detection limit per class x panel, also for null findings
    W = R[R.status == "ok"][["klasse", "panel", "n_satz", "n_panel",
                              "null_mittel", "null_sd", "mde80_anzahl", "OR_mde80"]]
    W.to_csv(AUS / "ws2_power.csv", index=False)

    # ---------------------------------------------------- |dWT| / direction
    log("\n" + "-" * 78)
    log("CONTINUOUS -- class against |dWT| resp. dWT direction (matched null)")
    dwt = K.set_index("ensembl")["dWT_med"]
    zellen = {}
    for z, idx in HG.groupby("zelle").groups.items():
        zellen[z] = HG.index.get_indexer(idx)
    HG_arr_idx = np.arange(len(HG))
    zellen_pos = {z: HG_arr_idx[HG.zelle.values == z] for z in HG.zelle.unique()}
    dwt_hg = dwt.reindex(HG.index).values

    dwt_zeilen = []
    rng_master = np.random.default_rng(SEED + 1)
    for kn, satz in sorted(klassen.items()):
        s = [g for g in satz if g in HG.index]
        n = len(s)
        if n < 8:
            dwt_zeilen.append(dict(klasse=kn, status="too small", n=n))
            log("%-16s too small (n=%d)" % (kn, n))
            continue
        idx_s = HG.index.get_indexer(s)
        obs_abs = float(np.mean(np.abs(dwt_hg[idx_s])))
        obs_dir = float(np.mean(dwt_hg[idx_s]))
        z_satz = HG.zelle.values[idx_s]

        nzieh = NZIEHUNGEN
        # vectorized per cell: all nzieh draws of one gene at once
        sum_abs = np.zeros(nzieh)
        sum_dir = np.zeros(nzieh)
        for z, cnt in zip(*np.unique(z_satz, return_counts=True)):
            pool = zellen_pos[z]
            draw_idx = rng_master.integers(0, len(pool), size=(nzieh, cnt))
            vals = dwt_hg[pool[draw_idx]]
            sum_abs += np.abs(vals).sum(axis=1)
            sum_dir += vals.sum(axis=1)
        null_abs = sum_abs / n
        null_dir = sum_dir / n
        mu_a, sd_a = null_abs.mean(), null_abs.std(ddof=1)
        mu_d, sd_d = null_dir.mean(), null_dir.std(ddof=1)
        z_abs = (obs_abs - mu_a) / sd_a
        z_dir = (obs_dir - mu_d) / sd_d
        p_abs = (1 + int((null_abs >= obs_abs).sum())) / (1 + nzieh)
        p_dir_ein = (1 + int((null_dir >= obs_dir).sum())) / (1 + nzieh)
        p_dir = float(min(1.0, 2 * min(p_dir_ein, 1 - p_dir_ein + 1 / (1 + nzieh))))
        dwt_zeilen.append(dict(
            klasse=kn, status="ok", n=n,
            mean_abs_dwt=obs_abs, null_abs_mittel=mu_a, null_abs_sd=sd_a,
            z_abs=z_abs, p_abs=p_abs,
            mean_dwt_richtung=obs_dir, null_dir_mittel=mu_d, null_dir_sd=sd_d,
            z_dir=z_dir, p_dir=p_dir))
        log("%-16s n=%3d |dWT| %.3f (null %.3f) z=%+.2f p=%.4g | dWT direction %+.3f (null %+.3f) z=%+.2f p=%.4g"
            % (kn, n, obs_abs, mu_a, z_abs, p_abs, obs_dir, mu_d, z_dir, p_dir))
    DW = pd.DataFrame(dwt_zeilen)
    n_tests_dwt = (DW.status == "ok").sum() * 2  # |dWT| and direction per class
    DW["n_tests_bonferroni_kontinuierlich"] = n_tests_dwt
    schwelle_dwt = 0.05 / n_tests_dwt if n_tests_dwt > 0 else np.nan
    DW["p_abs_bonferroni_signifikant"] = DW.status.eq("ok") & (DW.p_abs < schwelle_dwt)
    DW["p_dir_bonferroni_signifikant"] = DW.status.eq("ok") & (DW.p_dir < schwelle_dwt)
    DW.to_csv(AUS / "ws2_dwt_test.csv", index=False)
    log("\nContinuous: number of tests = %d (classes x {|dWT|,direction}), threshold p < %.5f"
        % (n_tests_dwt, schwelle_dwt))

    log("=" * 78)
    (AUS / "ws2_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()

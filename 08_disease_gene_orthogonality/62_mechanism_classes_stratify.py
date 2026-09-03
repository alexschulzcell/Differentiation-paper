# -*- coding: utf-8 -*-
"""
62_mechanism_classes_stratify.py -- WS2: stratification of the skeletal-dysplasia
disease genes by seven external mechanism classes (GO sets from
`60_mechanism_classes_go_build.R`), tested against the FIXED 173/147-gene program.

Guard compliance: NOTHING is changed about the convergence definition or the
differentiation program. What is tested is whether an external,
independently annotated subset of the disease genes (panel ∩ mechanism
class) is over- or underrepresented in the fixed program -- exactly the
matched draw mechanism from `08_disease_gene_orthogonality/11_human_genetics_anchor.py`, only with
swapped roles: instead of "program vs whole panel", "program vs
panel∩class" is drawn. The code in `gematcht()` is imported UNCHANGED, not
rewritten.

Additionally: contrast on the continuous quantity |dWT_med| resp. dWT_med
(signed) between (panel ∩ class) and the rest of the background, with the
existing function `_module.kontrast` (the same background-drawn null as
everywhere in the project).

Output: derived_data/followup/ws2_tests.csv, ws2_ueberlappung.csv,
ws2_positivkontrolle.csv, ws2_klassifikation_panels.csv
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

WURZEL = pathlib.Path(__file__).resolve().parents[1]
ALT_ANALYSE = WURZEL / "08_disease_gene_orthogonality"
NEU = WURZEL / "Neu"
AUS = NEU / "derived_data"
AUS.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ALT_ANALYSE))
from _marker import CHONDROGEN, OSTEOGEN  # noqa: E402
from _module import ERGEBNISSE, MODUL, kontrast  # noqa: E402
import importlib.util as _ilu  # noqa: E402

# 11_human_genetics_anchor.py loaded as a module (importable without hyphen
# despite the leading digit in the filename), in order to reuse `raster`,
# `gematcht`, `exonlaengen`, `sym2ens` UNCHANGED.
_spec = _ilu.spec_from_file_location("anker52", ALT_ANALYSE / "11_human_genetics_anchor.py")
anker52 = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(anker52)

MHUM = ERGEBNISSE / "M_humangenetik"
PANELE_WS2 = ["PA309", "NOSO", "NOSO_BREIT", "KLEIN", "KLEIN_BREIT", "PA1471"]
BONFERRONI = 84   # 7 classes x 6 panels x 2 statistics (categorical + continuous)


def log(msg=""):
    print(msg)


def main():
    log("=" * 78)
    log("WS2 -- mechanism classes of the skeletal dysplasia genes x fixed program")
    log("=" * 78)

    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    HG = anker52.raster(K)
    log("background (draw mechanism 11_human_genetics_anchor.py): %d genes" % len(HG))

    P = pd.read_csv(MHUM / "panels.csv")
    P = P[P.ensembl.notna()]
    panels = {k: set(v) for k, v in P.groupby("panel").ensembl}

    G = pd.read_csv(NEU / "derived_data" / "ws2_go_klassen.csv")
    klassen = {k: set(v) for k, v in G.groupby("klasse").ensembl}
    log("mechanism classes (GO sets, genome): " +
        ", ".join("%s %d" % (k, len(v)) for k, v in sorted(klassen.items())))

    programm = set(MODUL.ensembl)
    laesion_pfad = ERGEBNISSE / "M_patienten" / "laesionssatz_173.csv"
    laesion = set(pd.read_csv(laesion_pfad).ensembl) if laesion_pfad.exists() else set()

    # ---------------------------------------------------------- classification
    # panel x class -> intersection (disease genes of this panel carrying this
    # mechanism class). Multiple membership explicitly allowed.
    klass_panel = {}
    zeilen_klass = []
    for pn in PANELE_WS2:
        pa = panels.get(pn, set())
        for kl, kset in klassen.items():
            schnitt = pa & kset & set(HG.index)
            klass_panel[(pn, kl)] = schnitt
            zeilen_klass.append({"panel": pn, "klasse": kl,
                                  "n_panel_in_hg": len(pa & set(HG.index)),
                                  "n_klasse_genom": len(kset),
                                  "n_schnitt": len(schnitt)})
    KP = pd.DataFrame(zeilen_klass)
    KP.to_csv(AUS / "ws2_klassifikation_panels.csv", index=False)
    log("\nclassification written: ws2_klassifikation_panels.csv (%d rows)" % len(KP))

    # overlap matrix of the classes themselves (genome-wide, on the background)
    kls = sorted(klassen)
    UEB = pd.DataFrame(index=kls, columns=kls, dtype=float)
    for a in kls:
        for b in kls:
            sa, sb = klassen[a] & set(HG.index), klassen[b] & set(HG.index)
            UEB.loc[a, b] = len(sa & sb)
    UEB.to_csv(AUS / "ws2_ueberlappung.csv")
    log("overlap matrix of the classes (gene count, background) written: ws2_ueberlappung.csv")

    # ---------------------------------------------------------- positive control
    log("\n" + "-" * 78)
    log("POSITIVE CONTROL (reproduced from 11_human_genetics_anchor.py, unchanged)")
    eich = []

    def sym2ens_lokal(symbole):
        m = dict(zip(K.symbol.astype(str), K.ensembl))
        return {m[s] for s in symbole if s in m}

    marker = sym2ens_lokal(OSTEOGEN + CHONDROGEN)
    for pn in ("NOSO", "NOSO_BREIT", "PA309"):
        r = anker52.gematcht(marker, panels[pn], HG)
        r.update(teil="a_linienmarker", panel=pn)
        eich.append(r)
        if r["status"] == "ok":
            log("   (a) %-11s lineage markers: OR %.2f | z %+.2f | p %.4g"
                % (pn, r["OR_gematcht"], r["z"], r["p"]))
    a_ok = any(e.get("status") == "ok" and e["panel"] == "NOSO"
               and e["p"] < 0.05 and e["z"] > 0 for e in eich)

    Ggo = pd.read_csv(MHUM / "go_saetze.csv")
    S_DISTAL = set(Ggo.ensembl[Ggo.satz == "S_DISTAL"]) & set(HG.index)
    S_BIOSYN = set(Ggo.ensembl[Ggo.satz == "S_BIOSYN"]) & set(HG.index)
    from scipy import stats as _st
    b_ok = False
    for pn in ("NOSO", "NOSO_BREIT", "PA309"):
        pa = panels[pn] & set(HG.index)
        a1, a0 = len(S_DISTAL & pa), len(S_DISTAL - pa)
        b1, b0 = len(S_BIOSYN & pa), len(S_BIOSYN - pa)
        odds, pw = _st.fisher_exact([[a1, a0], [b1, b0]])
        eich.append(dict(teil="b_anker", panel=pn, status="ok",
                          OR_roh=float(odds), p_roh_fisher=float(pw)))
        log("   (b) %-11s secretion anchor: OR %.2f | p %.3g" % (pn, odds, pw))
        if pn == "NOSO":
            b_ok = bool(odds > 2 and pw < 1e-3)
    log("   -> positive control (a) %s | (b) %s"
        % ("PASSED" if a_ok else "FAILED",
           "PASSED" if b_ok else "FAILED"))
    tor = a_ok and b_ok
    pd.DataFrame(eich).to_csv(AUS / "ws2_positivkontrolle.csv", index=False)
    if not tor:
        log("\n*** GATE NOT PASSED -- the design does not carry. No number "
            "of this phase is reported as a finding. ***\n")

    # ---------------------------------------------------------- main computation
    log("\n" + "-" * 78)
    log("MAIN COMPUTATION -- 7 classes x 6 panels, categorical (matched) + continuous (dWT)")
    log("Bonferroni over %d tests: threshold p < %.5f" % (BONFERRONI, 0.05 / BONFERRONI))

    dWT = K.set_index("ensembl").dWT_med
    dWT_abs = dWT.abs()

    zeilen = []
    for pn in PANELE_WS2:
        for kl in kls:
            satzP = klass_panel[(pn, kl)]
            # categorical: is "panel ∩ class" enriched in the fixed program?
            r = anker52.gematcht(programm, satzP, HG)
            r.update(panel=pn, klasse=kl, kennzahl="kategorial_programm")
            zeilen.append(r)

            # continuous: is the |dWT| of (panel ∩ class) higher than the rest?
            rest = set(HG.index) - satzP
            rk = kontrast(dWT_abs, satzP, rest)
            rk.update(panel=pn, klasse=kl, kennzahl="stetig_abs_dWT")
            zeilen.append(rk)

            # direction: signed dWT
            rk2 = kontrast(dWT, satzP, rest)
            rk2.update(panel=pn, klasse=kl, kennzahl="stetig_dWT_richtung")
            zeilen.append(rk2)

    R = pd.DataFrame(zeilen)
    R.to_csv(AUS / "ws2_tests.csv", index=False)

    ok = R[R.status == "ok"].copy()
    n_kat = (R.kennzahl == "kategorial_programm").sum()
    n_stet = (R.kennzahl != "kategorial_programm").sum()
    log("\nComputed tests: %d categorical + %d continuous = %d (of which 'ok': %d)"
        % (n_kat, n_stet, len(R), len(ok)))

    schwelle = 0.05 / BONFERRONI
    sig = ok[ok.p < schwelle]
    log("\nAbove the Bonferroni threshold p < %.5f: %d of %d computed tests."
        % (schwelle, len(sig), len(ok)))
    for _, r in sig.iterrows():
        if r.kennzahl == "kategorial_programm":
            log("   [categorical] %s x %s: OR %.2f (MDE80 OR %.2f), z %+.2f, p %.3g, n_schnitt=%d"
                % (r.klasse, r.panel, r.OR_gematcht, r.get("OR_mde80", np.nan),
                   r.z, r.p, r.n_panel))
        else:
            log("   [%s] %s x %s: contrast %.3f (MDE80 %.3f), z %+.2f, p %.3g, n=%d"
                % (r.kennzahl, r.klasse, r.panel, r.kontrast, r.mde80, r.z, r.p, r.n_a))

    if len(sig) == 0:
        log("\n-> NO test above the Bonferroni threshold. No mechanism class "
            "shows a robust positive finding against the fixed program.")

    log("=" * 78)
    with open(AUS / "ws2_log.txt", "w", encoding="utf-8") as f:
        f.write("see stdout of this run")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
12_calibration_gene_space.py -- which gene space belongs under the calibration?

Occasion: `03_lineage_calibration/10_calibration_18_datasets.py` did not
reproduce its own stored output. The difference is NOT the annotation -- both
maps come from GENCODE v46 -- but a filter:

    derived_data/R_intern/R_interne_genkarte.csv keeps only genes with
    dWT_n >= 15 of 18 data sets (MIN_N in
    01_expression_landscape/11_internal_gene_map.py).

-----------------------------------------------------------------------------
THE QUESTION, PUT BIOLOGICALLY
-----------------------------------------------------------------------------
The MIN_N filter is right for its own purpose: a gene measurable in only three
data sets cannot be called "convergent across 18 data sets". The filter was
built to define the MODULE.

The calibration asks something else: **did THIS culture reach its lineage?**
That is a question about a single data set. And the filter does not hit it
neutrally -- it preferentially removes the TERMINAL differentiation markers:

    OSTEOGENIC   removed: SP7, BGLAP, IBSP, MEPE, DMP1, SOST
    CHONDROGENIC removed: SOX5, ACAN, COL2A1, COL9A1, COMP, HAPLN1, PRG4
    ADIPOGENIC   removed: FABP4, ADIPOQ, LEP, LPL, PLIN1, PLIN4, CIDEC, ...
    UNDIFFERENTIATED  nothing removed (10 of 10)

That is neither an accident nor a technicality. Terminal markers are off in
the undifferentiated arm and become measurable only where differentiation
actually runs; a filter of "measurable in at least 15 of 18 data sets"
therefore removes them systematically. **The positive control loses exactly
the genes whose presence would demonstrate arrival at the lineage** -- while
the undifferentiated set, which measures the departure from the starting
state, survives complete. The filter therefore acts asymmetrically against one
half of the calibration.

And the data sets do measure these genes: data set 15 measures 18 of 18
osteogenic markers where the internal map passes 12; data sets 3, 4, 8 and 9
measure 12 of 12 chondrogenic markers where the internal map passes 5.

-----------------------------------------------------------------------------
THE DECISION RULE -- before the run, and disclosed
-----------------------------------------------------------------------------
**A disclosure that belongs with it:** both outcomes are already known (3 of
18 under the internal map, 2 of 18 under the full one). The rule is therefore
not presented as blind. It is chosen on grounds that can be stated
independently of the outcome, and it goes against the interest of the paper --
the head count gets worse, not better. That is precisely why it can be
reported this way.

Two kinds of error have to be weighed against each other:

  (E1) The background contains genes whose dWT is essentially noise. The
       spread of the null then grows, every z shrinks, and the calibration
       fails for reasons that have nothing to do with biology.
  (E2) The marker set and the background exclude genes that this data set does
       measure and that are the textbook markers of its axis. The positive
       control is then deprived of its own evidence.

**The rule:** decompose which of the two interventions -- marker set or
background -- moves the z values. Four variants, the same statistic
(`00_shared/_module.py::kontrast`), the same seed:

    A  markers filtered   background filtered   = the stored version
    B  markers full       background filtered
    C  markers filtered   background full
    D  markers full       background full       = today's script run

  * If the MARKER SET carries the movement (B differs from A, C does not),
    E2 applies and the marker set must not be filtered.
  * If the BACKGROUND carries it (C differs from A, B does not), it has to be
    checked whether the additional genes are noise (E1). If they are not --
    comparable baseline level, comparable dWT spread -- there is no reason to
    exclude them.
  * The exchangeability condition of the permutation null requires in any case
    that markers and background come from the SAME population. The mixed
    variants (B, C) are therefore diagnosis, not candidates. Only A and D are
    candidates.

The head count that results is reported without looking at it beforehand.

Outputs  results/eichung_genraum.csv, results/eichung_genraum_log.txt
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pandas as pd

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
sys.path.insert(0, str(WURZEL / "00_shared"))
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import kontrast  # noqa: E402

SITZUNGEN = pathlib.Path(os.environ.get(
    "SCHERENPAPER_SITZUNGEN",
    str(WURZEL.parent / "backups" / "_backup_2026-08-19_vor_paperaufbau")))
GENE20D = SITZUNGEN / "20_Exploration" / "derived_data"
RES = WURZEL / "results"
RES.mkdir(parents=True, exist_ok=True)

MARKERSAETZE = {"OSTEOGEN": OSTEOGEN, "ADIPOGEN": ADIPOGEN,
                "MYOGEN": MYOGEN, "CHONDROGEN": CHONDROGEN, "NAIV": NAIV}
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def main() -> None:
    log("=" * 78)
    log("Which gene space belongs under the calibration? Rule in the script "
        "header.")
    log("=" * 78)

    es = pd.read_csv(WURZEL / "derived_data" / "R_intern" /
                     "ensembl_symbol_gencode46.csv")
    karte = dict(zip(es.ensembl, es.symbol))
    GK = pd.read_csv(WURZEL / "derived_data" / "R_intern" /
                     "R_interne_genkarte.csv")
    intern_ens = set(GK.ensembl)
    intern_sym = set(GK.symbol.dropna())

    KO = pd.read_csv(WURZEL / "derived_data" / "manuscript" / "f1_kohorte.csv")
    arm = dict(zip(KO.punkt, KO.arm))
    name = dict(zip(KO.punkt, KO.datensatz))

    zeilen = []
    for f in sorted(GENE20D.glob("20d_gene_*.csv")):
        G = pd.read_csv(f)
        pt = int(G.punkt.iloc[0])
        G["ens"] = [str(x).split(".")[0] for x in G.gen]
        G["symbol"] = [karte.get(e) for e in G.ens]
        G = G[G.symbol.notna() & G.dWT.notna()]

        voll = G.groupby("symbol").dWT.median()
        gef = G[G.ens.isin(intern_ens)].groupby("symbol").dWT.median()

        a = arm.get(pt, "")
        eigen = "OSTEOGEN" if a == "osteogen" else "CHONDROGEN"
        m_voll = {k: list(v) for k, v in MARKERSAETZE.items()}
        m_gef = {k: [g for g in v if g in intern_sym]
                 for k, v in MARKERSAETZE.items()}

        for var, msets, hg in (("A_marker_gef_hg_gef", m_gef, gef),
                               ("B_marker_voll_hg_gef", m_voll, gef),
                               ("C_marker_gef_hg_voll", m_gef, voll),
                               ("D_marker_voll_hg_voll", m_voll, voll)):
            sa = msets[eigen]
            sb = [g for k, v in msets.items() if k != eigen for g in v]
            # In the mixed variants the markers must be present in the
            # carrying space, otherwise one tests something other than intended.
            traeger = hg
            if var == "B_marker_voll_hg_gef":
                fehlend = [g for g in set(sa) | set(sb)
                           if g in voll.index and g not in hg.index]
                traeger = pd.concat([hg, voll.reindex(fehlend)])
            r = kontrast(traeger, sa, sb)
            zeilen.append({
                "punkt": pt, "datensatz": name.get(pt, ""), "arm": a,
                "variante": var, "n_hintergrund": len(traeger),
                "n_marker_a": r.get("n_a", np.nan),
                "n_marker_b": r.get("n_b", np.nan),
                "kontrast": r.get("kontrast", np.nan),
                "null_sd": r.get("null_sd", np.nan),
                "z": r.get("z", np.nan), "p": r.get("p", np.nan),
                "status": r.get("status", ""),
                "bestanden": bool(r.get("status") == "ok"
                                  and r.get("p", 1) < 0.05
                                  and r.get("kontrast", 0) > 0)})

    T = pd.DataFrame(zeilen)
    T.to_csv(RES / "eichung_genraum.csv", index=False)

    piv = T.pivot(index="punkt", columns="variante", values="z")
    piv["datensatz"] = [name.get(p, "") for p in piv.index]
    log("")
    log("z per point and variant")
    log("%2s %-28s %8s %8s %8s %8s" % ("Pt", "Dataset", "A", "B", "C", "D"))
    for p, r in piv.iterrows():
        log("%2d %-28s %8.2f %8.2f %8.2f %8.2f"
            % (p, str(r.datensatz)[:28], r.A_marker_gef_hg_gef,
               r.B_marker_voll_hg_gef, r.C_marker_gef_hg_voll,
               r.D_marker_voll_hg_voll))

    log("")
    log("How far does each intervention move the z values (median |Delta z|)?")
    d_marker = (piv.B_marker_voll_hg_gef - piv.A_marker_gef_hg_gef).abs()
    d_hg = (piv.C_marker_gef_hg_voll - piv.A_marker_gef_hg_gef).abs()
    log("  marker set swapped only (B - A): median %.3f | maximum %.3f"
        % (d_marker.median(), d_marker.max()))
    log("  background swapped only (C - A): median %.3f | maximum %.3f"
        % (d_hg.median(), d_hg.max()))

    log("")
    log("Head count per variant (passing out of 18):")
    for v in ("A_marker_gef_hg_gef", "B_marker_voll_hg_gef",
              "C_marker_gef_hg_voll", "D_marker_voll_hg_voll"):
        s = T[T.variante == v]
        log("  %-24s %d of %d" % (v, int(s.bestanden.sum()), len(s)))

    (RES / "eichung_genraum_log.txt").write_text(chr(10).join(LOG) + chr(10),
                                                 encoding="utf-8")
    log("-> %s" % (RES / "eichung_genraum.csv"))


if __name__ == "__main__":
    main()

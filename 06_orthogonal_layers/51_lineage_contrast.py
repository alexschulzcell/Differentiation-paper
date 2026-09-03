# -*- coding: utf-8 -*-
"""
51_lineage_contrast.py -- the lineage contrast, second attempt.

===========================================================================
THE DECISION RULE. It stands here because it was fixed BEFORE the first
number this script produced, and it is not changed after the run.
===========================================================================

WHY THERE IS A SECOND ATTEMPT
-----------------------------
In GSE332758 the lineage contrast (osteogenic minus adipogenic) was meant to
show that nothing of the module survives when one lineage is subtracted from
the other -- the positive-logic demonstration that the module is the SHARED
component of both lineages. The finding does not carry there, and for a
diagnosable reason: on the difference axis the adipogenic markers move
correctly (z -3.67 in the promoter window) but the osteogenic ones do not
move at all (z -0.47 to +1.66). The osteogenic axis of that data set carries
no lineage signal in chromatin; calibration L fails in all four windows. A
null result on such an axis means "not measurable".

No choice of marker set repairs that. What is needed is a data set in which
BOTH lineages carry their own signal.

THE RULE
--------
(1) Calibration D (differentiation axis, needs an undifferentiated arm):
    lineage markers minus undifferentiated markers, per axis. Only
    GSE151315 has it.
(2) Calibration L (lineage axis): osteogenic minus adipogenic markers on the
    difference axis log2(OB/AC). Testable in both cohorts.
(3) The module test on the difference axis is reported as a FINDING if
    calibration L passes in **at least 2 of the 4 windows** of a cohort
    (z >= 2 AND contrast > MDE80).
(4) If L passes, a module value BELOW its own detection limit is the result:
    the module is the shared component. A module value ABOVE the limit would
    refute the claim and would be reported just the same.
(5) If L fails, the level stays "not measurable" -- as in GSE332758 -- and
    the panel stays in the supplement.

WHERE THE PANEL GOES
--------------------
If L passes in at least 2 windows in at least one of the two cohorts, the
lineage contrast moves into the main figure. Otherwise it stays in the
supplement and the main text keeps its present, weaker wording.

===========================================================================

The metrics come unchanged from `00_shared/_module.py`
(`kontrast`, `konkordanz`) and `00_shared/_marker.py`. No
metric is implemented a second time.

Inputs    derived_data/B_atac2/<GSE>_matrix_<window>.csv
Outputs   results/linienkontrast_eichung.csv
          results/linienkontrast_modultest.csv
          results/linienkontrast_diagnose.csv
          results/linienkontrast_log.txt
Runtime   a few minutes
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
from _marker import ADIPOGEN, NAIV, OSTEOGEN          # noqa: E402
from _module import MODUL, konkordanz, kontrast, wilson  # noqa: E402

EIN = WURZEL / "derived_data" / "B_atac2"
AUS = WURZEL / "results"
AUS.mkdir(parents=True, exist_ok=True)
EPS = 0.05
FENSTER = ["P", "T10", "T50", "GB"]

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def achsen(M: pd.DataFrame, kohorte: str) -> pd.DataFrame:
    """Log2 axes per cohort. AC = adipogenic, OB = osteogenic, hMSC = undifferentiated."""
    sp = list(M.columns)
    ac = [c for c in sp if "_AC_" in c]
    ob = [c for c in sp if "_OB_" in c]
    naiv = [c for c in sp if "hMSC" in c]
    D = pd.DataFrame(index=M.index)
    if naiv:
        n = M[naiv].mean(axis=1)
        D["osteogen"] = np.log2((M[ob].mean(axis=1) + EPS) / (n + EPS))
        D["adipogen"] = np.log2((M[ac].mean(axis=1) + EPS) / (n + EPS))
        D["basis"] = n
    # The lineage contrast needs NO undifferentiated arm; it is measured directly.
    D["differenz"] = np.log2((M[ob].mean(axis=1) + EPS) /
                             (M[ac].mean(axis=1) + EPS))
    if not naiv:
        D["basis"] = M[ac + ob].mean(axis=1)
    log("    %s: naive %d, AC %d, OB %d columns" % (kohorte, len(naiv),
                                                    len(ac), len(ob)))
    return D


def main() -> None:
    log("=" * 78)
    log("51_lineage_contrast.py -- second attempt on an independent cohort")
    log("=" * 78)
    log("Decision rule see header comment; it was fixed before the run.")

    kohorten = sorted({p.name.split("_matrix_")[0]
                       for p in EIN.glob("*_matrix_*.csv")})
    if not kohorten:
        raise SystemExit("No matrices in %s -- run 22_ first." % EIN)

    eich, mod, diag = [], [], []
    sym_ri = dict(zip(MODUL.symbol, MODUL.ri))

    for ko in kohorten:
        log("\n" + "=" * 72)
        log("Cohort %s" % ko)
        log("=" * 72)
        daten = {}
        for f in FENSTER:
            p = EIN / ("%s_matrix_%s.csv" % (ko, f))
            if p.exists():
                daten[f] = achsen(pd.read_csv(p, index_col=0), ko)

        # ------------------------------------------------ CALIBRATION D
        hat_naiv = "osteogen" in next(iter(daten.values())).columns
        if hat_naiv:
            log("\n--- Calibration D: does the level resolve differentiation?")
            for f, D in daten.items():
                for a, marker in (("osteogen", OSTEOGEN),
                                  ("adipogen", ADIPOGEN)):
                    r = kontrast(D[a], marker, NAIV)
                    ok = r.get("status") == "ok" and r["kontrast"] >= r["mde80"] \
                        and r["z"] >= 2
                    log("  %-4s %-9s n %2d/%2d  contrast %+7.3f  MDE80 %6.3f  "
                        "z %+6.2f  p %8.4g  %s"
                        % (f, a, r.get("n_a", 0), r.get("n_b", 0),
                           r.get("kontrast", np.nan), r.get("mde80", np.nan),
                           r.get("z", np.nan), r.get("p", np.nan),
                           "PASSED" if ok else "failed"))
                    eich.append(dict(kohorte=ko, eichung="D", fenster=f,
                                     achse=a, bestanden=bool(ok),
                                     **{k: v for k, v in r.items()
                                        if k != "status"}))
        else:
            log("\n--- Calibration D: not testable, no undifferentiated "
                "arm in %s" % ko)

        # ------------------------------------------------ CALIBRATION L
        log("\n--- Calibration L: does the level keep the two lineages "
            "apart?")
        l_bestanden = 0
        for f, D in daten.items():
            r = kontrast(D["differenz"], OSTEOGEN, ADIPOGEN)
            ok = r.get("status") == "ok" and r["kontrast"] >= r["mde80"] \
                and r["z"] >= 2
            l_bestanden += int(ok)
            log("  %-4s differenz n %2d/%2d  contrast %+7.3f  MDE80 %6.3f  "
                "z %+6.2f  p %8.4g  %s"
                % (f, r.get("n_a", 0), r.get("n_b", 0),
                   r.get("kontrast", np.nan), r.get("mde80", np.nan),
                   r.get("z", np.nan), r.get("p", np.nan),
                   "PASSED" if ok else "failed"))
            eich.append(dict(kohorte=ko, eichung="L", fenster=f,
                             achse="differenz", bestanden=bool(ok),
                             **{k: v for k, v in r.items() if k != "status"}))

        log("\n  Calibration L passed in %d of %d windows -> lineage "
            "contrast %s"
            % (l_bestanden, len(daten),
               "IS MEASURABLE" if l_bestanden >= 2
               else "stays not measurable"))

        # ------------------------ diagnosis: where each marker set sits
        # The same decomposition that convicted GSE332758 -- here as a
        # cross-check, so that a pass of L cannot rest on one side alone.
        log("\n--- Diagnosis: position of the marker sets on the difference "
            "axis")
        for f, D in daten.items():
            hg = D["differenz"].dropna()
            rng = np.random.default_rng(20260823)
            v = hg.values
            zeile = []
            for name, satz in (("OSTEO", OSTEOGEN), ("ADIPO", ADIPOGEN),
                               ("NAIV", NAIV)):
                g = [x for x in satz if x in hg.index]
                if len(g) < 3:
                    continue
                null = np.array([v[rng.choice(len(v), len(g), replace=False)].mean()
                                 for _ in range(4000)])
                z = (hg[g].mean() - null.mean()) / null.std(ddof=1)
                zeile.append("%s n%d %+.3f (z %+5.2f)"
                             % (name, len(g), hg[g].mean(), z))
                diag.append(dict(kohorte=ko, fenster=f, achse="differenz",
                                 satz=name, n=len(g),
                                 mittel=float(hg[g].mean()), z=float(z)))
            log("  %-4s %s" % (f, " | ".join(zeile)))

        # --------- diagnosis on the SINGLE axes (only with an undifferentiated arm)
        # This decomposition answers the question that calibration D answers
        # only with yes or no: did the culture leave the undifferentiated
        # state, and did it reach its lineage? The two, separately.
        if hat_naiv:
            log("\n--- Diagnosis: position of the marker sets on the single "
                "axes")
            for f, D in daten.items():
                for a_ in ("osteogen", "adipogen"):
                    v = D[a_].dropna()
                    rng = np.random.default_rng(20260823)
                    arr = v.values
                    teil = []
                    for name, satz in (("OSTEO", OSTEOGEN), ("ADIPO", ADIPOGEN),
                                       ("NAIV", NAIV)):
                        g = [x for x in satz if x in v.index]
                        if len(g) < 3:
                            continue
                        null = np.array([arr[rng.choice(len(arr), len(g),
                                                        replace=False)].mean()
                                         for _ in range(4000)])
                        z = (v[g].mean() - null.mean()) / null.std(ddof=1)
                        teil.append("%s %+.3f (z %+5.2f)" % (name, v[g].mean(), z))
                        diag.append(dict(kohorte=ko, fenster=f, achse=a_,
                                         satz=name, n=len(g),
                                         mittel=float(v[g].mean()), z=float(z)))
                    log("  %-4s %-9s %s" % (f, a_, " | ".join(teil)))

        # ------------------------------------------------ MODULE TEST
        log("\n--- Module test")
        for f, D in daten.items():
            achsenliste = (["osteogen", "adipogen", "differenz"] if hat_naiv
                           else ["differenz"])
            for a in achsenliste:
                lic = [e for e in eich
                       if e["kohorte"] == ko and e["fenster"] == f
                       and e["achse"] == a]
                geeicht = bool(lic[0]["bestanden"]) if lic else False
                hg = D[a].dropna()
                m = hg[hg.index.isin(sym_ri)]
                erw = pd.Series({s: sym_ri[s] for s in m.index})
                for hname, sch in (("Hintergrund", None),
                                   ("H1 basisgeschichtet", D["basis"])):
                    r = konkordanz(m, erw, hintergrund=hg, schichtung=sch)
                    if r.get("status") != "ok":
                        continue
                    k = int((np.sign(m.values) ==
                             erw.reindex(m.index).values).sum())
                    lo, hi = wilson(k, len(m))
                    ueber = r["konkordanz"] >= r["konkordanz_mde80"]
                    log("  %-4s %-9s %-20s n %3d  C %.3f [%.3f-%.3f]  "
                        "MDE80 %.3f  z %+6.2f  p %8.4g  %s%s"
                        % (f, a, hname, r["n"], r["konkordanz"], lo, hi,
                           r["konkordanz_mde80"], r["konkordanz_z"],
                           r["konkordanz_p"],
                           "ABOVE" if ueber else "below",
                           "" if geeicht else "  [axis not calibrated]"))
                    mod.append(dict(kohorte=ko, fenster=f, achse=a, null=hname,
                                    geeicht=geeicht, ueber_grenze=bool(ueber),
                                    k=k, wilson_lo=lo, wilson_hi=hi,
                                    **{x: y for x, y in r.items()
                                       if x != "status"}))

    E = pd.DataFrame(eich); M = pd.DataFrame(mod); Dg = pd.DataFrame(diag)
    E.to_csv(AUS / "linienkontrast_eichung.csv", index=False)
    M.to_csv(AUS / "linienkontrast_modultest.csv", index=False)
    Dg.to_csv(AUS / "linienkontrast_diagnose.csv", index=False)

    # ------------------------------------------------------------ Verdict
    log("\n" + "=" * 78)
    log("VERDICT per the rule fixed in advance")
    log("=" * 78)
    for ko in kohorten:
        l = E[(E.kohorte == ko) & (E.eichung == "L")]
        n_ok = int(l.bestanden.sum())
        messbar = n_ok >= 2
        log("%-12s calibration L %d/%d passed -> lineage contrast %s"
            % (ko, n_ok, len(l),
               "MEASURABLE" if messbar else "not measurable"))
        if messbar:
            md = M[(M.kohorte == ko) & (M.achse == "differenz") &
                   (M.null == "H1 basisgeschichtet")]
            if len(md):
                log("             module on the difference axis: above its "
                    "own limit in %d of %d windows (z %+.2f .. %+.2f)"
                    % (int(md.ueber_grenze.sum()), len(md),
                       md.konkordanz_z.min(), md.konkordanz_z.max()))
                if md.ueber_grenze.sum() == 0:
                    log("             -> FINDING: the module is the shared "
                        "component of both lineages.")
                else:
                    log("             -> The claim 'shared component' is "
                        "REFUTED and is reported as such.")
    haupt = any(int(E[(E.kohorte == k) & (E.eichung == "L")].bestanden.sum()) >= 2
                for k in kohorten)
    log("\nPlacement: lineage contrast %s"
        % ("moves to the MAIN FIGURE (Fig. 2)" if haupt
           else "stays in the SUPPLEMENT (S6)"))

    log("\nnumpy %s | pandas %s | python %s"
        % (np.__version__, pd.__version__, sys.version.split()[0]))
    (AUS / "linienkontrast_log.txt").write_text("\n".join(LOG) + "\n",
                                                encoding="utf-8")


if __name__ == "__main__":
    main()

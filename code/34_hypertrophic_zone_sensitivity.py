# -*- coding: utf-8 -*-
"""
34_hypertrophic_zone_sensitivity.py -- does the fetal in vivo trend hang on the
single hypertrophic point?

Occasion  The postnatal computation (code/33_) showed that the textbook
          positive control breaks down at the hypertrophic terminal zone. In
          the fetal atlas the hypertrophic zone contributed 14 cells in ONE
          point of 66. Whether the module trend (rho 0.456, z +4.80, limit
          0.274, results/invivo_spendertest.csv) hangs on that point is an
          open question and belongs settled before submission.

Inputs    derived_data/followup/ws4_modulwert_je_probe.csv          (66 points)
          derived_data/followup/ws4_positivkontrolle_je_probe.csv   (132, filtered
                                                                     to the
                                                                     chondrogenic
                                                                     comparison)
          derived_data/followup/ws4_p3_panel_vs_modul.csv           (66 points)
          -- the stored per-sample CONTRAST VALUES themselves; nothing is
          recomputed, only the selection of points changes.

-----------------------------------------------------------------------------
THE DECISION RULE -- written here before the first number was computed, and
not changed afterwards.
-----------------------------------------------------------------------------
VARIANT A (THE RULE): every point with zone == "HyperChon" removed. Expected,
and secured by an assertion: exactly 1 of 66 points. The detection limit of
each variant is determined AFRESH from that variant's own null (project rule 1:
every selection of points has its own detection limit).

  * If the module trend STAYS above its new limit, the finding does NOT hang
    on the terminal point -- legend material for Figure 3C.
  * If it FALLS below, Figure 3C is reported as carried (in part) by the
    terminal point and the text and legend are changed.

THE SAME test runs for the positive control, because it licenses the axis: if
it falls below its new limit, the axis is not calibratable WITHOUT the
hypertrophic zone -- a substantial limitation that has to be reported.

DESCRIPTIVE (decides nothing):
  * VARIANT B, additionally without the prehypertrophic zone -- the mirror of
    the range that the postnatal atlas can support.
  * The contrast value of the hypertrophic point within the profile of all
    zone means.
  * The disease panel against the module in both variants.

NO REIMPLEMENTATION: test() is imported from code/20_in_vivo_donor_test.py,
all parameters at their defaults (zone permutation within the specimen,
20 000 draws, seed 20260823).

Outputs   results/invivo_hz_empfindlichkeit.csv
          results/invivo_hz_empfindlichkeit_log.txt
Runtime   a few minutes
"""
from __future__ import annotations

import importlib.util
import os
import pathlib

import pandas as pd

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
ERG = WURZEL / "derived_data" / "followup"
RES = WURZEL / "results"

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def spendertest():
    """The trend test from code/20_in_vivo_donor_test.py -- imported."""
    p = WURZEL / "code" / "20_in_vivo_donor_test.py"
    spec = importlib.util.spec_from_file_location("_spendertest", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.test


def main() -> None:
    test = spendertest()

    quellen = {
        "Modul (173 Gene)": ERG / "ws4_modulwert_je_probe.csv",
        "Positivkontrolle chondrogen-naiv": ERG /
            "ws4_positivkontrolle_je_probe.csv",
        "PA309 gegen Modul": ERG / "ws4_p3_panel_vs_modul.csv",
    }
    daten = {}
    for name, p in quellen.items():
        d = pd.read_csv(p)
        if "vergleich" in d.columns:
            d = d[d.vergleich.isin(["chondrogen_vs_naiv", "PA309_vs_Modul"])]
        assert len(d) == 66, (name, len(d))
        daten[name] = d

    # ---- variant A (the rule): without the hypertrophic zone ---------------
    zeilen = []
    log("=" * 78)
    log("Does the fetal trend hang on the single HyperChon point? Rule in "
        "the header.")
    log("=" * 78)

    hz_punkte = {}
    for name, d in daten.items():
        n_hz = int((d.zone == "HyperChon").sum())
        hz_punkte[name] = d.loc[d.zone == "HyperChon", ["zone", "probe",
                                                        "kontrast"]]
        log("%s: %d HyperChon points of %d" % (name, n_hz, len(d)))
        assert n_hz == 1, "Expected: exactly 1 HyperChon point (%s)" % name
        a = d[d.zone != "HyperChon"].copy()
        t = test(a, "%s [ohne HyperChon]" % name)
        zeilen.append(dict(t, variante="A_ohne_HyperChon", quelle=name))
        log("  variant A: rho %.4f | null %.4f +- %.4f | z %+.3f "
            "| new limit rho %.4f | above limit %s"
            % (t["rho"], t["null_mittel"], t["null_sd"], t["z"],
               t["mde80_rho"], t["ueber_mde80"]))

    # ---- variant B (descriptive): also without the prehypertrophic zone ----
    for name, d in daten.items():
        b = d[~d.zone.isin(["HyperChon", "PrehyperChon"])].copy()
        t = test(b, "%s [ohne Hyper+Praehyper]" % name)
        zeilen.append(dict(t, variante="B_ohne_Hyper_Praehyper", quelle=name))
        log("  variant B (descriptive): rho %.4f | z %+.3f | limit %.4f "
            "| above limit %s"
            % (t["rho"], t["z"], t["mde80_rho"], t["ueber_mde80"]))

    # ---- reference: the full run -------------------------------------------
    log("")
    log("Reference (full run, results/invivo_spendertest.csv):")
    ref = pd.read_csv(RES / "invivo_spendertest.csv")
    for _, r in ref.iterrows():
        log("  %-32s rho %.4f | limit %.4f | above limit %s"
            % (r.groesse[:32], r.rho, r.mde80_rho, r.ueber_mde80))

    # ---- the hypertrophic point within the profile -------------------------
    log("")
    log("The HyperChon contrast value within the zone profile (mean per "
        "zone):")
    for name, d in daten.items():
        prof = d.groupby("zone").kontrast.mean()
        hz = float(hz_punkte[name].kontrast.iloc[0])
        log("  %-32s HZ point %+.4f | zone means %s"
            % (name[:32], hz,
               ", ".join("%s %+.3f" % (z, v) for z, v in prof.items())))

    # ---- Entscheidung -------------------------------------------------------
    log("")
    log("=" * 78)
    mod = [z for z in zeilen if z["quelle"].startswith("Modul")
           and z["variante"] == "A_ohne_HyperChon"][0]
    pk = [z for z in zeilen if z["quelle"].startswith("Positivkontrolle")
          and z["variante"] == "A_ohne_HyperChon"][0]
    if bool(pk["ueber_mde80"]):
        pk_text = ("the positive control stays above its new limit "
                   "(rho %.4f > %.4f) -- the axis carries without the "
                   "terminal point too" % (pk["rho"], pk["mde80_rho"]))
    else:
        pk_text = ("the positive control falls below its new limit "
                   "(rho %.4f <= %.4f) -- without HyperChon the axis is "
                   "NOT calibratable" % (pk["rho"], pk["mde80_rho"]))
    if bool(mod["ueber_mde80"]):
        urteil = ("NOT HANGING ON THE TERMINAL POINT: without HyperChon too, "
                  "the module trend stays above its own new limit "
                  "(rho %.4f > %.4f, z %+.2f). %s."
                  % (mod["rho"], mod["mde80_rho"], mod["z"], pk_text))
    else:
        urteil = ("HANGING ON THE TERMINAL POINT: without HyperChon the "
                  "module trend falls below its new limit (rho %.4f <= %.4f, "
                  "z %+.2f) -- Fig. 3C and discussion §3 must be revisited. "
                  "%s."
                  % (mod["rho"], mod["mde80_rho"], mod["z"], pk_text))
    log(urteil)
    log("=" * 78)

    pd.DataFrame(zeilen).to_csv(
        RES / "invivo_hz_empfindlichkeit.csv", index=False)
    # A panel file under figures/data as well (same content), so that the
    # required values in 70_check_numbers.py have a script as their source and
    # no orphan file is created.
    pan = WURZEL / "figures" / "data"
    pan.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(zeilen).to_csv(pan / "TS13_in_vivo_hypertrophic_sensitivity.csv",
                                index=False)
    (RES / "invivo_hz_empfindlichkeit_log.txt").write_text(
        chr(10).join(LOG) + chr(10), encoding="utf-8")
    log("")
    log("-> results/invivo_hz_empfindlichkeit.csv | "
        "figures/data/TS13_in_vivo_hypertrophic_sensitivity.csv")


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
21_postnatal_growth_plate_test.py -- the postnatal anchor: module trend and
positive control on the postnatal growth-plate series, donor as the unit,
zone permutation as the null.

Inputs   results/gse288028_pseudobulk.csv.gz (from
         20_postnatal_growth_plate_annotation.py; one row per (zone, donor)
         with at least 5 cells, columns are gene symbols; the same column
         block as there: module and calibration markers plus 4 000 background
         genes, seed 20260822)
         data_raw/_referenz/wachstumsfuge_zonen/saetze.csv
Outputs  results/gse288028_spendertest.csv, results/gse288028_test_log.txt
Runtime  a few minutes (2 x 12 contrasts x 20 000 draws)

-----------------------------------------------------------------------------
THE DECISION RULE -- in the script header BEFORE the first number, and not
changed afterwards.
-----------------------------------------------------------------------------
> **The level counts only if the positive control passes** (a trend above its
> own detection limit). If it passes, the module trend counts as
> **replicated** when rho lies above its **own** detection limit --
> **regardless of whether it is as large as in the first atlas** (there
> rho 0.456, limit 0.274). If the positive control does not pass, the level is
> reported as **"not calibratable"** and changes nothing in the manuscript.
>
> If the module trend lies above its limit but with the **opposite sign**,
> that is a **contradiction** to the first atlas and is reported as one. The
> zone axis runs in the same direction in both atlases (immature to
> hypertrophic).
>
> The zones of the two atlases are not congruent. What is compared is the
> **common section**: the fetal atlas runs MesCond > ChondroProg > Resting >
> Prolif > Prehyper > Hyper, the postnatal growth plate RZ > PZ > (PHZ) > HZ,
> with common steps Resting<->RZ, Prolif<->PZ, Prehyper<->PHZ, Hyper<->HZ.
> MesCond and ChondroProg have no counterpart and are not compared. If a step
> is missing from the points, the ranks of the others are unchanged (Spearman
> is invariant under a strictly monotone rank transformation).

How the contradiction clause is implemented: the statistic is a rank
correlation, and a decline towards hypertrophy shows up as a negative rho. So
that such a finding cannot vanish quietly below the (upward-directed)
detection limit, the same test from 07_in_vivo_growth_plate/12_fetal_donor_trend_test.py also runs on
the mirrored contrasts. If |rho| then exceeds the limit of the mirrored run,
the decline is statistically real and is reported as a contradiction. That is
the same implementation, not a second test.

Disjointness check (BLOCKING): the intersection of the annotation set
(wachstumsfuge_zonen/saetze.csv, symbols resolved through the HGNC aliases)
with the calibration union (chondrogenic and undifferentiated from
00_shared/_marker.py). On a collision the CALIBRATION set is
reduced by the colliding genes, never the annotation set, and never by
measurability (project rule 4). The script stops if the check is violated or
if a reduced set falls below 3 genes.

The unit is the donor. Fewer than three donors: that is stated, and the MDE80
speaks. No second implementation: kontrast() from
00_shared/_module.py, the trend test imported from
07_in_vivo_growth_plate/12_fetal_donor_trend_test.py.
"""
from __future__ import annotations

import importlib.util
import os
import pathlib
import sys

import numpy as np
import pandas as pd

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
sys.path.insert(0, str(WURZEL / "00_shared"))
from _marker import CHONDROGEN, NAIV  # noqa: E402
from _module import MODUL, kontrast  # noqa: E402

RES = WURZEL / "results"
SATZ_DATEI = WURZEL / "derived_data" / "reference_tables" / "growth_plate_zone_markers.csv"

# Sensitivity mode (descriptive, decides nothing): with POSTNATALE_SUFFIX=_alle
# the same test runs on the all-samples matrix from
# 07_in_vivo_growth_plate/20_postnatal_growth_plate_annotation.py (POSTNATALE_PROBEN=alle).
SUFFIX = os.environ.get("POSTNATALE_SUFFIX", "")

ALIAS = {"PTHR1": "PTH1R"}      # HGNC official symbol; literature name Lee 1996
NZIEHUNGEN = 20000              # as code/27_ (default of _module.kontrast)
SEED_KONTRAST = None            # default of _module.kontrast (20260821)
RANG = {"RZ": 1, "PZ": 2, "PHZ": 3, "HZ": 4}

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def spendertest():
    """The trend test from 07_in_vivo_growth_plate/12_fetal_donor_trend_test.py -- imported."""
    p = WURZEL / "07_in_vivo_growth_plate" / "12_fetal_donor_trend_test.py"
    spec = importlib.util.spec_from_file_location("_spendertest", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.test


def main() -> None:
    test = spendertest()
    PB = pd.read_csv(RES / ("gse288028_pseudobulk%s.csv.gz" % SUFFIX))
    gen_spalten = [c for c in PB.columns
                   if c not in ("zone", "spender", "n_zellen")]
    PB["probe"] = PB["spender"]          # test() expects a probe column

    # ---- disjointness check (blocking) --------------------------------------
    saetze = pd.read_csv(SATZ_DATEI)
    annot_gene = sorted({ALIAS.get(g, g) for g in saetze.gen})
    kollision = sorted(set(annot_gene) & set(CHONDROGEN) | set(annot_gene)
                       & set(NAIV))
    kollision = sorted(set(kollision))
    naiv_red = [g for g in NAIV if g not in kollision]
    chondrogen_red = [g for g in CHONDROGEN if g not in kollision]
    log("Disjointness check: annotation %d genes, calibration union %d, "
        "intersection %s" % (len(annot_gene),
                             len(set(CHONDROGEN) | set(NAIV)),
                             kollision or "EMPTY"))
    log("Calibration sets after reduction: CHONDROGEN %d | NAIV %d -> %d"
        % (len(chondrogen_red), len(NAIV), len(naiv_red)))
    if len(naiv_red) < 3 or len(chondrogen_red) < 3:
        raise SystemExit("Reduced calibration sets too small -- abort.")

    modul_up = MODUL.loc[MODUL.ri > 0, "symbol"].tolist()
    modul_dn = MODUL.loc[MODUL.ri < 0, "symbol"].tolist()

    def punkte(a, b) -> pd.DataFrame:
        zeilen = []
        for _, r in PB.iterrows():
            res = kontrast(r[gen_spalten], a, b,
                           nziehungen=NZIEHUNGEN)   # Default-Seed 20260821
            res.update({"zone": r.zone, "spender": r.spender,
                        "probe": r.spender,
                        "n_zellen": int(r.n_zellen)})
            zeilen.append(res)
        return pd.DataFrame(zeilen)

    def pruefe(df: pd.DataFrame, name: str, satz_a, satz_b):
        fehlen_a = [g for g in satz_a if g not in gen_spalten]
        fehlen_b = [g for g in satz_b if g not in gen_spalten]
        log("")
        log("%s: %d/%d set-a genes, %d/%d set-b genes in the block"
            % (name, len(satz_a) - len(fehlen_a), len(satz_a),
               len(satz_b) - len(fehlen_b), len(satz_b)))
        kt = punkte([g for g in satz_a if g not in fehlen_a],
                    [g for g in satz_b if g not in fehlen_b])
        if (kt.status != "ok").any():
            raise SystemExit("Contrast not computable for all points.")
        t = test(kt[["zone", "spender", "probe", "kontrast"]], name,
                 rang=RANG, spezimen_fn=lambda s: s)
        log("  rho %.4f | null %.4f +- %.4f | z %+.3f | limit rho %.4f"
            "| above limit %s"
            % (t["rho"], t["null_mittel"], t["null_sd"], t["z"],
               t["mde80_rho"], t["ueber_mde80"]))
        return t

    t_pk = pruefe(PB, "Positivkontrolle chondrogen-naiv",
                  chondrogen_red, naiv_red)
    t_mo = pruefe(PB, "Modul (173 Gene)", modul_up, modul_dn)

    # ---- contradiction clause: mirrored contrasts ---------------------------
    log("")
    log("Mirrored contrasts (--kontrast), for the contradiction clause:")
    sp = PB.copy()
    sp[gen_spalten] = -sp[gen_spalten]

    def spiegel(t):
        kt_rows = []
        for _, r in sp.iterrows():
            res = kontrast(r[gen_spalten], modul_up, modul_dn,
                           nziehungen=NZIEHUNGEN)   # Default-Seed 20260821
            res.update({"zone": r.zone, "spender": r.spender,
                        "probe": r.spender})
            kt_rows.append(res)
        kt = pd.DataFrame(kt_rows)
        return test(kt[["zone", "spender", "probe", "kontrast"]],
                    "Modul gespiegelt", rang=RANG, spezimen_fn=lambda s: s)

    t_sp = spiegel(t_mo)
    log("  rho %.4f | z %+.3f | limit %.4f | above limit %s"
        % (t_sp["rho"], t_sp["z"], t_sp["mde80_rho"], t_sp["ueber_mde80"]))
    abfall_echt = bool(t_sp["rho"] > t_sp["mde80_rho"])

    n_spender = PB["spender"].nunique()
    log("")
    log("Donors: %d | points: %d (%s)"
        % (n_spender, len(PB),
           ", ".join("%s:%s" % (r.spender, r.zone) for _, r in PB.iterrows())))
    if n_spender < 3:
        log("WARNING: fewer than three donors -- the permutation of zones "
            "within the donor is weak; the MDE80 speaks.")

    # ---- Entscheidungsregel -------------------------------------------------
    log("")
    log("=" * 78)
    pk_ok = bool(t_pk["ueber_mde80"])
    mo_ok = bool(t_mo["ueber_mde80"])
    vorzeichen_richtig = t_mo["rho"] > 0
    if not pk_ok:
        urteil = "NOT CALIBRATABLE"
        begruendung = ("the positive control lies below its own detection "
                       "limit (rho %.4f <= %.4f)" %
                       (t_pk["rho"], t_pk["mde80_rho"]))
    elif mo_ok and abfall_echt and not vorzeichen_richtig:
        urteil = "CONTRADICTION"
        begruendung = ("module trend significantly DECLINING towards "
                       "hypertrophic")
    elif mo_ok and vorzeichen_richtig:
        urteil = "REPLICATED"
        begruendung = ("positive control passes (rho %.4f > %.4f) and the "
                       "module trend lies above its own limit "
                       "(rho %.4f > %.4f, z %+.2f)"
                       % (t_pk["rho"], t_pk["mde80_rho"], t_mo["rho"],
                          t_mo["mde80_rho"], t_mo["z"]))
    else:
        urteil = "NOT REPLICATED"
        begruendung = ("positive control passes, but the module trend stays "
                       "below its own detection limit "
                       "(rho %.4f <= %.4f)" % (t_mo["rho"], t_mo["mde80_rho"]))
    log("VERDICT per the rule in the script header: %s -- %s."
        % (urteil, begruendung))
    log("(For scale: first atlas rho 0.456, limit 0.274.)")
    log("=" * 78)

    aus = pd.DataFrame([
        dict(t_pk, ebene="postnatale_fuge", groesse="Positivkontrolle "
             "chondrogen-naiv (CHONDROGEN %d gegen NAIV_red %d Gene)"
             % (len(chondrogen_red), len(naiv_red))),
        dict(t_mo, ebene="postnatale_fuge", groesse="Modul (173 Gene)"),
        dict(t_sp, ebene="postnatale_fuge", groesse="Modul gespiegelt "
             "(--kontrast, Widerspruchsklausel)",
             abfall_signifikant=abfall_echt),
    ])
    aus.to_csv(RES / ("gse288028_spendertest%s.csv" % SUFFIX), index=False)
    (RES / ("gse288028_test_log%s.txt" % SUFFIX)).write_text(
        chr(10).join(LOG) + chr(10), encoding="utf-8")
    log("")
    log("-> results/gse288028_spendertest%s.csv" % SUFFIX)


if __name__ == "__main__":
    main()

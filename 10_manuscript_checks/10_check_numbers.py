# -*- coding: utf-8 -*-
"""
10_check_numbers.py -- every number of the manuscript against its panel file.

Purpose   The project rule is that every number in every figure is checked
          against its source file before the figure counts as finished. This
          script makes that verifiable rather than assertable: the values
          quoted in the manuscript and in the legends stand here as required
          values and are compared with what figures/data/ actually holds.

          It exits with code 1 as soon as one number does not match, so it is
          usable as a test run before every submission and in continuous
          integration.

Tolerance absolute 0.006 or relative 0.4 %, appropriate to the rounding used
          in the text.

Inputs    figures/data/*.csv
Outputs   results/numbers_check.txt, and the exit code
Runtime   seconds
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
D = WURZEL / "figures" / "data"
RES = WURZEL / "results"
RES.mkdir(parents=True, exist_ok=True)

ZEILEN: list[str] = []
PRUEFUNGEN: list[tuple[str, bool]] = []


def L(n: str) -> pd.DataFrame:
    return pd.read_csv(D / (n + ".csv"))


def one(x) -> float:
    v = np.asarray(x).ravel()
    assert v.size == 1, ("ambiguous", v)
    return float(v[0])


def pruefe(name: str, ist, soll: float) -> None:
    ist = float(ist)
    ok = abs(ist - soll) < max(0.006, abs(soll) * 0.004)
    PRUEFUNGEN.append((name, ok))
    s = "%s %-44s got %10.4f  want %s" % ("OK  " if ok else "FAIL", name,
                                          ist, soll)
    print(s)
    ZEILEN.append(s)


def main() -> None:
    # ------------------------------------------------------------- FIGURE 1
    e = L("F1D_calibration_per_dataset")
    pruefe("calibration 18: passed", e.passed.sum(), 2)
    pruefe("calibration 18: not calibratable", (~e.calibratable).sum(), 1)
    # No data set reaches its own MDE80 on the calibration (that would be
    # about z 2.8); the two that pass do so on the preregistered rule z >= 2.
    pruefe("calibration 18: above own MDE80", (e.contrast > e.detection_limit).sum(), 0)
    z = L("F1E_calibration_per_cell")
    pruefe("calibration cells: passed", z.passed.sum(), 7)
    sd = L("F1C_screen_diagnoses")
    pruefe("screening: candidates",
      one(sd[sd.step == "with a diagnosis axis"].n), 50)
    pruefe("screening: without undifferentiated arm",
      one(sd[sd.step == "of these: no undifferentiated arm (A1)"].n), 46)

    # ------------------------------------------------------------- FIGURE 2
    b = L("F2B_module_per_dataset")
    pruefe("module above own limit", b.above_detection_limit.sum(), 18)
    pruefe("module z min", b.concordance_z.min(), 5.25)
    pruefe("module z max", b.concordance_z.max(), 13.10)
    a = L("F2A_cross_arm_concordance")
    pruefe("cross-concordance rho", a.rho_module_genes[0], 0.622)
    pruefe("cross-concordance z", a.z[0], 7.03)
    c = L("F2C_pooled_by_calibration")
    pruefe("pooled n failed",
      one(c[c.label == "failed"].n_datasets), 16)
    pruefe("pooled n passed",
      one(c[c.label == "passed"].n_datasets), 2)
    pruefe("pooled failed z",
      one(c[c.label == "failed"].concordance_z), 13.13)
    pruefe("pooled passed z",
      one(c[c.label == "passed"].concordance_z), 12.79)
    g = L("F2D_go_composition")
    for kat, soll in [("matrix remodelling", 7.92), ("matrix components", 4.03),
                      ("cell-cycle exit", 2.11), ("secretory machinery", 0.32)]:
        pruefe("GO OR " + kat, one(g[g.category == kat].odds_ratio), soll)
    pruefe("GO n_module (not 173!)", one(g.n_module_genes.unique()), 147)
    # Gene sets v2: broad, independently curated sets alongside the narrow
    # ones. The decision rule is in 06_orthogonal_layers/61_gene_set_enrichment.R.
    v2 = L("F2D_gene_sets_v2")
    br = v2[v2.variant == "broad"]
    for kat, soll in [("matrix remodelling", 3.62), ("matrix components", 3.38),
                      ("secretory machinery", 0.32), ("cell-cycle exit", 1.80)]:
        pruefe("v2 broad OR " + kat, one(br[br.category == kat].odds_ratio), soll)
    pruefe("v2 broad CI low matrix remodelling",
      one(br[br.category == "matrix remodelling"].ci_low), 1.137)
    pruefe("v2 broad CI high matrix remodelling",
      one(br[br.category == "matrix remodelling"].ci_high), 8.913)
    pruefe("v2 broad CI low matrix components",
      one(br[br.category == "matrix components"].ci_low), 1.201)
    pruefe("v2 broad CI high matrix components",
      one(br[br.category == "matrix components"].ci_high), 7.733)
    # Cell-cycle exit: the confidence interval MUST include 1 -- that is
    # exactly why the number left the running text (case b). If this check
    # fails, the text has to be revisited.
    pruefe("v2 cell-cycle CI includes 1",
      float(one(br[br.category == "cell-cycle exit"].ci_low) < 1.0), 1.0)
    pruefe("v2 n_module equal in both sets",
      float(len(v2.n_module_genes.unique())), 1.0)
    pruefe("v2 n_module", one(v2.n_module_genes.unique()), 147)
    # No category survives Bonferroni -- as stated in the text and the legend.
    pruefe("v2 categories under Bonferroni alpha",
      float((br.p_bonferroni < 0.05 / 6).sum()), 0.0)
    pruefe("GO categories under Bonferroni alpha",
      float((g.p_bonferroni < 0.05 / 6).sum()), 0.0)
    emp = v2[v2.variant == "empfindlichkeit"]
    pruefe("v2 sensitivity OR matrix components",
      one(emp[emp.category == "matrix components"].odds_ratio), 1.55)
    pruefe("v2 sensitivity OR matrix remodelling",
      one(emp[emp.category == "matrix remodelling"].odds_ratio), 1.45)
    pruefe("v2 sensitivity never reaches criterion",
      float(emp.ci_excludes_1.sum()), 0.0)
    at = L("F2E_atac_per_axis")
    h1 = at[at.null_model == "H1 baseline-stratified"]
    for ax, soll in [("adipogenic", 3), ("osteogenic", 4), ("lineage contrast", 0)]:
        pruefe("ATAC H1 above limit " + ax,
          h1[h1.axis == ax].threshold_reached.sum(), soll)
    hg = at[at.null_model == "background"]
    pruefe("ATAC background adipogenic",
      hg[hg.axis == "adipogenic"].threshold_reached.sum(), 2)
    pruefe("ATAC calibration adipogenic",
      h1[h1.axis == "adipogenic"].calibration_passed.sum(), 4)
    pruefe("ATAC calibration osteogenic",
      h1[h1.axis == "osteogenic"].calibration_passed.sum(), 0)
    pruefe("ATAC calibration L (difference)",
      h1[h1.axis == "lineage contrast"].calibration_passed.sum(), 0)

    # --- S9C and S9D: the decomposition in the second cohort (GSE151315).
    # This used to be panel F2F and now stands in the supplement; the
    # required values are unchanged, only their place has moved.
    f = L("S9C_second_cohort_decomposition")
    ost = f[f.axis == "osteogenic"]
    pruefe("S9 undifferentiated markers, osteo axis, max z", ost[ost.quantity == "Undifferentiated"].z.max(),
      -1.46)
    pruefe("S9 undifferentiated markers, osteo axis, min z", ost[ost.quantity == "Undifferentiated"].z.min(),
      -2.80)
    pruefe("S9 osteogenic markers, max z", ost[ost.quantity == "Osteogenic"].z.max(), -1.16)
    pruefe("S9 adipogenic markers on osteo axis, min z",
      ost[ost.quantity == "Adipogenic"].z.min(), 1.70)
    pruefe("S9 module on osteo axis, min z", ost[ost.quantity == "module"].z.min(), 3.10)
    pruefe("S9 module on osteo axis, max z", ost[ost.quantity == "module"].z.max(), 7.00)
    adi = f[f.axis == "adipogenic"]
    pruefe("S9 module on adipo axis, max z", adi[adi.quantity == "module"].z.max(), 7.02)
    s9c = L("S9C_decoupling_second_cohort")
    pruefe("S9 module above own limit", s9c.above_detection_limit.sum(), 8)
    pruefe("S9 calibrations passing", s9c.calibration_passed.sum(), 0)

    # --- F2F: the decomposition across all 18 data sets, preregistered in
    # preregistrations/PRAEREG_F2F.md. The primary run is the unfiltered gene
    # space -- the one the corrected calibration uses.
    zb = L("F2F_decomposition_balance")
    zi = zb[zb.gene_map == "full"]
    pruefe("F2F decomposable (naive markers drop)", one(zi.n_decomposable), 10)
    pruefe("F2F decomposition confirmed", one(zi.n_confirmed), 8)
    pruefe("F2F other case (lineage reached)", one(zi.n_other_case), 2)
    pruefe("F2F decomposition refuted", one(zi.n_refuted), 0)
    pruefe("F2F module above own limit", one(zi.n_module_above_limit), 18)
    zg = zb[zb.gene_map == "filtered"]
    pruefe("F2F filtered confirmed", one(zg.n_confirmed), 7)
    pruefe("F2F filtered refuted", one(zg.n_refuted), 0)
    f18 = L("F2F_decomposition_eighteen")
    # The load-bearing statement: NO data set with a failed calibration
    # reaches its lineage, and ALL with a passed one do.
    eig = f18[f18.quantity == "own-lineage markers"]
    pruefe("F2F lineage reached with failed calibration (tautological)",
      float((eig[eig.calibration == "failed"].z > 2).sum()), 0.0)
    # NOTE: these two checks are TAUTOLOGICAL -- criterion (ii) is the same
    # test statistic at the same threshold as the calibration. They stand
    # here as consistency guards (the two computations must not drift apart),
    # NOT as a finding. The legend says so too.
    pruefe("F2F lineage reached with passed calibration (tautological)",
      float((eig[eig.calibration == "passed"].z > 2).sum()), 2.0)
    # And the composition itself (2 passed, 15 failed, 1 not calibratable) --
    # per DATA SET, not per row: three marker contrasts share the calibration
    # of one data set. These counts also stand in the panel labelling of F2F
    # and in the legend.
    mk18 = f18[f18.quantity != "the module"].drop_duplicates("dataset")
    pruefe("F2F passed", float((mk18.calibration == "passed").sum()), 2.0)
    pruefe("F2F failed", float((mk18.calibration == "failed").sum()), 15.0)
    pruefe("F2F not calibratable", float((mk18.calibration == "not calibratable").sum()),
      1.0)
    s9a = L("S9A_calibration_L_three_cohorts")
    pruefe("calibration L passing, 3 cohorts", s9a.passed.sum(), 0)
    pruefe("calibration L attempts", len(s9a), 12)

    # --- F2G/H/I: held-out validation and robustness of the programme,
    # 05_programme_validation/10_heldout_and_robustness.py.
    g = L("F2G_leave_one_study_out")
    above = g.above_mde80.astype(str).isin(["True", "TRUE"])
    pruefe("F2G held-out above own limit (study)", float(above.sum()), 14.0)
    pruefe("F2G held-out datasets", float(len(g)), 18.0)
    pruefe("F2G held-out median z (study)", float(g.z.median()), 4.66)
    # F2A also carries the matched-null control (same panel as the unmatched null)
    h = L("F2A_matched_nulls_summary")
    hm = h[h.null_type == "matched"]
    pruefe("F2A matched null observed rho", one(hm.rho_observed), 0.635)
    pruefe("F2A matched null z", one(hm.z_sd_units), 6.19)
    pruefe("F2A matched null mean", one(hm.null_mean), 0.136)
    pruefe("F2A matched genes", one(hm.n_genes), 147)
    hu = h[h.null_type == "unmatched"]
    pruefe("F2A unmatched null z", one(hu.z_sd_units), 6.96)
    dr = L("F2H_dropout")
    full_rho = one(dr[dr.scheme == "full"].rho)
    pruefe("F2H full cross-arm rho", full_rho, 0.622)
    pruefe("F2H strongest 20% removed rho",
           one(dr[(dr.scheme == "drop_top_abs_dwt") & (dr.removed_frac == 0.20)].rho), 0.573)
    pruefe("F2H most-expressed 20% removed rho",
           one(dr[(dr.scheme == "drop_top_expression") & (dr.removed_frac == 0.20)].rho), 0.671)
    jk = L("F2H_jackknife")
    pruefe("F2H jackknife min rho", float(jk.rho_without_gene.min()), 0.616)
    pruefe("F2H jackknife max rho", float(jk.rho_without_gene.max()), 0.633)
    # F2I external validation on independent differentiation datasets
    ex = L("F2I_external_validation")
    above = ex.above_mde80.astype(str).isin(["True", "TRUE"])
    pruefe("F2I external datasets", float(len(ex)), 4.0)
    pruefe("F2I external above own limit", float(above.sum()), 3.0)
    zby = ex.set_index("lineage").z
    pruefe("F2I osteogenic z", float(zby.get("osteogenic")), 10.10)
    pruefe("F2I adipogenic z", float(zby.get("adipogenic")), 8.22)
    pruefe("F2I vascular z", float(zby.get("vascular calcification")), 8.99)
    pruefe("F2I chondrogenic z (below limit)", float(zby.get("chondrogenic")), 1.69)

    # ------------------------------------------------------------- FIGURE 3
    t = L("F3C_trend_test_donor")
    m = t[t.quantity == "module (173 genes)"]
    pruefe("in vivo module rho", one(m.rho), 0.456)
    pruefe("in vivo module z", one(m.z), 4.80)
    pruefe("in vivo module MDE80", one(m.detection_limit_rho), 0.274)
    pruefe("in vivo n specimens", one(m.n_specimens), 16)
    # Gene decomposition: the decision rule is in the header of
    # 07_in_vivo_growth_plate/13_fetal_gene_decomposition.py; 10 % is the rule, the rest is
    # descriptive.
    gz = L("TS11_in_vivo_gene_decomposition")
    r10 = gz[gz.fraction_removed == 0.10]
    pruefe("gene decomposition 10 %: rho", one(r10.rho), 0.144)
    pruefe("gene decomposition 10 %: z", one(r10.z), 2.74)
    pruefe("gene decomposition 10 %: own limit", one(r10.detection_limit_rho), 0.125)
    pruefe("gene decomposition 10 %: removed genes", one(r10.n_removed), 18)
    # This check carries the statement "not carried by a few genes". If it
    # fails, the legend of Figure 3C has to be revisited.
    pruefe("gene decomposition 10 %: above own limit",
      float(bool(one(r10.above_detection_limit.astype(str).str.lower() == "true"))), 1.0)
    # And this one carries the counter-statement that stands beside it.
    r20 = gz[gz.fraction_removed == 0.20]
    pruefe("gene decomposition 20 %: below own limit",
      float(bool(one(r20.above_detection_limit.astype(str).str.lower() == "false"))), 1.0)
    # Hypertrophic-zone sensitivity (code/34_): the anchor does NOT hang on
    # the single hypertrophic point. Both guards point in the load-bearing
    # direction -- if one fails, the legend of Figure 3C has to be revisited.
    hz = L("TS13_in_vivo_hypertrophic_sensitivity")
    ha = hz[hz.variant == "A without HyperChon"]
    hm = ha[ha.source == "module (173 genes)"]
    hp = ha[ha.source == "positive control: chondrogenic vs undifferentiated"]
    pruefe("HZ test: module rho without HyperChon", one(hm.rho), 0.430)
    pruefe("HZ test: module z without HyperChon", one(hm.z), 4.65)
    pruefe("HZ test: module new own limit", one(hm.detection_limit_rho), 0.250)
    pruefe("HZ test: module above new limit",
      float(bool(one(hm.above_detection_limit.astype(str).str.lower() == "true"))), 1.0)
    pruefe("HZ test: PC above new limit",
      float(bool(one(hp.above_detection_limit.astype(str).str.lower() == "true"))), 1.0)
    pruefe("HZ test: points without HyperChon", one(hm.n_samples), 65)
    # And the legend of Figure 3C carries the sentence. If this fails, either
    # the sentence is gone or a number has changed -- both are to be checked
    # against results/invivo_hz_empfindlichkeit.csv.
    # The caption checks need the manuscript sources; in a companion checkout
    # without manuscript/ they are skipped, and the data checks above still
    # cover every number.
    cap_datei = WURZEL / "manuscript" / "CAPTIONS_MAIN.md"
    if cap_datei.exists():
        cap = cap_datei.read_text(encoding="utf-8")
        for zahl in ["+0.430", "0.250", "+4.65", "+0.880", "0.615"]:
            pruefe("caption 3C contains HZ-sensitivity value %s" % zahl,
              float(zahl in cap), 1.0)
        # Legend F2F: the data-driven counts stand there as well.
        pruefe("caption F2F counter passed 2",
          float("green = passed, 2;" in cap), 1.0)
        pruefe("caption F2F counter not calibratable 1",
          float("grey = not calibratable, 1" in cap), 1.0)
        # Legend F6: the forest plot points at its table.
        # BMC Genomics cites supplementary material by additional-file name;
        # the tables are Additional file 2, sheets S1-S14.
        for teil in ["Additional file 2: Table S7",
                     "Additional file 2: Table S14"]:
            pruefe("caption F6 refers to %s" % teil, float(teil in cap), 1.0)
        capsup = (WURZEL / "manuscript" / "CAPTIONS_SUPPLEMENT.md").read_text(
            encoding="utf-8")
        pruefe("supplement legends carry table S14",
          float("**Table S14." in capsup), 1.0)

    gr = L("TS11b_in_vivo_gene_ranking")
    pruefe("gene ranking rows", float(len(gr)), 173.0)
    pk = L("F3B_positive_control_per_sample")
    pruefe("PC above own limit", pk.above_detection_limit.sum(), 64)
    pruefe("PC z > 2", pk.fraction_z_above_2.sum(), 65)
    pruefe("PC median z", pk.z.median(), 10.92)
    zz = L("F3D_cells_per_zone")
    pruefe("HyperChon samples",
      one(zz[zz.zone == "HyperChon"].n_samples), 1)
    pruefe("HyperChon cells",
      one(zz[zz.zone == "HyperChon"].n_cells), 8)

    # --- The postnatal growth plate: NOT CALIBRATABLE. The decision rule is
    # in the header of 07_in_vivo_growth_plate/21_postnatal_growth_plate_test.py; the positive
    # control lies below its own detection limit, so by project rule 1 the
    # level carries neither a positive nor a negative finding and the
    # manuscript is unchanged (the "single-source" limitation stands
    # literally). This file deliberately lives under results/ and not under
    # figures/data/: there is no panel for a level that is not calibratable.
    # The checks are guards in both directions: were the positive control to
    # pass on a re-run, the level would have to be decided afresh and the
    # text revisited.
    gse = pd.read_csv(RES / "gse288028_spendertest.csv")
    gp = gse[gse.groesse.str.startswith("Positivkontrolle")]
    pruefe("postnatal plate: PC rho", one(gp.rho), -0.033)
    pruefe("postnatal plate: PC limit", one(gp.mde80_rho), 0.371)
    pruefe("postnatal plate: PC NOT above limit",
      float(not bool(one(gp.ueber_mde80.astype(str).str.lower()
                         == "true"))), 1.0)
    gm = gse[gse.groesse == "Modul (173 Gene)"]
    pruefe("postnatal plate: module NOT above limit",
      float(not bool(one(gm.ueber_mde80.astype(str).str.lower()
                         == "true"))), 1.0)
    gs = gse[gse.groesse.str.startswith("Modul gespiegelt")]
    pruefe("postnatal plate: decline NOT significant",
      float(bool(one(gs.abfall_signifikant.astype(str).str.lower()
                     == "false"))), 1.0)
    pruefe("postnatal plate: donors", one(gm.n_spezimen), 4)
    pruefe("postnatal plate: points", one(gm.n_proben), 12)

    # Sensitivity of the fetal anchor to the hypertrophic zone: see the TS13
    # block above (code/34_, the rule in its header; result: it does NOT hang
    # on the terminal point). Here only the all-samples sensitivity of the
    # postnatal growth plate: code/32_ and 33_ with POSTNATALE_PROBEN=alle
    # and POSTNATALE_SUFFIX=_alle. Descriptive; the verdict of the primary run
    # ("not calibratable") is unchanged. If the positive control fails again
    # on a re-run there is nothing to do; if it PASSES, this check fires and
    # the level would have to be decided afresh.
    ga = pd.read_csv(RES / "gse288028_spendertest_alle.csv")
    gap = ga[ga.groesse.str.startswith("Positivkontrolle")]
    pruefe("postnatal plate all samples: PC NOT above limit",
      float(not bool(one(gap.ueber_mde80.astype(str).str.lower()
                         == "true"))), 1.0)
    pruefe("postnatal plate all samples: points", one(gap.n_proben), 13)

    # ------------------------------------------------------------- FIGURE 4
    k = L("F4B_complementarity")
    d1 = k[(k.side == "disease genes PA309") & (k.compared_to == "S_DISTAL")]
    pruefe("PA309 distal OR", one(d1.odds_ratio_matched), 2.84)
    pruefe("PA309 distal z", one(d1.z), 6.30)
    pruefe("PA309 distal limit", one(d1.odds_ratio_detection_limit), 1.79)
    pr = k[(k.side == "programme") & (k.compared_to == "S_DISTAL")]
    pruefe("programme distal OR", one(pr.odds_ratio_matched), 0.81)
    pruefe("programme distal limit", one(pr.odds_ratio_detection_limit), 2.34)
    v = L("F4C_mode_of_inheritance")
    lo = v[v.variable == "loeuf"]
    pruefe("LOEUF mono", one(lo.median_monoallelic), 0.283)
    pruefe("LOEUF bi", one(lo.median_biallelic), 0.826)
    # NOTE: n_bi is NOT the same number for both axes. LOEUF is measurable in
    # 247 of the biallelic genes, |dWT| in 249.
    pruefe("n mono", one(lo.n_monoallelic), 120)
    pruefe("n bi (LOEUF)", one(lo.n_biallelic), 247)
    pruefe("n bi (|dWT|)", one(v[v.variable == "dWT_abs"].n_biallelic), 249)
    pruefe("dWT P mono/bi", one(v[v.variable == "dWT_abs"].p_value), 0.86)
    dd = L("F4E_dynamics_axis")
    pruefe("programme dWT z", one(dd[dd.gene_set == "programme (173 genes)"].z), 18.10)
    pruefe("PA309 dWT z", one(dd[dd.gene_set == "PanelApp 309"].z), -0.71)
    pruefe("PA309 dWT limit", one(dd[dd.gene_set == "PanelApp 309"].detection_limit_delta), 0.073)
    pruefe("GWAS dWT limit", one(dd[dd.gene_set == "height GWAS"].detection_limit_delta), 0.016)
    pruefe("cell cycle dWT z", one(dd[dd.gene_set == "cell cycle (neg. control)"].z), 3.72)
    # F4G equivalence: for every disease/height set the 95% CI upper bound of the
    # observed effect lies BELOW the MDE80, i.e. a detectable effect is excluded;
    # the cell-cycle positive control clears its limit.
    neg = dd[dd.gene_set.isin(["PanelApp 309", "Nosology (core)", "short stature",
                               "height GWAS"])].copy()
    neg["ci_hi"] = neg.delta_observed + 1.96 * neg.null_sd
    pruefe("F4G disease sets with CI upper < MDE80",
           float((neg.ci_hi < neg.detection_limit_delta).sum()), 4.0)
    cc = dd[dd.gene_set == "cell cycle (neg. control)"].iloc[0]
    pruefe("F4G cell-cycle control above its limit",
           float(cc.delta_observed > cc.detection_limit_delta), 1.0)
    c2 = L("F4D_constraint_publication_matched")
    pruefe("GWAS LOEUF z", one(c2[c2.gene_set == "height GWAS"].z), -4.67)
    pruefe("GWAS LOEUF limit", one(c2[c2.gene_set == "height GWAS"].detection_limit_delta), 0.024)
    pruefe("PA309 LOEUF z", one(c2[c2.gene_set == "PanelApp 309"].z), 1.83)
    kl = L("F4F_mechanism_classes")
    pruefe("class tests created", len(kl), 49)
    pruefe("of which computable", kl.computable.sum(), 40)
    pruefe("above Bonferroni", kl.above_bonferroni.sum(), 0)
    pruefe("above own limit", kl.above_detection_limit.sum(), 4)
    gp = L("F4F_pooled_per_panel")
    pruefe("GWAS pooled OR", one(gp[gp.panel == "height GWAS"].odds_ratio_matched), 1.00)
    pruefe("GWAS pooled limit", one(gp[gp.panel == "height GWAS"].odds_ratio_detection_limit), 1.59)
    pruefe("no panel above own limit", gp.above_detection_limit.sum(), 0)

    # ------------------------------------------------------------- FIGURE 5
    f5 = L("F5A_panel_vs_module_per_zone")
    pruefe("prehyper contrast",
      one(f5[f5.zone == "PrehyperChon"].contrast_median), 0.168)
    pruefe("prehyper z", one(f5[f5.zone == "PrehyperChon"].z_median), 4.98)

    # ------------------------------------------------------------- FIGURE 6
    # F6 is a forest and range plot; the numbers are in
    # figures/data/F6_levels_forest.csv (derived from the same source files as
    # the panels, nothing entered by hand), and the text table is Table S14 in
    # the supplement. The checks are built so that they fire when a row
    # changes direction (estimate moving above or below its own limit).
    eb = L("F6_levels_forest")
    tt = L("TS14_levels_book")
    pruefe("F6: 15 levels, forest and table identical",
      float(len(eb) == len(tt) == 14 and list(eb.level) == list(tt.level)),
      1.0)
    r = eb[eb.level == "Human fetal growth plate, in vivo (E-MTAB-8813)"].iloc[0]
    pruefe("F6 fetal anchor estimate", r.estimate_min, 0.456)
    pruefe("F6 fetal anchor limit", r.detection_limit_min, 0.274)
    r = eb[eb.level == "Donors within one experiment"].iloc[0]
    pruefe("F6 donors S1", r.estimate_min, 0.3489)
    pruefe("F6 donors limit", r.detection_limit_min, 0.3438)
    r = eb[eb.level == "Chromatin H3K27ac (GSE129031)"].iloc[0]
    pruefe("F6 H3K27ac min", r.estimate_min, 0.5896)
    pruefe("F6 H3K27ac max", r.estimate_max, 0.6358)
    pruefe("F6 H3K27ac limit min", r.detection_limit_min, 0.5778)
    r = eb[eb.level == "Adipogenic axis, ATAC (GSE332758)"].iloc[0]
    pruefe("F6 ATAC adipo range", float(r.estimate_max - r.estimate_min), 0.058)
    r = eb[eb.level == "Lineage contrast, ATAC (GSE332758)"].iloc[0]
    # Per WINDOW below its own limit (a range version would be imprecise:
    # the limits vary from window to window).
    atd = at[(at["null_model"] == "H1 baseline-stratified") & (at.axis == "lineage contrast")]
    pruefe("F6 ATAC difference BELOW own limit",
      float((atd.concordance < atd.concordance_detection_limit).all()), 1.0)
    r = eb[eb.level == "Promoter methylome 450K (GSE129266)"].iloc[0]
    pruefe("F6 450K below own limit", float(r.estimate_min < r.detection_limit_min), 1.0)
    r = eb[eb.level == "External cohorts, fixed 173-gene set (Table S12)"].iloc[0]
    pruefe("F6 S12 share", r.estimate_min, 0.6821)
    pruefe("F6 S12 share above limit", float(r.estimate_min > r.detection_limit_min), 1.0)
    r = eb[eb.level == "Human genetics, 7 gene panels"].iloc[0]
    pruefe("F6 GWAS pooled OR", r.estimate_min, 1.00)
    pruefe("F6 GWAS limit", r.detection_limit_max, 1.59)
    r = eb[eb.level == "Disease genes on the dynamics axis"].iloc[0]
    pruefe("F6 dynamics-axis PA309 z", r.estimate_min, -0.71)
    r = eb[eb.level == "Gene constraint, publication-matched"].iloc[0]
    pruefe("F6 constraint PA309 z", r.estimate_min, 1.83)
    r = eb[eb.level == "Undifferentiated state as a predictor (day zero)"].iloc[0]
    pruefe("F6 day-zero z", r.estimate_min, -1.68)
    r = eb[eb.level == "Perturbation datasets, transcriptome (18 series)"].iloc[0]
    pruefe("F6 module z min", r.estimate_min, 5.25)
    pruefe("F6 module z max", r.estimate_max, 13.10)
    # Text table: the head count of the calibration and the derived H3K27ac
    # range are pinned.
    pruefe("TS14 row 1 calibration 2/18",
      float(tt.iloc[0].calibration.startswith("2/18")), 1.0)
    h3row = tt[tt.level.str.startswith("Chromatin H3K27ac")].iloc[0]
    for teil in ["9/9", "+2.97..+3.63"]:
        pruefe("TS14 H3K27ac contains %s" % teil,
          float(teil in (h3row.finding + " " + h3row.value)), 1.0)

    # ---------------------------------------------------------- Supplement
    s3 = L("S3A_external_triangulation")
    f = s3[s3.has_detection_limit]
    pruefe("S3 primary z", f.iloc[0].z, 2.84)
    pruefe("S3 primary MDE80", f.iloc[0].detection_limit, 0.679)
    s3b = L("S3B_triangulation_per_study")
    pruefe("S3 studies above own limit", s3b.above_detection_limit.sum(), 0)
    s4 = L("S4B_patient_concordance")
    pruefe("S4 cohorts above own limit", s4.above_detection_limit.sum(), 0)
    tg = L("TS4_module_genes")
    pruefe("module genes total", len(tg), 173)
    pruefe("of which ri = +1", (tg.direction_ri > 0).sum(), 129)
    pruefe("of which ri = -1", (tg.direction_ri < 0).sum(), 44)
    ts8 = L("TS8_preregistrations")
    pruefe("preregistration documents (34 + PRAEREG_F2F)", len(ts8), 35)

    # ------------------------------------------------------------- Balance
    schlecht = [n for n, ok in PRUEFUNGEN if not ok]
    kopf = "=== %d of %d checks passed ===" % (
        len(PRUEFUNGEN) - len(schlecht), len(PRUEFUNGEN))
    print("\n" + kopf)
    ZEILEN.append("")
    ZEILEN.append(kopf)
    if schlecht:
        print("FAILED:", schlecht)
        ZEILEN.append("FAILED: " + ", ".join(schlecht))
    ZEILEN.append("numpy %s | pandas %s | python %s"
                  % (np.__version__, pd.__version__, sys.version.split()[0]))
    (RES / "numbers_check.txt").write_text("\n".join(ZEILEN) + "\n",
                                            encoding="utf-8")
    if schlecht:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

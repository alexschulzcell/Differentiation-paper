# -*- coding: utf-8 -*-
"""
10_panel_data_main.py -- one CSV per PANEL, not per figure.

Purpose  Writes every data file for the main and supplementary figures into
         figures/data/. Each CSV holds only what that one panel draws -- no
         left-over columns. The naming scheme is F2E_atac_per_axis.csv,
         S8C_study_intensity.csv and so on.

Rule     NO second implementation of any metric. This script recomputes
         nothing; it reads the output files of the reference implementations
         (`00_shared/_module.py`,
         `02_matrix_programme_derivation/12_metric_reference.R`) and
         reshapes them. The single exception is
         `results/invivo_spendertest.csv` from
         `07_in_vivo_growth_plate/12_fetal_donor_trend_test.py`.

Inputs   derived_data/**, derived_data/followup/**,
         derived_data/manuscript/**, results/**
Outputs  figures/data/*.csv  (one CSV per panel)
         results/paneldaten_log.txt
Runtime  a few seconds
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
import _display  # noqa: E402  presentation layer: everything shown is English

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
ERG = WURZEL / "derived_data"
NEU = WURZEL / "derived_data" / "followup"
PDAT = WURZEL / "derived_data" / "manuscript"
RES = WURZEL / "results"
AUS = WURZEL / "figures" / "data"
AUS.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def schreib(df: pd.DataFrame, name: str, was: str) -> None:
    """Write one panel CSV and log it in a single line.

    What leaves here is English: _display.englisch() translates the display
    values and the column headers before the file is written.
    """
    p = AUS / f"{name}.csv"
    df = _display.englisch(df)
    df.to_csv(p, index=False)
    log(f"  {name:38s} {len(df):5d} rows x {df.shape[1]:2d} cols  {was}")


# =============================================================================
# FIGURE 1 -- material, the screens, and the calibration almost nothing passes
# =============================================================================
def f1() -> None:
    log("\nFig. 1 -- material, screens, calibration")

    # --- F1B  screen of the perturbation data sets: 89 -> 18
    # The head counts come from the screening record; drawn as a waterfall.
    b = pd.DataFrame([
        {"stufe": "screened series", "n": 89, "rang": 1},
        {"stufe": "no undifferentiated arm (A1)", "n": -36, "rang": 2},
        {"stufe": "other exclusions", "n": -35, "rang": 3},
        {"stufe": "included datasets", "n": 18, "rang": 4},
    ])
    schreib(b, "F1B_screen_perturbation", "waterfall 89 -> 18")

    # --- F1C  diagnosis screen 1424 -> 50 -> 46 without an undifferentiated arm
    s = pd.read_csv(ERG / "M_diagnosen" / "sichtung.csv")
    # The candidate count is formed from the exclusion codes: 127 hand checks
    # minus the 77 with code A2 ("no diagnosis as the lesion axis") = 50
    # candidates. The boolean column `diagnose` is wider (53) and is NOT used
    # -- otherwise the denominator of the 92 % figure would be wrong.
    n_hand = int(len(s))
    n_a2 = int((s.code == "A2").sum())
    n_kandidat = n_hand - n_a2
    n_a1 = int((s.code == "A1").sum())
    c = pd.DataFrame([
        {"stufe": "GEO series screened", "n": 1424, "ebene": "Entitaet"},
        {"stufe": "hand-checked", "n": n_hand, "ebene": "Entitaet"},
        {"stufe": "with a diagnosis axis", "n": n_kandidat, "ebene": "Entitaet"},
        {"stufe": "of these: no undifferentiated arm (A1)", "n": n_a1, "ebene": "Entitaet"},
        {"stufe": "series screened by design", "n": 22, "ebene": "Anlage"},
        {"stufe": "complete 2x2", "n": 2, "ebene": "Anlage"},
    ])
    schreib(c, "F1C_screen_diagnoses",
            f"{n_hand} hand-checked, {n_kandidat} candidates, {n_a1} without "
            f"an undifferentiated arm = {100 * n_a1 / n_kandidat:.0f} %")

    # --- S7  exclusion codes of the hand screen (for S7 and the table)
    codes = (s[s.urteil == "AUS"].groupby("code")
             .agg(n=("gse", "size"), beispiel=("begruendung", "first"))
             .reset_index().sort_values("n", ascending=False))
    # Presentation layer: the reasons are German in the source file; for the
    # English panel S7 they are translated here (shortened to the 40
    # characters the panel carries; the full wording is in Table S2 and in the
    # legend of S7).
    codes["beispiel"] = codes["beispiel"].replace({
        "keine Diagnose als Laesionsachse": "no diagnosis as the lesion axis",
        "kein naiver Arm": "no undifferentiated arm",
        "humaner OPLL-Arm ohne Differenzierungsachse; die Tag-0/14-Achse "
        "liegt in MC3T3-E1 -- Maus, immortalisiert":
            "human OPLL arm; day 0/14 axis is mouse",
        "Retinaorganoide, Mueller-Glia -- keine skelettale Entitaet":
            "retina organoids; not a skeletal entity",
    })
    schreib(codes, "S7A_exclusion_codes",
            "exclusion codes of the diagnosis screen")

    # --- F1D  the calibration per data set: 2 of 18 (in the corrected gene
    # space, see 03_lineage_calibration/12_calibration_gene_space.py)
    e = pd.read_csv(ERG / "M_kalibrierung" / "eichung_achtzehn.csv")
    d = e[["punkt", "datensatz", "arm", "z", "mde80", "kontrast", "p",
           "status", "bestanden", "n_gene_messbar"]].copy()
    d["eichbar"] = d.status == "ok"
    d = d.sort_values("z", na_position="first")
    schreib(d, "F1D_calibration_per_dataset",
            f"passed {int(d.bestanden.sum())}/{len(d)}, "
            f"not calibratable {int((~d.eichbar).sum())}")

    # --- F1E  the same, resolved by donor: 7 of 14
    z = pd.read_csv(ERG / "M_donoren" / "eichung.csv")
    dz = z[["zelle", "spender", "studie", "achse", "z", "mde80", "kontrast",
            "p", "bestanden", "laesion", "e2", "herkunft"]].copy()
    dz = dz.sort_values("z")
    schreib(dz, "F1E_calibration_per_cell",
            f"passed {int(dz.bestanden.sum())}/{len(dz)}")


# =============================================================================
# FIGURE 2 -- the main finding
# =============================================================================
def f2() -> None:
    log("\nFig. 2 -- line independence and decoupling")

    # --- F2A  cross-arm concordance, osteogenic against chondrogenic
    a = pd.read_csv(NEU / "ws6_p1_arm_kreuzkonkordanz.csv")
    schreib(a, "F2A_cross_arm_concordance", "rho +0.622, z +7.03")

    # --- F2B  module concordance per data set, coloured by calibration status
    m = pd.read_csv(NEU / "ws6_p1p2_modul_je_datensatz.csv")
    e = pd.read_csv(ERG / "M_kalibrierung" / "eichung_achtzehn.csv")[
        ["punkt", "status"]].rename(columns={"status": "eichstatus"})
    b = m.merge(e, on="punkt", how="left")
    b = b[["punkt", "datensatz", "arm", "n", "konkordanz", "konkordanz_mde80",
           "konkordanz_null", "konkordanz_null_sd", "konkordanz_z",
           "konkordanz_p", "bestanden", "eichstatus"]].copy()
    b["ueber_grenze"] = b.konkordanz > b.konkordanz_mde80
    b["eichung"] = np.where(b.eichstatus != "ok", "not calibratable",
                            np.where(b.bestanden, "passed", "failed"))
    b = b.sort_values("konkordanz_z")
    schreib(b, "F2B_module_per_dataset",
            f"above own limit {int(b.ueber_grenze.sum())}/{len(b)}, "
            f"z {b.konkordanz_z.min():+.2f}..{b.konkordanz_z.max():+.2f}")

    # --- F2C  pooled: failed against passed
    c = pd.read_csv(NEU / "ws6_p2_gepoolt_bestanden_vs_durchgefallen.csv")
    c = c[["label", "n_datensaetze", "konkordanz", "konkordanz_null",
           "konkordanz_null_sd", "konkordanz_z", "konkordanz_p",
           "konkordanz_mde80"]]
    schreib(c, "F2C_pooled_by_calibration",
            "failed n%d vs passed n%d"
            % (int(c[c.label == "durchgefallen"].n_datensaetze.iloc[0]),
               int(c[c.label == "bestanden"].n_datensaetze.iloc[0])))

    # --- F2D  composition: external GO sets
    d = pd.read_csv(NEU / "ws6_p3_go_annotation.csv")
    d = d.rename(columns={"kategorie": "kategorie_de"})
    d["category"] = d.kategorie_de.map({
        "1_Matrixbestandteile": "matrix components",
        "2_Matrixremodellierung": "matrix remodelling",
        "3_Sekretionsmaschine": "secretory machinery",
        "4_TGFb_BMP": "TGFb / BMP",
        "5_Hypoxie_Stress": "hypoxia / stress",
        "6_Zellzyklusausstieg": "cell-cycle exit"})
    d = d[["category", "kategorie_de", "k", "n_modul", "n_satz_im_hg",
           "n_hintergrund", "OR", "lo", "hi", "p", "p_bonferroni"]]
    d = d.sort_values("OR", ascending=False)
    schreib(d, "F2D_go_composition",
            "n_modul = 147 of the 173 genes lie in the GO background")

    # --- F2D  broad, independently curated gene sets from
    # `06_orthogonal_layers/61_gene_set_enrichment.R`; the sets were fixed before the run. The narrow
    # and the broad set stand SIDE BY SIDE -- nothing is replaced without both
    # staying visible.
    v2 = pd.read_csv(RES / "gensaetze_v2_anreicherung.csv")
    v2["category"] = v2.kategorie.map({
        "1_Matrixbestandteile": "matrix components",
        "2_Matrixremodellierung": "matrix remodelling",
        "3_Sekretionsmaschine": "secretory machinery",
        "4_TGFb_BMP": "TGFb / BMP",
        "5_Hypoxie_Stress": "hypoxia / stress",
        "6_Zellzyklusausstieg": "cell-cycle exit"})
    v2["set"] = v2.variante.map({"eng": "narrow (GO)",
                                 "breit": "broad (v2)",
                                 "empfindlichkeit": "sensitivity (Reactome)"})
    v2 = v2.rename(columns={"kategorie": "kategorie_de"})
    v2 = v2[["category", "kategorie_de", "variante", "set", "quelle",
             "bestandteile", "k", "n_modul", "n_satz_im_hg", "n_hintergrund",
             "OR", "lo", "hi", "p", "p_bonferroni", "ki_schliesst_1_aus"]]
    reihe = v2[v2.variante == "breit"].sort_values("OR", ascending=False)
    v2["category"] = pd.Categorical(v2.category, list(reihe.category),
                                    ordered=True)
    v2 = v2.sort_values(["category", "variante"])
    schreib(v2, "F2D_gene_sets_v2",
            "narrow | broad | sensitivity, n_modul = "
            f"{int(v2.n_modul.iloc[0])} of the 173")

    # --- F2E  ATAC per axis, with the calibration status. Both nulls are
    # carried; H1 (stratified by baseline accessibility) is the primary one,
    # because it is the harder of the two.
    mt = pd.read_csv(ERG / "B_atac" / "B_atac_modultest_final.csv")
    ei = pd.read_csv(ERG / "B_atac" / "B_atac_eichung_je_achse.csv")
    ei_k = ei[["fenster", "achse", "eichung", "bestanden", "z", "kontrast",
               "mde80"]].rename(columns={
                   "bestanden": "eichung_bestanden", "z": "eichung_z",
                   "kontrast": "eichung_kontrast", "mde80": "eichung_mde80",
                   "eichung": "eichungstyp"})
    e2 = mt.merge(ei_k, on=["fenster", "achse"], how="left")
    e2 = e2[["fenster", "achse", "null", "n", "konkordanz", "konkordanz_mde80",
             "konkordanz_null", "konkordanz_null_sd", "konkordanz_z",
             "konkordanz_p", "schwelle_erreicht", "eichungstyp",
             "eichung_bestanden", "eichung_z", "eichung_kontrast",
             "eichung_mde80", "urteil"]]
    schreib(e2, "F2E_atac_per_axis",
            "GSE332758, 4 windows x 3 axes x 2 nulls")

    # Cross-check: the earlier file must be contained row for row
    alt = pd.read_csv(NEU / "ws8_atac_linienunabhaengigkeit.csv")
    neu = e2[e2.null == "Hintergrund"]
    j = alt.merge(neu, left_on=["fenster", "achse"],
                  right_on=["fenster", "achse"], how="left")
    ok = np.allclose(j.konkordanz_x, j.konkordanz_y) and \
        np.allclose(j.konkordanz_mde80_x, j.konkordanz_mde80_y) and \
        np.allclose(j.konkordanz_z_x, j.konkordanz_z_y)
    log(f"     cross-check against the ws8 orphan file: "
        f"{'sign-identical' if ok else 'MISMATCH'}")
    if not ok:
        raise SystemExit("F2E deviates from ws8_atac_linienunabhaengigkeit.csv")

    # --- F2F  the three-way decomposition across all 18 data sets.
    # Preregistered in preregistrations/PRAEREG_F2F.md; computed by
    # 04_programme_decomposition/10_decomposition_18_datasets.py. The primary run is the FULL gene
    # space -- the same one the corrected calibration uses. The filtered
    # version stays in the source file and in Table S10, not in the main
    # figure.
    z18 = pd.read_csv(RES / "zerlegung_achtzehn.csv")
    z18 = z18[z18.genkarte == "voll"].copy()
    lang = []
    for _, r in z18.iterrows():
        eich = ("not calibratable" if r.eichung_status_abgelegt != "ok"
                else "passed" if r.eichung_bestanden_abgelegt else "failed")
        for schl, gr in (("naiv_z", "naive markers"),
                         ("eigen_z", "own-lineage markers"),
                         ("ander_z", "adipogenic markers"),
                         ("modul_z", "the module")):
            lang.append({"punkt": int(r.punkt), "datensatz": r.datensatz,
                         "arm": r.arm, "groesse": gr, "z": r[schl],
                         "eichung": eich, "urteil": r.urteil,
                         "zerlegbar": bool(r.i_naiv_faellt)})
    f2f = pd.DataFrame(lang)
    schreib(f2f, "F2F_decomposition_eighteen",
            "18 datasets x 4 quantities, full gene space")

    b = pd.read_csv(RES / "zerlegung_achtzehn_bilanz.csv")
    schreib(b, "F2F_decomposition_balance",
            "tally per gene map, preregistered rule")


# =============================================================================
# FIGURE 3 -- in vivo, the human fetal growth plate
# =============================================================================
def f3() -> None:
    log("\nFig. 3 -- in-vivo anchor")
    ZONEN = ["MesCond", "ChondroProg", "RestingChon", "ProlifChon",
             "PrehyperChon", "HyperChon"]

    # --- F3A  the atlas: cell counts per zone and sample
    a = pd.read_csv(NEU / "ws4_proben_je_zone.csv")
    a["spezimen"] = a.probe.str.extract(r"^(Pcw[0-9.]+(?:_s[0-9])?)")
    a["stadium"] = a.probe.str.extract(r"^(Pcw[0-9.]+)")
    a["zonenrang"] = a.zone.map({z: i + 1 for i, z in enumerate(ZONEN)})
    schreib(a, "F3A_atlas_zones",
            f"{a.zone.nunique()} zones, {a.probe.nunique()} samples, "
            f"{a.spezimen.nunique()} specimens, {int(a.n_zellen.sum())} cells")

    # --- F3B  the positive control per sample
    pk = pd.read_csv(NEU / "ws4_positivkontrolle_je_probe.csv")
    pk = pk[pk.vergleich == "chondrogen_vs_naiv"].copy()
    pk["spezimen"] = pk.probe.str.extract(r"^(Pcw[0-9.]+(?:_s[0-9])?)")
    pk["zonenrang"] = pk.zone.map({z: i + 1 for i, z in enumerate(ZONEN)})
    pk["ueber_mde80"] = pk.kontrast > pk.mde80
    pk["z_ueber_2"] = pk.z > 2
    schreib(pk[["zone", "zonenrang", "probe", "spezimen", "kontrast", "mde80",
                "z", "p", "ueber_mde80", "z_ueber_2"]],
            "F3B_positive_control_per_sample",
            f"above own limit {int(pk.ueber_mde80.sum())}/{len(pk)}, "
            f"z>2 in {int(pk.z_ueber_2.sum())}/{len(pk)}, "
            f"median z {pk.z.median():+.2f}")

    # --- F3C  module value along the axis, per specimen
    mw = pd.read_csv(NEU / "ws4_modulwert_je_probe.csv")
    mw["spezimen"] = mw.probe.str.extract(r"^(Pcw[0-9.]+(?:_s[0-9])?)")
    mw["zonenrang"] = mw.zone.map({z: i + 1 for i, z in enumerate(ZONEN)})
    mw = mw.merge(a[["zone", "probe", "n_zellen"]], on=["zone", "probe"],
                  how="left")
    schreib(mw[["zone", "zonenrang", "probe", "spezimen", "n_zellen",
                "kontrast", "mde80", "z", "p"]],
            "F3C_module_per_sample", "one line per specimen")

    # --- F3C  the trend test, donor-stratified (the number for the legend)
    t = pd.read_csv(RES / "invivo_spendertest.csv")
    schreib(t, "F3C_trend_test_donor",
            "zone permutation within the specimen")

    # --- F3D  honesty row: cell count per zone
    d = (a.groupby(["zone", "zonenrang"])
         .agg(n_proben=("probe", "size"), n_zellen=("n_zellen", "sum"),
              n_zellen_min=("n_zellen", "min"))
         .reset_index().sort_values("zonenrang"))
    schreib(d, "F3D_cells_per_zone",
            "HyperChon carries only one evaluable sample")


# =============================================================================
# FIGURE 4 -- where the disease genes are
# =============================================================================
def f4() -> None:
    log("\nFig. 4 -- the other layer")

    # --- F4A  positive controls of the set-up
    a = pd.read_csv(ERG / "M_humangenetik" / "eichung_A.csv")
    lin = a[a.teil == "a_linienmarker"][
        ["panel", "n_panel", "n_satz", "beobachtet", "null_mittel",
         "OR_gematcht", "OR_mde80", "z", "p"]].copy()
    lin["teil"] = "lineage markers in panel"
    ank = a[a.teil == "b_anker"][["panel", "n_panel", "n_satz", "beobachtet",
                                  "OR_roh", "p_roh_fisher"]].copy()
    ank["teil"] = "distal vs biosynthetic secretion"
    schreib(lin, "F4A_positive_control_lineage_markers",
            f"OR {lin.OR_gematcht.min():.1f}-{lin.OR_gematcht.max():.1f}")
    schreib(ank, "F4A_positive_control_anchor",
            f"OR {ank.OR_roh.min():.2f}-{ank.OR_roh.max():.2f}")

    # --- F4B  localisation: programme against disease genes
    b = pd.read_csv(NEU / "ws7_komplementaritaet.csv")
    b = b[["gegen", "seite", "n_satz", "n_panel", "beobachtet", "null_mittel",
           "null_sd", "OR_gematcht", "OR_mde80", "z", "p"]].copy()
    b["ueber_grenze"] = b.OR_gematcht > b.OR_mde80
    b["satz"] = b.seite.map({
        "Programm": "programme (173 genes)",
        "Krankheitsgene_PA309": "PanelApp 309",
        "Krankheitsgene_NOSO": "Nosology (core)",
        "Krankheitsgene_NOSO_BREIT": "Nosology (broad)"})
    b["kompartiment"] = b.gegen.map({"S_DISTAL": "distal secretion",
                                     "S_BIOSYN": "biosynthetic secretion"})
    schreib(b, "F4B_complementarity", "8 Tests, Bonferroni alpha 0.00625")

    # --- F4C  gene dosage by mode of inheritance
    c = pd.read_csv(NEU / "ws1_p5_vererbungsmodus.csv")
    schreib(c, "F4C_mode_of_inheritance", "PanelApp 309, mono vs bi")

    # --- F4D  the constraint contrast, publication-matched
    p6 = pd.read_csv(NEU / "ws1_p6_publikationsmatching.csv")
    c2 = p6[p6.lauf == "P3_loeuf_mit_pub"][
        ["satz", "n_satz", "beobachtet", "null_mittel", "null_sd", "z", "p",
         "mde80_delta", "delta_beobachtet", "matching"]].copy()
    c2["ueber_grenze"] = c2.delta_beobachtet.abs() > c2.mde80_delta
    schreib(c2, "F4D_constraint_publication_matched",
            "LOEUF after matching expression, length and publication count")

    # --- F4D  the axis that does not separate them
    d = p6[p6.lauf == "P2_dynamik_mit_pub"][
        ["satz", "n_satz", "beobachtet", "null_mittel", "null_sd", "z", "p",
         "mde80_delta", "delta_beobachtet", "matching"]].copy()
    d["ueber_grenze"] = d.delta_beobachtet.abs() > d.mde80_delta
    schreib(d, "F4E_dynamics_axis",
            "|dWT| at matched expression, length and publication count")

    # --- F4E  49 tests, none above Bonferroni
    k = pd.read_csv(NEU / "ws2_klassentest.csv")
    kk = k[k.klasse != "_GESAMT"].copy()
    kk["rechenbar"] = kk.status == "ok"
    kk["ueber_grenze"] = (kk.OR_gematcht > kk.OR_mde80) & kk.rechenbar
    alpha = 0.05 / int(kk.rechenbar.sum())
    kk["ueber_bonferroni"] = (kk.p < alpha) & (kk.OR_gematcht > 1) & kk.rechenbar
    schreib(kk[["panel", "klasse", "n_klasse_panel", "n_satz", "beobachtet",
                "null_mittel", "OR_gematcht", "OR_mde80", "z", "p",
                "rechenbar", "ueber_grenze", "ueber_bonferroni"]],
            "F4F_mechanism_classes",
            f"{len(kk)} tests, {int(kk.rechenbar.sum())} computable, "
            f"alpha {alpha:.5f}, above Bonferroni "
            f"{int(kk.ueber_bonferroni.sum())}, above own limit "
            f"{int(kk.ueber_grenze.sum())}")

    # --- F4F  the pooled test per gene panel
    g = k[k.klasse == "_GESAMT"][
        ["panel", "n_satz", "n_panel", "beobachtet", "null_mittel",
         "OR_gematcht", "OR_mde80", "z", "p"]].copy()
    g["ueber_grenze"] = g.OR_gematcht > g.OR_mde80
    schreib(g, "F4F_pooled_per_panel", "pooled enrichment, 7 panels")


# =============================================================================
# FIGURE 5 -- both layers meet at the prehypertrophic transition
# =============================================================================
def f5() -> None:
    log("\nFig. 5 -- the prehypertrophic transition")
    ZONEN = ["MesCond", "ChondroProg", "RestingChon", "ProlifChon",
             "PrehyperChon", "HyperChon"]
    rang = {z: i + 1 for i, z in enumerate(ZONEN)}

    # --- F5A  disease genes against the module, per zone
    p3 = pd.read_csv(NEU / "ws4_p3_panel_vs_modul.csv")
    p3["zonenrang"] = p3.zone.map(rang)
    p3["spezimen"] = p3.probe.str.extract(r"^(Pcw[0-9.]+(?:_s[0-9])?)")
    a = (p3.groupby(["zone", "zonenrang"])
         .agg(n_proben=("probe", "size"), kontrast_median=("kontrast", "median"),
              kontrast_lo=("kontrast", lambda s: s.quantile(0.25)),
              kontrast_hi=("kontrast", lambda s: s.quantile(0.75)),
              z_median=("z", "median"), mde80_median=("mde80", "median"))
         .reset_index().sort_values("zonenrang"))
    schreib(a, "F5A_panel_vs_module_per_zone",
            "median and quartiles per zone")
    schreib(p3[["zone", "zonenrang", "probe", "spezimen", "kontrast", "mde80",
                "z", "p"]],
            "F5A_panel_vs_module_per_sample", "raw points for F5A")

    # --- F5B  both curves overlaid
    mw = pd.read_csv(NEU / "ws4_modulwert_je_probe.csv")
    mw["zonenrang"] = mw.zone.map(rang)
    kurve_m = (mw.groupby(["zone", "zonenrang"])
               .agg(median=("kontrast", "median"),
                    lo=("kontrast", lambda s: s.quantile(0.25)),
                    hi=("kontrast", lambda s: s.quantile(0.75)),
                    n=("probe", "size")).reset_index())
    kurve_m["kurve"] = "programme (173 genes)"
    kurve_p = a.rename(columns={"kontrast_median": "median",
                                "kontrast_lo": "lo", "kontrast_hi": "hi",
                                "n_proben": "n"})[
        ["zone", "zonenrang", "median", "lo", "hi", "n"]].copy()
    kurve_p["kurve"] = "PanelApp 309 minus programme"
    b = pd.concat([kurve_m, kurve_p], ignore_index=True).sort_values(
        ["kurve", "zonenrang"])
    schreib(b, "F5B_both_curves", "shared prehypertrophic maximum")


# =============================================================================
# FIGURE 6 -- levels and detection limits
# =============================================================================
def f6() -> None:
    log("\nFig. 6 -- levels and detection-limit book (table + forest)")

    # Figure F6 draws a forest and range plot (20_figures_main.R::bau_f6)
    # from the NUMERIC file F6_levels_forest.csv; the text table goes into the
    # supplement as Table S14. Both files carry the same levels in the same
    # order. The numbers of the numeric file are derived programmatically
    # from the same source files as the panels -- nothing is entered by hand.
    z = pd.DataFrame([
        # level, calibration, calibration value, limit, finding, value, verdict, status
        dict(ebene="Perturbation datasets, transcriptome (18 series)",
             eichung="2/18 passed (z>=2)", grenze="per dataset",
             befund="18/18 above own limit, z +5.25..+13.10", wert="z +5.25..+13.10",
             urteil="carries", art="explorativ"),
        dict(ebene="Donors within one experiment",
             eichung="7/14 cells passed", grenze="S1 0.344",
             befund="S1 0.349, z +3.00 (+4.51 without\nmodule-forming cells)",
             wert="z +3.00",
             urteil="carries", art="vorregistriertes Follow-up"),
        dict(ebene="Human fetal growth plate, in vivo (E-MTAB-8813)",
             eichung="64/66 samples above own limit,\nmedian z +10.92",
             grenze="rho 0.274", befund="rho 0.456, z +4.80", wert="z +4.80",
             urteil="carries", art="explorativ"),
        # NOTE: this level also has a rank-based statistic, which is larger.
        # This figure carries the concordance throughout, because only the
        # concordance has an MDE80 (project rule 3).
        dict(ebene="Chromatin H3K27ac (GSE129031)",
             eichung="passed", grenze="C 0.578..0.620",
             befund="above own limit 9/9,\nz +2.97..+3.63",
             wert="z +2.97..+3.63", urteil="carries", art="explorativ"),
        dict(ebene="Adipogenic axis, ATAC (GSE332758)",
             eichung="4/4 windows (z +3.73..+5.01)",
             grenze="per window",
             befund="above limit in 3/4 windows (H1 null),\n2/4 (background), z +2.38..+4.76",
             wert="z +2.38..+4.76", urteil="carries", art="explorativ"),
        dict(ebene="Osteogenic axis, ATAC (GSE332758)",
             eichung="0/4 windows (z -0.54..+1.66)",
             grenze="per window", befund="above limit in 4/4, z +3.50..+4.51",
             wert="z +3.50..+4.51",
             urteil="decoupling observation, not a calibrated module result",
             art="explorativ"),
        dict(ebene="Lineage contrast, ATAC (GSE332758)",
             eichung="0/4 (calibration L,\nz +0.93..+2.26)",
             grenze="per window", befund="0/4 above limit, z +0.85..+1.65",
             wert="z +0.85..+1.65",
             urteil="NOT MEASURABLE - no calibrated axis", art="explorativ"),
        dict(ebene="Promoter methylome 27K (GSE33896)",
             eichung="markers weak", grenze="not reachable",
             befund="null on both axes, z -0.46 / -1.76", wert="z -0.46 / -1.76",
             urteil="carries nothing", art="explorativ"),
        dict(ebene="Promoter methylome 450K (GSE129266)",
             eichung="passed", grenze="K 0.312", befund="K 0.298, below own limit", wert="below own limit", urteil="carries nothing", art="explorativ"),
        dict(ebene="External cohorts, fixed 173-gene set (Table S12)",
             eichung="no per-study control", grenze="0.679",
             befund="pooled share 0.682, z +2.84, p 0.0015", wert="z +2.84, p 0.0015",
             urteil="carries as synthesis only (0/11 studies reach own limit)",
             art="vorregistriertes Follow-up"),
        dict(ebene="Human genetics, 7 gene panels",
             eichung="OR 16.8-50.9 (lineage markers),\n3.70-5.75 (secretion anchor)",
             grenze="OR 1.59 (GWAS)", befund="pooled OR 1.00 (GWAS), 0/40 tests\nover Bonferroni",
             wert="0/40 over Bonferroni",
             urteil="good negative result", art="konfirmatorisch"),
        dict(ebene="Disease genes on the dynamics axis",
             eichung="programme z +18.10 in the same run",
             grenze="0.016-0.113 |dWT|", befund="PA309 z -0.71, n.s.", wert="n.s.", urteil="good negative result", art="explorativ"),
        dict(ebene="Gene constraint, publication-matched",
             eichung="GWAS z -4.67", grenze="0.024-0.139 LOEUF",
             befund="PA309 z +1.83, n.s.", wert="n.s.",
             urteil="good negative result", art="explorativ"),
        dict(ebene="Undifferentiated state as a predictor (day zero)",
             eichung="n = 7 calibrated cells", grenze="not reachable",
             befund="z -1.68, p 0.084", wert="p 0.084", urteil="fails",
             art="explorativ"),
    ])
    schreib(z, "TS14_levels_book",
            "one row per level, including the levels without a finding "
            "(supplement table to forest panel F6)")

    # ---- The numeric file for the forest figure F6. Every value is read
    # from the same source files as the panels; the columns are:
    #   einheit       the facet of the figure ("unit-free statistic",
    #                 "odds ratio", "z (contrasts)")
    #   est_min/max   the estimate range per level (a point if min == max)
    #   grenze_min/max  the level's own detection limit in the same unit
    #                   (NaN = not expressible in that unit; the legend then
    #                   carries the limit as text)
    #   null_min/max  optionally the null mean +/- SD (donor level only)
    # No estimate is converted: every row stands in its own unit, readable on
    # the axis of its facet.
    b2 = pd.read_csv(AUS / "F2B_module_per_dataset.csv")
    st = pd.read_csv(RES / "statistik.csv")
    st1 = st[st.groesse == "Programm"].iloc[0]
    t3 = pd.read_csv(AUS / "F3C_trend_test_donor.csv")
    t3m = t3[t3.quantity == "module (173 genes)"].iloc[0]
    h3 = pd.read_csv(ERG / "B_atac" / "B3_GSE129031_modultest.csv")
    h3 = h3[(h3.achse == "chondrogen") & (h3.geeicht.astype(str) == "True")]
    at = pd.read_csv(AUS / "F2E_atac_per_axis.csv")
    at = at[at["null_model"] == "H1 baseline-stratified"]
    m27 = pd.read_csv(ERG / "A_dnam" / "A_dnam_GSE33896_modultest.csv")
    m45 = pd.read_csv(ERG / "A_dnam" / "A_dnam450_GSE129266_modultest.csv"
                      ).iloc[0]
    s12 = pd.read_csv(PDAT / "f6_s12_fixed173_summary.csv").iloc[0]
    gp4 = pd.read_csv(AUS / "F4F_pooled_per_panel.csv")
    gp4 = gp4[gp4.panel == "height GWAS"].iloc[0]
    d4d = pd.read_csv(AUS / "F4E_dynamics_axis.csv")
    c42 = pd.read_csv(AUS / "F4D_constraint_publication_matched.csv")
    d0 = pd.read_csv(AUS / "S8A_day_zero_falls.csv")
    d0p = d0[d0.variant == "primary (7 calibrated cells)"].iloc[0]

    def zeile(ebene, einheit, est_min, est_max, gmin=np.nan, gmax=np.nan,
              urteil="carries", n_min=np.nan, n_max=np.nan):
        return dict(ebene=ebene, einheit=einheit, est_min=float(est_min),
                    est_max=float(est_max), grenze_min=float(gmin),
                    grenze_max=float(gmax), urteilklasse=urteil,
                    null_min=float(n_min), null_max=float(n_max))

    zf = pd.DataFrame([
        zeile("Perturbation datasets, transcriptome (18 series)", "z (contrasts)",
              b2.concordance_z.min(), b2.concordance_z.max(),
              urteil="carries"),
        zeile("Donors within one experiment", "unit-free statistic",
              st1.S1_beobachtet, st1.S1_beobachtet, st1.S1_mde80,
              st1.S1_mde80, "carries",
              st1.S1_null_mittel - st1.S1_null_sd,
              st1.S1_null_mittel + st1.S1_null_sd),
        zeile("Human fetal growth plate, in vivo (E-MTAB-8813)", "unit-free statistic",
              t3m.rho, t3m.rho, t3m.detection_limit_rho, t3m.detection_limit_rho),
        zeile("Chromatin H3K27ac (GSE129031)", "unit-free statistic",
              h3.konkordanz.min(), h3.konkordanz.max(),
              h3.konkordanz_mde80.min(), h3.konkordanz_mde80.max()),
        zeile("Adipogenic axis, ATAC (GSE332758)", "unit-free statistic",
              at[at.axis == "adipogenic"].concordance.min(),
              at[at.axis == "adipogenic"].concordance.max(),
              at[at.axis == "adipogenic"].concordance_detection_limit.min(),
              at[at.axis == "adipogenic"].concordance_detection_limit.max()),
        zeile("Osteogenic axis, ATAC (GSE332758)", "unit-free statistic",
              at[at.axis == "osteogenic"].concordance.min(),
              at[at.axis == "osteogenic"].concordance.max(),
              at[at.axis == "osteogenic"].concordance_detection_limit.min(),
              at[at.axis == "osteogenic"].concordance_detection_limit.max(),
              urteil="observation only"),
        zeile("Lineage contrast, ATAC (GSE332758)", "unit-free statistic",
              at[at.axis == "lineage contrast"].concordance.min(),
              at[at.axis == "lineage contrast"].concordance.max(),
              at[at.axis == "lineage contrast"].concordance_detection_limit.min(),
              at[at.axis == "lineage contrast"].concordance_detection_limit.max(),
              urteil="not measurable"),
        zeile("Promoter methylome 27K (GSE33896)", "z (contrasts)",
              m27.konkordanz_z.min(), m27.konkordanz_z.max(),
              urteil="carries nothing"),
        zeile("Promoter methylome 450K (GSE129266)", "unit-free statistic",
              m45.konkordanz, m45.konkordanz, m45.konkordanz_mde80,
              m45.konkordanz_mde80, urteil="carries nothing"),
        zeile("External cohorts, fixed 173-gene set (Table S12)", "unit-free statistic",
              s12.pooled_share, s12.pooled_share, s12.mde80, s12.mde80,
              urteil="observation only"),
        zeile("Human genetics, 7 gene panels", "odds ratio",
              gp4.odds_ratio_matched, gp4.odds_ratio_matched, gp4.odds_ratio_detection_limit, gp4.odds_ratio_detection_limit,
              urteil="good negative"),
        zeile("Disease genes on the dynamics axis", "z (contrasts)",
              d4d[d4d.gene_set == "PanelApp 309"].z.iloc[0],
              d4d[d4d.gene_set == "PanelApp 309"].z.iloc[0], urteil="good negative"),
        zeile("Gene constraint, publication-matched", "z (contrasts)",
              c42[c42.gene_set == "PanelApp 309"].z.iloc[0],
              c42[c42.gene_set == "PanelApp 309"].z.iloc[0], urteil="good negative"),
        zeile("Undifferentiated state as a predictor (day zero)", "z (contrasts)",
              d0p.amplitude_z, d0p.amplitude_z, urteil="fails"),
    ])
    # Order and verdicts as in the text table; the check that both files carry
    # the same levels is in 10_check_numbers.py.
    schreib(zf, "F6_levels_forest",
            "one row per level with estimate range and own limit "
            "(basis of the forest panel F6)")


# =============================================================================
# SUPPLEMENT
# =============================================================================
def supplement() -> None:
    log("\nSupplement")

    # --- S2  every calibration and every limit in one place. The row labels
    # of figure S2 carry the display name including the accession, not the
    # machine key.
    e18 = pd.read_csv(ERG / "M_kalibrierung" / "eichung_achtzehn.csv")
    e18 = e18.assign(
        ebene="18 perturbation\ndatasets",
        einheit=e18.datensatz.map(lambda v: _display.DATENSATZ.get(v, v)))[
        ["ebene", "einheit", "arm", "kontrast", "mde80", "z", "p", "status",
         "bestanden"]]
    ez = pd.read_csv(ERG / "M_donoren" / "eichung.csv")
    ez = ez.assign(
        ebene="14 donor cells",
        einheit=ez.zelle.map(lambda v: _display.ZELLE.get(v, v)))[
        ["ebene", "einheit", "achse", "kontrast", "mde80", "z", "p", "status",
         "bestanden"]].rename(columns={"achse": "arm"})
    ea = pd.read_csv(ERG / "B_atac" / "B_atac_eichung_je_achse.csv")
    # "D" is the axis calibration (undifferentiated against differentiated),
    # "L" the lineage calibration (osteogenic against adipogenic) -- spelled
    # out rather than coded.
    eart = ea.eichung.map({"D": "differentiation", "L": "lineage"})
    ea = ea.assign(ebene="ATAC\n(GSE332758)",
                   einheit="window " + ea.fenster + ", " + eart)[
        ["ebene", "einheit", "achse", "kontrast", "mde80", "z", "p",
         "bestanden"]].rename(columns={"achse": "arm"})
    ea["status"] = "ok"
    s2 = pd.concat([e18, ez, ea], ignore_index=True)
    schreib(s2, "S2A_all_calibrations",
            "every calibration in the project, one row")

    # --- S3  external triangulation, three versions
    fixed = pd.read_csv(PDAT / "f6_s12_fixed173_summary.csv")
    pooled = pd.read_csv(PDAT / "f6_pooled_summary.csv")
    study = pd.read_csv(PDAT / "f6_study_level_summary.csv")
    s3 = pd.DataFrame([
        dict(fassung="S12, fixed 173-gene set (primary)",
             ebene="observations", n=int(fixed.n_studies.iloc[0]),
             statistik=float(fixed.pooled_share.iloc[0]),
             z=float(fixed.z.iloc[0]), p=float(fixed.p_perm_two_sided.iloc[0]),
             mde80=float(fixed.mde80.iloc[0]), hat_grenze=True),
        dict(fassung="pooled directional replication",
             ebene="observations", n=int(pooled.n_datasets.iloc[0]),
             statistik=float(pooled.pooled_share.iloc[0]),
             z=float(pooled.z.iloc[0]), p=float(pooled.p_perm_two_sided.iloc[0]),
             mde80=float(pooled.mde80.iloc[0]), hat_grenze=True),
        dict(fassung="study-level synthesis (mean of per-study z)",
             ebene="studies", n=int(study.n_studies.iloc[0]),
             statistik=float(study.observed_stat.iloc[0]),
             z=float(study.z.iloc[0]), p=float(study.p_upper.iloc[0]),
             mde80=np.nan, hat_grenze=False),
    ])
    schreib(s3, "S3A_external_triangulation",
            "three versions; only the first has a detection limit")
    by = pd.read_csv(PDAT / "f6_s12_fixed173_by_study.csv")
    schreib(by, "S3B_triangulation_per_study",
            f"none of the {len(by)} studies reaches its own MDE80")

    # --- S5  robustness and leave-one-out
    schreib(pd.read_csv(ERG / "M_donoren" / "auslassung.csv"),
            "S5A_leave_one_out", "leave-one-out over the 7 calibrated cells")
    schreib(pd.read_csv(ERG / "M_donoren" / "zirkularitaet.csv"),
            "S5A_circularity", "without the module-forming cells")
    schreib(pd.read_csv(ERG / "M_kalibrierung" / "eichung_empfindlichkeit.csv"),
            "S5C_calibration_sensitivity",
            "preregistered vs markers reachable only in vitro")

    # --- S6  orthogonal layers
    mt = pd.read_csv(ERG / "B_atac" / "B_atac_modultest_final.csv")
    schreib(mt, "S6A_atac_complete", "ATAC, all windows, axes, nulls")
    schreib(pd.read_csv(ERG / "B_atac" / "B3_GSE129031_modultest.csv"),
            "S6B_h3k27ac", "H3K27ac GSE129031")
    schreib(pd.read_csv(ERG / "A_dnam" / "A_dnam_GSE33896_modultest.csv"),
            "S6C_methylome_27k", "GSE33896 -- null finding with a limit")
    schreib(pd.read_csv(ERG / "A_dnam" / "A_dnam450_GSE129266_modultest.csv"),
            "S6D_methylome_450k", "GSE129266 -- below its own limit")

    # --- S8  day zero falls, and the publication matching
    d0s = pd.read_csv(NEU / "ws3_p1_p2_korrelation_z.csv")
    # Presentation layer: variant names for the English panel.
    d0s["variante"] = d0s["variante"].replace({
        "primaer_n7_geeicht": "primary (7 calibrated cells)",
        "sensitivitaet_n14_alle": "sensitivity (all 14 cells)"})
    schreib(d0s, "S8A_day_zero_falls", "naive state as predictor")
    schreib(pd.read_csv(NEU / "ws3_p4_eichung_naivunterschied.csv"),
            "S8A_day_zero_calibration",
            "the 7 failed cells at day 0")
    u = pd.read_csv(NEU / "ws1_p6_untersuchungsintensitaet.csv")
    schreib(u, "S8C_study_intensity",
            "gene2pubmed, retrieved 2026-08-22")
    p6 = pd.read_csv(NEU / "ws1_p6_publikationsmatching.csv")
    schreib(p6, "S8D_publication_matching",
            "P1, P2, P3 each without and with publications")


def main() -> None:
    log("=" * 78)
    log("10_panel_data_main.py -- one CSV per panel to figures/data/")
    log("=" * 78)
    # supplement() BEFORE f6(): the levels book reads S8A_day_zero_falls.csv,
    # and supplement() writes it. Previously f6() read the file of the
    # previous run -- unnoticeable as long as nothing changed.
    f1(); f2(); f3(); f4(); f5(); supplement(); f6()
    n = len(list(AUS.glob("*.csv")))
    log(f"\n{n} panel CSVs in {AUS}")
    log(f"\nnumpy {np.__version__} | pandas {pd.__version__} | "
        f"python {sys.version.split()[0]}")
    (RES / "panel_data_log.txt").write_text("\n".join(LOG) + "\n",
                                            encoding="utf-8")


if __name__ == "__main__":
    main()

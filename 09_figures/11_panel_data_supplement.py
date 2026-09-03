# -*- coding: utf-8 -*-
"""
11_panel_data_supplement.py -- panel data for the supplementary figures and for the
                         supplementary tables.

Purpose  Adds to `10_panel_data_main.py` the parts that are only needed in the
         supplement: the scale critique (S1), patient against control (S4),
         robustness (S5), the screen in detail (S7), day zero (S8) -- and the
         supplementary tables S1 to S14.

Rule     As in 10_panel_data_main.py: no second implementation of any metric.
         Files are only read and reshaped.

Inputs   derived_data/**, derived_data/followup/**,
         derived_data/manuscript/**, derived_data/reference_tables/**
Outputs  figures/data/S*.csv, figures/data/TS*.csv
         results/supplementdaten_log.txt
Runtime  a few seconds
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import _display  # noqa: E402  presentation layer: everything shown is English

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
ERG = WURZEL / "derived_data"
NEU = WURZEL / "derived_data" / "followup"
PDAT = WURZEL / "derived_data" / "manuscript"
TAB = WURZEL / "derived_data" / "reference_tables"
RES = WURZEL / "results"
AUS = WURZEL / "figures" / "data"
AUS.mkdir(parents=True, exist_ok=True)

SPEZIMEN = r"^(Pcw[0-9.]+(?:_s[0-9])?)"
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def schreib(df: pd.DataFrame, name: str, was: str) -> None:
    # Display values and column headers go out in English.
    df = _display.englisch(df)
    df.to_csv(AUS / (name + ".csv"), index=False)
    log("  %-38s %5d rows x %2d cols  %s"
        % (name, len(df), df.shape[1], was))


# =============================================================================
# S1 -- scale critique, unabridged
# =============================================================================
def s1() -> None:
    log("\nS1 -- scale critique")
    schreib(pd.read_csv(PDAT / "f5_ausgangsexpression.csv"),
            "S1A_baseline_expression",
            "median z vs median baseline shift per GO set")
    schreib(pd.read_csv(PDAT / "f3_empfindlichkeit.csv"),
            "S1B_sensitivity",
            "which additive offset the statistic recovers")
    schreib(pd.read_csv(PDAT / "f3_neutralkontrast.csv"),
            "S1C_neutral_contrast", "null rate of the statistic, self test")
    schreib(pd.read_csv(PDAT / "f4_konvergenzkurve.csv"),
            "S1D_convergence_curve",
            "convergent genes vs noise floor -- the step function")


# =============================================================================
# S4 -- patient against control
# =============================================================================
def s4() -> None:
    log("\nS4 -- patient vs control (downgraded)")
    e = pd.read_csv(PDAT / "m3_eichung.csv")
    schreib(e[["gse", "entitaet", "gewebe", "satz_a", "n_patienten",
               "n_kontrollen", "kontrast", "mde80", "z", "p", "bestanden"]],
            "S4A_patient_calibration",
            "gate B is withdrawn -- two cohorts were calibrated against NAIVE")
    pt = pd.read_csv(PDAT / "m3_patienten.csv")
    pt = pt[pt.satz == "Programm"].copy()
    pt["ueber_grenze"] = pt.konkordanz > pt.konkordanz_mde80
    schreib(pt[["gse", "entitaet", "gewebe", "n_patienten", "n_kontrollen",
                "n", "konkordanz", "konkordanz_mde80", "konkordanz_z",
                "konkordanz_p", "eichung_bestanden", "ueber_grenze"]],
            "S4B_patient_concordance",
            "%d of %d cohorts above their own limit"
            % (int(pt.ueber_grenze.sum()), len(pt)))


# =============================================================================
# S5 -- robustness and leave-one-out
# =============================================================================
def s5() -> None:
    log("\nS5 -- robustness (supplement)")
    schreib(pd.read_csv(ERG / "M_donoren" / "selbsttest.csv"),
            "S5B_self_test", "null rate of the metric in the donor set-up")


# =============================================================================
# S7 and S8 -- the screen in detail, and day zero
# =============================================================================
def s7_s8() -> None:
    log("\nS7 / S8 -- screening and day zero")
    schreib(pd.read_csv(ERG / "M_diagnosen" / "sichtung.csv"),
            "S7B_screen_diagnoses_full", "all 127 hand checks")
    schreib(pd.read_csv(ERG / "M_donoren" / "sichtung_hand.csv"),
            "S7C_screen_by_design", "screening along the design")
    schreib(pd.read_csv(NEU / "ws3_zellen_tabelle.csv"),
            "S8B_cells_day_zero", "module value at day 0 per cell")


# =============================================================================
# the supplementary tables
# =============================================================================
def tabellen() -> None:
    log("\nSupplementary tables TS1..TS13")

    schreib(pd.read_csv(TAB / "S7_kohorte_18_datensaetze.csv"),
            "TS1_eighteen_datasets", "accessions and design")
    schreib(pd.read_csv(TAB / "S1_sichtung_alle_datensaetze.csv"),
            "TS2_screen_exclusion_codes", "all screened series")

    e18 = pd.read_csv(ERG / "M_kalibrierung" / "eichung_achtzehn.csv")
    e18 = e18.assign(ebene="dataset", einheit=e18.datensatz)
    ez = pd.read_csv(ERG / "M_donoren" / "eichung.csv")
    ez = ez.assign(ebene="donor cell", einheit=ez.zelle, arm=ez.achse)
    sp = ["ebene", "einheit", "arm", "n_a", "n_b", "kontrast", "null_mittel",
          "null_sd", "z", "p", "mde80", "status", "bestanden"]
    schreib(pd.concat([e18[sp], ez[sp]], ignore_index=True),
            "TS3_calibrations",
            "study level %d/%d, donor level %d/%d passed"
            % (int(e18.bestanden.sum()), len(e18),
               int(ez.bestanden.sum()), len(ez)))

    g = pd.read_csv(TAB / "S5_konvergente_gene.csv")
    schreib(g, "TS4_module_genes",
            "%d genes, %d up and %d down"
            % (len(g), int((g.ri > 0).sum()), int((g.ri < 0).sum())))

    a = pd.read_csv(NEU / "ws4_proben_je_zone.csv")
    mw = pd.read_csv(NEU / "ws4_modulwert_je_probe.csv")
    pk = pd.read_csv(NEU / "ws4_positivkontrolle_je_probe.csv")
    pk = pk[pk.vergleich == "chondrogen_vs_naiv"]
    t = a.merge(mw[["zone", "probe", "kontrast", "mde80", "z", "p"]],
                on=["zone", "probe"], how="left")
    t = t.merge(pk[["zone", "probe", "kontrast", "mde80", "z", "p"]],
                on=["zone", "probe"], how="left",
                suffixes=("_module", "_positivkontrolle"))
    t["spezimen"] = t.probe.str.extract(SPEZIMEN, expand=False)
    schreib(t, "TS5_in_vivo_per_zone_and_specimen",
            "cell count, module value and positive control per (zone, sample)")

    schreib(pd.read_csv(NEU / "ws2_klassifikation_panels.csv"),
            "TS6_panels_and_classes", "mapping of panel genes to classes")

    # TS7 -- everything that has a detection limit, in one table
    teile = []
    konk = [("module per dataset", NEU / "ws6_p1p2_modul_je_datensatz.csv"),
            ("ATAC module test", ERG / "B_atac" / "B_atac_modultest_final.csv"),
            ("H3K27ac module test",
             ERG / "B_atac" / "B3_GSE129031_modultest.csv"),
            ("methylome 27K", ERG / "A_dnam" / "A_dnam_GSE33896_modultest.csv"),
            ("methylome 450K",
             ERG / "A_dnam" / "A_dnam450_GSE129266_modultest.csv")]
    for name, quelle in konk:
        d = pd.read_csv(quelle)
        if "konkordanz" not in d.columns:
            log("     skipped (no concordance column): " + name)
            continue
        lab = [c for c in ["datensatz", "fenster", "achse", "null", "ebene"]
               if c in d.columns]
        sp2 = [c for c in ["konkordanz", "konkordanz_null",
                           "konkordanz_null_sd", "konkordanz_z",
                           "konkordanz_p", "konkordanz_mde80"]
               if c in d.columns]
        x = d[lab + sp2].copy()
        x.insert(0, "block", name)
        if "konkordanz_mde80" in x.columns:
            x["ueber_grenze"] = x.konkordanz > x.konkordanz_mde80
        teile.append(x)
    zteste = [("gene-set matching", NEU / "ws1_p6_publikationsmatching.csv",
               "z", "mde80_delta"),
              ("mechanism classes", NEU / "ws2_klassentest.csv", "z",
               "OR_mde80"),
              ("complementarity", NEU / "ws7_komplementaritaet.csv", "z",
               "OR_mde80")]
    for name, quelle, zsp, msp in zteste:
        d = pd.read_csv(quelle)
        lab = [c for c in ["satz", "lauf", "panel", "klasse", "seite", "gegen"]
               if c in d.columns]
        x = d[lab + [zsp, msp, "p"]].copy()
        x.columns = lab + ["statistik_z", "grenze", "p"]
        x.insert(0, "block", name)
        teile.append(x)
    schreib(pd.concat(teile, ignore_index=True), "TS7_all_statistics",
            "every test statistic with its detection limit")

    # TS8 -- the preregistrations and their status
    zeilen = []
    for f in sorted((WURZEL / "preregistrations").glob("*.md")):        # Some of the preregistrations are stored in cp1252 (Windows editor).
        # errors="replace" turned that into U+FFFD and made the titles
        # untranslatable, so utf-8 is tried first and cp1252 second.
        try:
            kopf = f.read_text(encoding="utf-8")[:400]
        except UnicodeDecodeError:
            kopf = f.read_text(encoding="cp1252")[:400]
        art = ("protocol" if f.name.startswith("PROTOKOLL")
               else "preregistration" if f.name.startswith("PRAEREG")
               else "document")
        zeilen.append({"datei": f.name, "art": art,
                       "titel": titel_aus_kopf(kopf)})
    schreib(pd.DataFrame(zeilen), "TS8_preregistrations",
            "%d documents in preregistrations/" % len(zeilen))

    # TS9 -- Figure 2D in full: the narrow set, the broad set and the second
    # independent source side by side, with the source and the members of
    # each set. This is where the number stands that the decision rule
    # removed from the running text (cell-cycle exit, case b).
    v2 = pd.read_csv(RES / "gensaetze_v2_anreicherung.csv")
    u = pd.read_csv(RES / "gensaetze_v2_urteil.csv")[["kategorie", "fall"]]
    v2 = v2.merge(u, on="kategorie", how="left")
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
    v2 = v2[["category", "set", "quelle", "bestandteile", "k", "n_modul",
             "n_satz_gesamt", "n_satz_im_hg", "n_hintergrund", "OR", "lo",
             "hi", "p", "p_bonferroni", "ki_schliesst_1_aus", "fall"]]
    schreib(v2, "TS9_gene_sets_v2",
            "6 categories x up to 3 sets, decision rule "
            "06_orthogonal_layers/61_gene_set_enrichment.R")

    # TS9b -- which module genes carry which category, per set
    schreib(pd.read_csv(RES / "gensaetze_v2_modulgene.csv"),
            "TS9b_gene_sets_v2_module_genes",
            "module genes per category and set")

    # TS10 -- the three-way decomposition, each of the 18 data sets
    # individually, under BOTH gene maps. Preregistered in
    # preregistrations/PRAEREG_F2F.md.
    z18 = pd.read_csv(RES / "zerlegung_achtzehn.csv")
    schreib(z18, "TS10_decomposition_eighteen",
            "18 datasets x 2 gene maps, with rule application per point")

    # TS11 -- does a minority of large excursions carry the in vivo trend?
    # The decision rule is in the header of
    # 07_in_vivo_growth_plate/13_fetal_gene_decomposition.py.
    gz = pd.read_csv(RES / "invivo_genzerlegung.csv")
    schreib(gz, "TS11_in_vivo_gene_decomposition",
            "trend after removing the largest |Delta|, 0/5/10/20/30 %")

    # TS11b -- the ranking itself, every module gene with its Delta
    gk = pd.read_csv(NEU / "ws4_p2_gen_konkordanz.csv").copy()
    gk["abs_delta"] = gk.delta_hyper_minus_mes.abs()
    gk = gk.sort_values("abs_delta", ascending=False)
    gk["rang_abs_delta"] = range(1, len(gk) + 1)
    schreib(gk, "TS11b_in_vivo_gene_ranking",
            "173 module genes by |Delta| hyper minus mesCond")

    # TS12 -- the primary publications of every GEO series used. The two
    # bibliographic tables are frozen in derived_data/reference_tables/
    # (geo_primary_publications*.csv); no fetch step runs here.
    # TS13 -- does the fetal in vivo trend hang on the single hypertrophic
    # point? The decision rule is in the header of
    # 07_in_vivo_growth_plate/14_hypertrophic_zone_sensitivity.py.
    hz = pd.read_csv(RES / "invivo_hz_empfindlichkeit.csv")
    schreib(hz, "TS13_in_vivo_hypertrophic_sensitivity",
            "trend without HyperChon (rule) and without Hyper+Prehyper "
            "(descriptive)")

    # TS12 draws on two sources: what GEO links itself
    # (geo_primary_publications.csv) and the four series whose publication
    # GEO does NOT link and which were documented by hand. Two of them are
    # published and GEO simply does not link them; two are genuinely
    # unpublished and are cited by their accession, which the Cell Press
    # reference policy allows.
    pa = pd.read_csv(WURZEL / "derived_data" / "reference_tables" /
                     "geo_primary_publications.csv")
    hand = pd.read_csv(WURZEL / "derived_data" / "reference_tables" /
                       "geo_primary_publications_manual.csv")
    pa = pa[~pa.gse.isin(hand.gse)]
    pa = pd.concat([pa, hand], ignore_index=True).sort_values("gse")
    # Hand corrections where the automatic GEO/PubMed pull stays incomplete:
    # MDPI, Life Science Alliance and AAAS put an article number where PubMed
    # leaves the pages field empty, and the Chu et al. volume exists only at
    # the journal itself. Keyed by DOI, so a re-pull (code/28_geo_primary_
    # publications.py) cannot silently lose them.
    korrekturen = {
        # Hernandez 2021, Biology 10(8), article 802
        "10.3390/biology10080802": {"seiten": "802"},
        # Schoenmaker 2023, Int J Mol Sci 24(7), article 6822
        "10.3390/ijms24076822": {"seiten": "6822"},
        # Sun 2024, Life Sci Alliance 7(3), eLocator e202302219
        "10.26508/lsa.202302219": {"seiten": "e202302219"},
        # Chu 2026, Sci Transl Med 18, eLocator eadw3590 (issue 845)
        "10.1126/scitranslmed.adw3590": {"band": "18",
                                         "seiten": "eadw3590"},
    }
    for doi, fix in korrekturen.items():
        treffer = pa.doi == doi
        if not treffer.any():
            raise ValueError(f"reference correction without a row: {doi}")
        for spalte, wert in fix.items():
            pa[spalte] = pa[spalte].astype(object)   # float column -> text
            pa.loc[treffer, spalte] = wert
    ohne = int(pa.status.str.startswith("unpublished").sum())
    schreib(pa, "TS12_primary_publications",
            "%d series, %d publication rows, %d unpublished "
            "(cited via the accession)"
            % (pa.gse.nunique(), int(len(pa)) - ohne, ohne))


def titel_aus_kopf(kopf: str) -> str:
    """The title line of a preregistration document.

    The translated documents carry a blockquote note above the title
    ("Translated from the German original ..."); the title is the first
    line that is neither empty nor part of that note.
    """
    for zeile in kopf.splitlines():
        z = zeile.strip()
        if not z or z.startswith(">"):
            continue
        return z.lstrip("# ").strip()
    return ""


def main() -> None:
    log("=" * 78)
    log("11_panel_data_supplement.py -- supplement panels and tables")
    log("=" * 78)
    s1(); s4(); s5(); s7_s8(); tabellen()
    log("\nnumpy %s | pandas %s | python %s"
        % (np.__version__, pd.__version__, sys.version.split()[0]))
    (RES / "supplement_data_log.txt").write_text("\n".join(LOG) + "\n",
                                                 encoding="utf-8")


if __name__ == "__main__":
    main()

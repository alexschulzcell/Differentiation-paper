"""
_display.py -- the presentation layer: everything that goes into a figure
panel or into a supplementary table is ENGLISH.

This module is the ONE place where that translation is defined. It is called
by 10_panel_data_main.py, 11_panel_data_supplement.py and 12_panel_data_second_cohort.py from their
schreib() function, so that figures/data/ carries English display text and
English column headers and nothing else. The upstream analysis scripts are
left untouched; their outputs under derived_data/ keep the short internal
names, and column_glossary.csv (written by 09_figures/_glossary.py from the same
mapping) gives the English equivalent of each of them.

Three conventions are fixed here:

  1. THE ACCESSION RULE: the information first, then the accession in
     parentheses. Canonically "ATAC (GSE332758)", never "GSE332758 (ATAC)"
     and never "ATAC GSE332758". mit_gse() produces that form.

  2. THE CASE OF IDENTIFIERS: spelled out and lower case inside a running
     label ("osteogenic axis"), title case as a standalone identifier
     ("Osteogenic"). No abbreviations -- no "Osteo", no "OSTEO", no
     "osteog", no "chondr".

  3. "undifferentiated" rather than "naive", without exception.
"""
from __future__ import annotations

import re

import pandas as pd


# ------------------------------------------------------- the accession rule
def mit_gse(info: str, akzession: str | None) -> str:
    """The canonical form: information, then the accession in parentheses.

    Without an accession (own data, not yet deposited) the label reads
    "(this study)" -- a blank would be dishonest.
    """
    if akzession is None or not str(akzession).strip() or \
            str(akzession).startswith("GSE_"):
        # No such case is left in the project; should one arise, it should
        # be conspicuous rather than disappear quietly.
        return f"{info} (accession pending)"
    return f"{info} ({akzession})"


# --------------------------------------------------- the eighteen data sets
# Display names of the eighteen perturbation data sets. EVERY row carries its
# accession -- none of the eighteen is unpublished. Three accessions were
# placeholders in the source files until they were resolved:
#   GSE_LAMA5      -> E-MTAB-16566 (ArrayExpress; our own urine-derived stem
#                     cell series, three states: undifferentiated, osteogenic,
#                     chondrogenic)
#   GSE_SERPINA3ch -> GSE247491
#   GSE_MIR181     -> GSE184087
DATENSATZ = {
    "LAMA5-KO (chondr)":             "LAMA5-KO, chondrogenic (E-MTAB-16566)",
    "LAMA5-KO (osteo)":              "LAMA5-KO, osteogenic (E-MTAB-16566)",
    "FN1 FNC123R":                   "FN1 C123R (GSE251698)",
    "FN1 FNC231W":                   "FN1 C231W (GSE251698)",
    "SERPINA3-KD":                   "SERPINA3-KD, chondrogenic (GSE247491)",
    "MIR181A1HG-KD":                 "MIR181A1HG-KD (GSE184087)",
    "LINC01638-KD (GSE227512)":      "LINC01638-KD (GSE227512)",
    "ARSB / MPS VI (GSE218101)":     "ARSB / MPS VI (GSE218101)",
    "ACVR1 / FOP (GSE221128)":       "ACVR1 / FOP (GSE221128)",
    "RNF4-KD (GSE205432)":           "RNF4-KD (GSE205432)",
    "RB1 +/- (GSE245585)":           "RB1 +/- (GSE245585)",
    "Nonunion (GSE226565)":          "Nonunion (GSE226565)",
    "SERPINA3-KD osteo (GSE247528)": "SERPINA3-KD, osteogenic (GSE247528)",
    "ERCC6L2-KD (GSE190542)":        "ERCC6L2-KD (GSE190542)",
    "TP53 LFS (GSE102732)":          "TP53 LFS (GSE102732)",
    "YAP/TAZ-KD (GSE137035)":        "YAP/TAZ-KD (GSE137035)",
    "RB1-mut isogen (GSE145235)":    "RB1-mut isogenic (GSE145235)",
    "RB1-del isogen (GSE145235)":    "RB1-del isogenic (GSE145235)",
}

# ------------------------------------------------------- the fourteen cells
# Donor-resolved cells. The raw names are machine keys
# (SERPINA3_D2_chondr); the display names the cell, the axis and the accession.
ZELLE = {
    "SERPINA3_D1_chondr": "SERPINA3 donor 1, chondrogenic (GSE247491)",
    "SERPINA3_D2_chondr": "SERPINA3 donor 2, chondrogenic (GSE247491)",
    "SERPINA3_D3_chondr": "SERPINA3 donor 3, chondrogenic (GSE247491)",
    "SERPINA3_D1_osteog": "SERPINA3 donor 1, osteogenic (GSE247528)",
    "SERPINA3_D2_osteog": "SERPINA3 donor 2, osteogenic (GSE247528)",
    "SERPINA3_D3_osteog": "SERPINA3 donor 3, osteogenic (GSE247528)",
    "LAMA5_USC_chondr":   "LAMA5-KO USC, chondrogenic (E-MTAB-16566)",
    "LAMA5_USC_osteog":   "LAMA5-KO USC, osteogenic (E-MTAB-16566)",
    "GSE218101_Line1_chondr": "ARSB line 1, chondrogenic (GSE218101)",
    "GSE218101_Line2_chondr": "ARSB line 2, chondrogenic (GSE218101)",
    "GSE218101_Line3_chondr": "ARSB line 3, chondrogenic (GSE218101)",
    "GSE218101_Line4_chondr": "ARSB line 4, chondrogenic (GSE218101)",
    "GSE221128_FOP_chondr":   "ACVR1 / FOP, chondrogenic (GSE221128)",
    "GSE245585_WT1_osteog":   "RB1 +/- WT1, osteogenic (GSE245585)",
}

# ------------------------------------------------------------- accessions
# Three placeholders stood in the source files and are resolved here; the
# `gse` column carries real accessions from this point on.
QUELLE = {
    "suche": "preregistered search",
    "vorrangig+suche": "priority list + preregistered search",
    "nachrangig+suche": "secondary list + preregistered search",
}

AKZESSION = {
    "GSE_LAMA5": "E-MTAB-16566",
    "GSE_SERPINA3ch": "GSE247491",
    "GSE_MIR181": "GSE184087",
}


# ------------------------------------------------------ panels and gene sets
SATZ = {
    "PA309": "PanelApp 309",
    "PA1471": "PanelApp 1471",
    "NOSO": "Nosology (core)",
    "NOSO_BREIT": "Nosology (broad)",
    "KLEIN": "short stature",
    "KLEIN_BREIT": "short stature (broad)",
    "GWAS": "height GWAS",
    "ZELLZYKLUS_NK": "cell cycle (neg. control)",
    "PROGRAMM": "programme (173 genes)",
    # the seven mechanism classes, in the two spellings that occur upstream
    "CILIUM": "cilium",
    "Ziliopathie": "cilium",
    "ECM_STRUKTUR": "ECM structure",
    "Struktur_ECM": "ECM structure",
    "GLYKO_LINKER": "glycosylation / linker",
    "Glykosylierung": "glycosylation / linker",
    "LYSOSOM": "lysosome",
    "Lysosom": "lysosome",
    "SIGNAL_FWBI": "signalling FGF/WNT/BMP/IHH",
    "Signal_FGF_WNT_BMP_IHH": "signalling FGF/WNT/BMP/IHH",
    "TF_DNABINDEND": "transcription factors",
    "TF_DNAbindend": "transcription factors",
    "VESIKEL_SEKRET": "vesicle / secretion",
    "Vesikel_Sekretion": "vesicle / secretion",
}

# ------------------------------------------- diagnosis axes of the screen
# The `achse` column of S7B carries search axes, some of them joined by "|".
# They are therefore translated token by token, not as a whole string.
DIAGNOSE = {
    "OI": "OI",
    "MPS": "MPS",
    "FGFR3": "FGFR3",
    "SHOX": "SHOX",
    "FOP": "FOP",
    "CCD_CMPD": "CCD / CMPD",
    "PSACH_MED": "PSACH / MED",
    "FREITEXT_DYSPLASIE": "free text: dysplasia",
    "FREITEXT_KLEINWUCHS": "free text: short stature",
    "FREITEXT_WACHSTUMSFUGE": "free text: growth plate",
}


def diagnoseachse(v: object) -> object:
    if not isinstance(v, str):
        return v
    return " | ".join(DIAGNOSE.get(t, t) for t in v.split("|"))


# --------------------------------------------------- general value glossary
# Exact matches across all text columns. Values only; the column names are
# translated separately, by SPALTEN further down, so that a value and a
# column name can never be confused for one another.
WERTE = {
    # axes and states
    "osteogen": "osteogenic",
    "adipogen": "adipogenic",
    "chondrogen": "chondrogenic",
    "myogen": "myogenic",
    "naiv": "undifferentiated",
    "differenz": "lineage contrast",
    "OSTEO": "Osteogenic",
    "ADIPO": "Adipogenic",
    "NAIV": "Undifferentiated",
    "OSTEOGEN": "Osteogenic",
    "ADIPOGEN": "Adipogenic",
    "CHONDROGEN": "Chondrogenic",
    "osteogen (Zielachse)": "osteogenic (target axis)",
    "myogen (Linienkontrolle)": "myogenic (lineage control)",
    # Eichungsurteile
    "bestanden": "passed",
    "durchgefallen": "failed",
    "nicht eichbar": "not calibratable",
    "zu wenige Gene": "too few genes",
    "falls": "fails",
    # Gensatzvarianten
    "eng": "narrow",
    "breit": "broad",
    # Sichtungsebenen
    "Entitaet": "Entity",
    "Anlage": "Predisposition",
    # Nullmodelle
    "Hintergrund": "background",
    "H1 basisgeschichtet": "H1 baseline-stratified",
    "basisgeschichtet": "H1 baseline-stratified",
    # Cohort identifiers turned round to follow the accession rule: the assay
    # first, then the accession in parentheses.
    "GSE332758 (ATAC)": "ATAC (GSE332758)",
    "GSE151311 (ATAC)": "ATAC (GSE151311)",
    "GSE151315 (H3K27ac)": "H3K27ac (GSE151315)",
    "GSE129031 (H3K27ac)": "H3K27ac (GSE129031)",
    # entities and tissues of the patient cohorts (S4)
    "Turner-Syndrom": "Turner syndrome",
    "Enchondromatose (Ollier)": "enchondromatosis (Ollier)",
    "akromele Dysplasie (ADAMTSL2/FBN1)": "acromelic dysplasia (ADAMTSL2 / FBN1)",
    "Osteogenesis imperfecta": "osteogenesis imperfecta",
    "EDS-HT / JHS": "EDS-HT / JHS",
    "FOP (ACVR1)": "FOP (ACVR1)",
    "primaere Osteoblasten": "primary osteoblasts",
    "Knorpel / Enchondrom": "cartilage / enchondroma",
    "dermale Fibroblasten": "dermal fibroblasts",
    "Vollblut": "whole blood",
    "CD14+-Monozyten": "CD14+ monocytes",
    "Fibroblasten": "fibroblasts",
    # quantities of the decomposition and of the trend tests
    "Modul": "module",
    "Modul (173 Gene)": "module (173 genes)",
    "PA309 gegen Modul": "disease genes vs programme",
    "Positivkontrolle chondrogen-naiv":
        "positive control: chondrogenic vs undifferentiated",
    "Zerlegung bestaetigt": "decomposition confirmed",
    "Ebene fuer diese Achse nicht geeicht":
        "level not calibrated for this axis",
    # side, verdict, role and variant labels that arrive in German
    "Programm": "programme",
    "Krankheitsgene_PA309": "disease genes PA309",
    "Krankheitsgene_NOSO": "disease genes NOSO",
    "Krankheitsgene_NOSO_BREIT": "disease genes NOSO broad",
    "unter Schwelle": "below threshold",
    "voll": "full",
    "gefiltert": "filtered",
    # case labels of the external gene-set analysis (61_gene_set_enrichment.R)
    "(a) bestaetigt": "(a) confirmed",
    "(b) abgeschwaecht": "(b) weakened",
    "(c) umgekehrt oder null": "(c) reversed or null",
    "war schon null": "already null",
    # ATAC verdict label
    "SCHWELLE ERREICHT": "threshold reached",
    # role and variant labels of the in vivo gene decomposition
    "Ausgangswert": "baseline",
    "Regel": "rule",
    "deskriptiv": "descriptive",
    "A_ohne_HyperChon": "A without HyperChon",
    "B_ohne_Hyper_Praehyper": "B without Hyper+Prehyper",
    "alle geeichten Zellen": "all calibrated cells",
    "nur modulbildend": "module-forming cells only",
    "nur NICHT modulbildend": "non-module-forming cells only",
    "-- (alle)": "none left out (reference)",
    # leave-one-out variants and runs (51_ and 50_ already translate these
    # while assembling; this is only a fallback)
    "primaer_n7_geeicht": "primary (7 calibrated cells)",
    "sensitivitaet_n14_alle": "sensitivity (all 14 cells)",
    # exclusion reasons of the hand screen
    "keine Diagnose als Laesionsachse": "no diagnosis as the lesion axis",
    "kein naiver Arm": "no undifferentiated arm",
    "no naive arm": "no undifferentiated arm",
    "no naive arm (A1)": "no undifferentiated arm (A1)",
    "of these: no naive arm (A1)": "of these: no undifferentiated arm (A1)",
    # variable names that appear as a VALUE in a table (S8A)
    "naiv_modulwert_z": "undifferentiated module value (z)",
    "naiv_modulwert_roh": "undifferentiated module value (raw)",
    "naiv_marker_achsendiff": "undifferentiated marker axis difference",
    "naiv_marker_eigen": "undifferentiated own markers",
    "naiv_marker_fremd": "undifferentiated other markers",
    # marker and curve names of the decomposition
    "naive markers": "undifferentiated markers",
    "naive state as a predictor (day zero)":
        "undifferentiated state as a predictor (day zero)",
    # diagnosis axes of the screen (S7B)
    "FREITEXT_WACHSTUMSFUGE": "free text: growth plate",
    "FREITEXT_KLEINWUCHS": "free text: short stature",
    "FREITEXT_DYSPLASIE": "free text: dysplasia",
    "CCD_CMPD": "CCD / CMPD",
    "PSACH_MED": "PSACH / MED",
    "OI": "OI",
}

# Whole sentences from the source files. The screening reasons, the titles of
# the preregistrations and the lesion descriptions are free text and are
# translated here in full -- they are delivered as supplementary tables.
PROSA: dict[str, str] = {
    'kein differenzierter Arm':
        'no differentiated arm',
    'kein naiver Arm':
        'no undifferentiated arm',
    'nur induzierte Chondrozyten; kein naiver Arm':
        'induced chondrocytes only; no undifferentiated arm',
    'nur Tag 7; kein naiver Arm':
        'day 7 only; no undifferentiated arm',
    'siCPM nur Tag 7; kein naiver Arm':
        'siCPM, day 7 only; no undifferentiated arm',
    'shCTR9 nur Tag 7; kein naiver Arm':
        'shCTR9, day 7 only; no undifferentiated arm',
    'COL10A1-KO nur Tag 42; kein naiver Arm':
        'COL10A1-KO, day 42 only; no undifferentiated arm',
    'Knockdown ohne naiven Arm':
        'knockdown without an undifferentiated arm',
    'NAT10-KD ohne naiven Arm':
        'NAT10-KD without an undifferentiated arm',
    'PDK4-siRNA ohne naiven Arm':
        'PDK4-siRNA without an undifferentiated arm',
    'siUHRF1 ohne naiven Arm':
        'siUHRF1 without an undifferentiated arm',
    'siRUNX2 ohne naiven Arm':
        'siRUNX2 without an undifferentiated arm',
    'siUSP34 ohne naiven Arm':
        'siUSP34 without an undifferentiated arm',
    'siCHD7 ohne naiven Arm':
        'siCHD7 without an undifferentiated arm',
    'shVAPB ohne naiven Arm':
        'shVAPB without an undifferentiated arm',
    'shCald1 ohne naiven Arm':
        'shCald1 without an undifferentiated arm',
    'lncRNA-MRF-KD ohne naiven Arm':
        'lncRNA-MRF-KD without an undifferentiated arm',
    'miRNA-Ueberexpression ohne naiven Arm':
        'miRNA overexpression without an undifferentiated arm',
    'Osteoblasten ohne naiven Arm':
        'osteoblasts without an undifferentiated arm',
    'CAVD-Gewebe ohne naiven Arm':
        'CAVD tissue without an undifferentiated arm',
    'Plaquegewebe ohne naiven Arm':
        'plaque tissue without an undifferentiated arm',
    'Klappengewebe ohne naiven Arm':
        'valve tissue without an undifferentiated arm',
    'hMSC-Knockdowns nur differenziert; kein naiver Arm':
        'hMSC knockdowns differentiated only; no undifferentiated arm',
    'Darmorganoide; kein naiver gegen differenzierten Arm':
        'intestinal organoids; no undifferentiated vs differentiated arm',
    'STAT5A-KO ohne osteogenen/chondrogenen Differenzierungsarm':
        'STAT5A-KO without an osteogenic or chondrogenic differentiation arm',
    'LGMN-KD/OE ohne Differenzierungsarm':
        'LGMN knockdown / overexpression without a differentiation arm',
    'Chondroblastom ohne Differenzierungsarm':
        'chondroblastoma without a differentiation arm',
    'GCT-Linien ohne Differenzierungsarm':
        'giant cell tumour lines without a differentiation arm',
    'Chondrosarkomlinien shBCAT1 ohne Differenzierung':
        'chondrosarcoma lines, shBCAT1, without differentiation',
    'nur iPSC; kein Differenzierungsarm':
        'iPSC only; no differentiation arm',
    'KO nur ohne Differenzierung; MSC-Verlauf ohne KO; Zelle leer':
        'knockout without differentiation and MSC time course without knockout; the 2x2 cell is empty',
    'CALD1-siRNA im Basalmedium fehlt; Zelle leer':
        'CALD1-siRNA in basal medium missing; the 2x2 cell is empty',
    'Zelle WT-differenziert fehlt vollstaendig':
        'the wild-type differentiated cell is missing entirely',
    'Klappengewebe ohne Perturbationsarm':
        'valve tissue without a perturbation arm',
    'Spendergewebe ohne Perturbationsarm':
        'donor tissue without a perturbation arm',
    'Bandscheibenzellen ohne Perturbationsarm':
        'intervertebral disc cells without a perturbation arm',
    'myeloide Zellen ohne Perturbationsarm':
        'myeloid cells without a perturbation arm',
    'Reporter-Sortierung ohne Perturbationsarm':
        'reporter sorting without a perturbation arm',
    'Sortierung ohne Perturbationsarm':
        'cell sorting without a perturbation arm',
    'Transplantation ohne Perturbationsarm':
        'transplantation without a perturbation arm',
    'Reprogrammierungsverlauf ohne Perturbationsarm':
        'reprogramming time course without a perturbation arm',
    'hFOB-Differenzierung ohne Perturbationsarm':
        'hFOB differentiation without a perturbation arm',
    'Geruest ohne Perturbationsarm; n = 1':
        'scaffold experiment without a perturbation arm; n = 1',
    'Tumorkohorte ohne Perturbations- und Differenzierungsarm':
        'tumour cohort without a perturbation or differentiation arm',
    'Tag 3 gegen Tag 7 beide bereits induziert':
        'day 3 vs day 7, both already induced',
    'Tag 2 und Tag 7 beide bereits induziert':
        'day 2 and day 7, both already induced',
    'TRPV4; Tag 28 und Tag 56 beide bereits induziert':
        'TRPV4; day 28 and day 56, both already induced',
    'alle Proben undifferenziert':
        'all samples undifferentiated',
    'Perturbation nur im differenzierten Arm; zudem n = 1':
        'perturbation in the differentiated arm only; also n = 1',
    'n = 1 je Gruppe':
        'n = 1 per group',
    'n = 1 je Zelle':
        'n = 1 per 2x2 cell',
    'WT-Linie nur n = 1':
        'wild-type line has n = 1 only',
    'Maus':
        'mouse',
    'Einzelzelldaten':
        'single-cell data',
    'sc/snRNA-seq':
        'sc/snRNA-seq',
    'Einzelzell-Multiom (GEX+ATAC)':
        'single-cell multiome (GEX + ATAC)',
    'Einzelzell-CRISPRi (STING-seq)':
        'single-cell CRISPRi (STING-seq)',
    'Einzelzell-/Einzelkernformat, Muskel':
        'single-cell / single-nucleus format, muscle',
    'MPRA; keine Matrix auf Genebene':
        'MPRA; no gene-level matrix',
    'SubSeries von GSE202147; derselbe Grund':
        'SubSeries of GSE202147; same reason',
    'wie GSE226406; ueberwiegend ChIP':
        'as GSE226406; predominantly ChIP',
    'bereits gesichtet; humaner Teil ohne Perturbation':
        'already screened; the human arm carries no perturbation',
    'bereits gesichtet; nur DE-Tabellen; am 2026-08-18 an der Dateiliste bestaetigt':
        'already screened; differential-expression tables only; confirmed against the file list on 2026-08-18',
    'bereits gesichtet; Gewebequelle ist keine Perturbation':
        'already screened; tissue source is not a perturbation',
    'bereits eingeschlossen (FN1 C123R / C231W)':
        'already included (FN1 C123R / C231W)',
    'bereits eingeschlossen (RB1 +/-)':
        'already included (RB1 +/-)',
    'bereits eingeschlossen (SERPINA3-KD chondrogen)':
        'already included (SERPINA3-KD, chondrogenic)',
    'bereits eingeschlossen (LINC01638-KD)':
        'already included (LINC01638-KD)',
    'bereits eingeschlossen (ARSB MPS VI)':
        'already included (ARSB / MPS VI)',
    'bereits eingeschlossen (RNF4-KD)':
        'already included (RNF4-KD)',
    'bereits eingeschlossen (MIR181A1HG-KD)':
        'already included (MIR181A1HG-KD)',
    'SERPINA3-KD osteogen; Tag 0 gegen Tag 3/7; n 3/3':
        'SERPINA3-KD, osteogenic; day 0 vs day 3/7; n 3/3',
    'ERCC6L2-KD; Expansion gegen Osteogenmedium; n 2/2':
        'ERCC6L2-KD; expansion vs osteogenic medium; n 2/2',
    'YAP/TAZ-siRNA in hFOB; 0 nM gegen 5 nM BMP2; n 3/3':
        'YAP/TAZ-siRNA in hFOB; 0 nM vs 5 nM BMP2; n 3/3',
    'TP53 (Li-Fraumeni) gegen WT; D0 MSC gegen D7/14/17; n 2/2':
        'TP53 (Li-Fraumeni) vs wild type; day 0 MSC vs day 7/14/17; n 2/2',
    'RB1-mut und RB1-del je gegen die genkorrigierte Linie; zwei Punkte; n 3/3':
        'RB1-mut and RB1-del, each against the gene-corrected line; two datasets; n 3/3',
    'keine Diagnose als Laesionsachse':
        'no diagnosis as the lesion axis',
    'iPSC gegen Sklerotom, nur Wildtyp -- keine Diagnose als Laesionsachse':
        'iPSC vs sclerotome, wild type only -- no diagnosis as the lesion axis',
    'Retinaorganoide, Mueller-Glia -- keine skelettale Entitaet':
        'retina organoids, Mueller glia -- not a skeletal entity',
    'Werner-Syndrom als WRN-Knockdown in hESC modelliert -- Engineering-Eingriff, keine Diagnose; Linie statt Patient':
        'Werner syndrome modelled as WRN knockdown in hESC -- an engineered lesion, not a diagnosis; a line, not a patient',
    'dieselbe Studie, ChIP-Arm; WRN-Knockdown in hESC. Der Tag-0/14-Arm des Begleitdatensatzes ist MC3T3-E1':
        'same study, ChIP arm; WRN knockdown in hESC. The day 0/14 arm of the companion dataset is MC3T3-E1',
    'humaner OPLL-Arm ohne Differenzierungsachse; die Tag-0/14-Achse liegt in MC3T3-E1 -- Maus, immortalisiert':
        'human OPLL arm without a differentiation axis; the day 0/14 axis is in MC3T3-E1 -- mouse, immortalised',
    'FOP, iMSC Tag 0/6, FOP gegen resFOP -- vollstaendiges 2x2, aber BEREITS einer der achtzehn Punkte (ACVR1 / FOP)':
        'FOP, iMSC day 0/6, FOP vs resFOP -- a complete 2x2, but ALREADY one of the eighteen datasets (ACVR1 / FOP)',
    'MPS VI, Tag 0/14, 4 Patientenlinien gegen isogene Korrektur -- vollstaendiges 2x2, aber BEREITS einer der achtzehn Punkte (ARSB / MPS VI)':
        'MPS VI, day 0/14, four patient lines vs isogenic correction -- a complete 2x2, but ALREADY one of the eighteen datasets (ARSB / MPS VI)',
    'MPS VI, 4 Patientenlinien mit isogener Korrektur, iPS -> Tag 14 chondrogen':
        'MPS VI, four patient lines with isogenic correction, iPS -> day 14 chondrogenic',
    'Skelettmuskel, keine osteogene/chondrogene Achse':
        'skeletal muscle, no osteogenic or chondrogenic axis',
    'myogene Spezifizierung, keine skelettale Achse':
        'myogenic specification, no skeletal axis',
    'Myotone Dystrophie, iPSC, keine skelettale Achse':
        'myotonic dystrophy, iPSC, no skeletal axis',
    'familiaeres Krebssyndrom, keine skelettale Achse':
        'familial cancer syndrome, no skeletal axis',
    'FOP-Monozyten, entzuendliche Signatur, keine Differenzierungsachse':
        'FOP monocytes, inflammatory signature, no differentiation axis',
    'iPSC-abgeleitete kraniofaziale MSC, keine Laesion und keine Kontrollgruppe im 2x2-Sinn':
        'iPSC-derived craniofacial MSC, no lesion and no control group in the 2x2 sense',
    'Somitenuhr, Entwicklungsmodell ohne Laesion':
        'somite clock, a developmental model without a lesion',
    'Implantate aus stabilem/hypertrophem Knorpel, keine Laesion':
        'implants of stable or hypertrophic cartilage, no lesion',
    'Dermis- und Fettgewebs-Stromazellen, gesunde Spender, keine Laesion':
        'dermal and adipose stromal cells, healthy donors, no lesion',
    'ZNF145-Ueberexpression in gesunden MSC, kein Patientendefekt, kein 2x2':
        'ZNF145 overexpression in healthy MSC, no patient defect, no 2x2',
    'RB1-Patientenlinie gegen WT, MSC, Osteogenese Tag 0/7/14/21':
        'RB1 patient line vs wild type, MSC, osteogenesis day 0/7/14/21',
    'Rhabdomyosarkom, Tumorentitaet':
        'rhabdomyosarcoma, a tumour entity',
    'FSHD2-Myotuben, Muskel; zudem Einzelkernformat (A9)':
        'FSHD2 myotubes, muscle; also single-nucleus format (A9)',
    'FSHD-Muskelfasern':
        'FSHD muscle fibres',
    'DMD, myogene Kulturen':
        'DMD, myogenic cultures',
    'DMD, Muskel':
        'DMD, muscle',
    'LGMDR21, Skelettmuskel':
        'LGMDR21, skeletal muscle',
    'LAMA2, Muskelstammzellen':
        'LAMA2, muscle stem cells',
    'Mesoangioblasten, Muskel':
        'mesoangioblasts, muscle',
    'vorab benannt; 4 Zellen (Line #1-#4), chondrogen':
        'named in advance; four cells (lines 1-4), chondrogenic',
    'vorab benannt; FOP/resFOP = EINE Linie -> 1 Zelle, ex1-3 sind Replikate':
        'named in advance; FOP/resFOP are ONE line -> one cell, ex1-3 are replicates',
    'vorab benannt; 3 Spender, chondrogen; E2 nicht erfuellt (siRNA) -> nur dWT-tragend':
        'named in advance; three donors, chondrogenic; criterion E2 not met (siRNA) -> carries dWT only',
    'vorab benannt; 3 Spender, osteogen; E2 nicht erfuellt (siRNA) -> nur dWT-tragend':
        'named in advance; three donors, osteogenic; criterion E2 not met (siRNA) -> carries dWT only',
    'eigene Daten; WT1-3/KO9/46/75 sind Klone EINER Linie -> 1 Spender, 2 Achsen':
        'own data; WT1-3 and KO9/46/75 are clones of ONE line -> one donor, two axes',
    'kein naiver Arm: Tag 34 gegen Tag 44, beide differenzierter Knorpel; zudem Klone einer Linie':
        'no undifferentiated arm: day 34 vs day 44, both differentiated cartilage; also clones of one line',
    "kein naiver Arm: nur differenzierte Knorpelpellets; 'biologische Replikate' sind keine Spender":
        'no undifferentiated arm: differentiated cartilage pellets only; biological replicates are not donors',
    'nicht zerlegbar':
        'not decomposable',
    'Entkopplung, anderer Fall (Linie erreicht)':
        'decoupling, the other case (lineage reached)',
    'keine verknuepfte Publikation in GEO':
        'no linked publication in GEO',
    'ok':
        'ok',
    'Modul voll (173 Gene)':
        'full module (173 genes)',
    'Modul ohne die 5 % groessten |Delta| (164 Gene)':
        'module without the top 5 % |Delta| (164 genes)',
    'Modul ohne die 10 % groessten |Delta| (155 Gene)':
        'module without the top 10 % |Delta| (155 genes)',
    'Modul ohne die 20 % groessten |Delta| (138 Gene)':
        'module without the top 20 % |Delta| (138 genes)',
    'Modul ohne die 30 % groessten |Delta| (121 Gene)':
        'module without the top 30 % |Delta| (121 genes)',
    'Modul (173 Gene) [ohne HyperChon]':
        'module (173 genes), hypertrophic zone excluded',
    'Modul (173 Gene) [ohne Hyper+Praehyper]':
        'module (173 genes), hypertrophic and prehypertrophic zones excluded',
    'PA309 gegen Modul [ohne HyperChon]':
        'disease genes vs programme, hypertrophic zone excluded',
    'PA309 gegen Modul [ohne Hyper+Praehyper]':
        'disease genes vs programme, hypertrophic and prehypertrophic zones excluded',
    'Positivkontrolle chondrogen-naiv [ohne HyperChon]':
        'positive control: chondrogenic vs undifferentiated, hypertrophic zone excluded',
    'Positivkontrolle chondrogen-naiv [ohne Hyper+Praehyper]':
        'positive control: chondrogenic vs undifferentiated, hypertrophic and prehypertrophic zones excluded',
    'ARSB-Patientenmutation, isogen korrigiert':
        'ARSB patient mutation, isogenically corrected',
    'RB1 +/- (Retinoblastom-Patientenlinie)':
        'RB1 +/- (retinoblastoma patient line)',
    'SERPINA3-siRNA (Engineering, kein Patientendefekt)':
        'SERPINA3-siRNA (engineered, not a patient defect)',
    'LAMA5-KO in gesunder Linie (Engineering)':
        'LAMA5 knockout in a healthy line (engineered)',
    'ACVR1 R206H, isogen korrigiert (resFOP)':
        'ACVR1 R206H, isogenically corrected (resFOP)',
    'Punkt der achtzehn':
        'one of the eighteen datasets',
    'zwei Punkte der achtzehn':
        'two of the eighteen datasets',
    'Abb. S3C, studienweise':
        'Figure S3C, study-wise',
    'PRAEREG_F2F – die Dreifachzerlegung als Hypothese':
        'PREREG_F2F -- the three-way decomposition as a hypothesis',
    'Vorregistrierung der Gesamtstudie – Typologie des Matrixproduktionsversagens':
        'Preregistration of the whole study -- a typology of matrix-production failure',
    'Vorregistrierung FOLLOWUP – Donor-resolved Differenzierungsatlas und unabhaengige isogene 2x2-Läsionsstudie':
        'Preregistration FOLLOW-UP -- donor-resolved differentiation atlas and an independent isogenic 2x2 lesion study',
    'Vorregistrierung M-A – humangenetischer Anker, mit voller Kraft':
        'Preregistration M-A -- the human-genetics anchor, at full power',
    'Vorregistrierung M-B – Individualität patientenweise messen':
        'Preregistration M-B -- measuring individuality patient by patient',
    'Vorregistrierung M-D – die donoraufgelöste Schere':
        'Preregistration M-D -- the donor-resolved scissors',
    'Vorregistrierung M-E – ist die Schere ein Befund oder ein Größenunterschied?':
        'Preregistration M-E -- are the scissors a finding or a difference in sample size?',
    'Vorregistrierung S1 – induktionsangeglichene Null für `FN1`':
        'Preregistration S1 -- an induction-matched null for FN1',
    'Vorregistrierung S2 – Positivkontrolle Plasmazelldifferenzierung':
        'Preregistration S2 -- plasma-cell differentiation as a positive control',
    'Vorregistrierung S4 – Die distale Achse und der sekretorische TF-Arm':
        'Preregistration S4 -- the distal axis and the secretory transcription-factor arm',
    'Vorregistrierung S5 – Konvergenzachsen':
        'Preregistration S5 -- convergence axes',
    'Vorregistrierung S6 – Entdeckung und Validierung, getrennte Hälften':
        'Preregistration S6 -- discovery and validation in separate halves',
    'Vorregistrierung S7 – Zelltyp-Frage und 19. Punkt':
        'Preregistration S7 -- the cell-type question and a possible nineteenth dataset',
    'Vorregistrierung S8 – orthogonale Einzelzellpruefung der Papierclaims':
        'Preregistration S8 -- an orthogonal single-cell test of the published claims',
    'Vorregistrierung S9 – neuer R2-Kandidat fuer zelltypspezifische Laesionsantwort':
        'Preregistration S9 -- a new R2 candidate for a cell-type-specific lesion response',
    'Vorregistrierung S10 – R2 mit gepinnter HCA-Primärreferenz (C8-Kollaps als Ausschluss)':
        'Preregistration S10 -- R2 with a pinned HCA primary reference (C8 collapse as an exclusion)',
    'Vorregistrierung S12 – orthogonale Datensatz-Triangulation':
        'Preregistration S12 -- orthogonal dataset triangulation',
    'Protokoll M-A – der humangenetische Anker, mit voller Kraft':
        'Protocol M-A -- the human-genetics anchor, at full power',
    'Protokoll M-B – Individualität patientenweise gemessen':
        'Protocol M-B -- individuality measured patient by patient',
    'Protokoll M-C – die Schere an echten Diagnosen':
        'Protocol M-C -- the scissors tested on real diagnoses',
    'Protokoll M-D – die Schere donoraufgelöst gesucht':
        'Protocol M-D -- the scissors sought donor by donor',
    'Protokoll M-E – die Schere ist ein Größenunterschied':
        'Protocol M-E -- the scissors are a difference in sample size',
    'Protokoll M-GESAMT – das Entscheidungstor':
        'Protocol M-OVERALL -- the decision gate',
    'Protokoll M-GESAMT v2 – das Entscheidungstor nach Phase D':
        'Protocol M-OVERALL v2 -- the decision gate after phase D',
    'Protokoll – orthogonale Messebenen (Chromatin, Methylom)':
        'Protocol -- orthogonal measurement levels (chromatin, methylome)',
    'S7 Schritt A – Vorregistrierung vor der Rechnung':
        'S7 step A -- preregistration before the computation',
    'S7 Nachtrag 1 – Freigabe der neuen scRNA-GEO-Sichtung':
        'S7 addendum 1 -- release of the new scRNA GEO screen',
    'S7 Nachtrag 2 – Kontextauswertung `GSE337700`':
        'S7 addendum 2 -- context analysis of GSE337700',
    'S8 Schritt A – Freigabe der orthogonalen Einzelzellpruefung':
        'S8 step A -- release of the orthogonal single-cell test',
    'S8 Schritt E – Bericht und Entscheidung':
        'S8 step E -- report and decision',
    'S9 Schritt A – Freigabe des R2-Impact-Scans':
        'S9 step A -- release of the R2 impact scan',
    'S9 Schritt E – Bericht und Halt nach S9-AB3':
        'S9 step E -- report and halt after S9-AB3',
    'S10 Schritt A – Freigabe der R2-Rechnung mit gepinnter Referenz':
        'S10 step A -- release of the R2 computation with a pinned reference',
}

# Rule-based clean-up: whatever is left as an abbreviation or as a residual
# German word inside a display string.
MUSTER: list[tuple[str, str]] = [
    (r"\bFREITEXT_", "free text: "),
    (r"\bLinie (\d)", r"line \1"),
    (r"\bLinie\b", "line"),
    (r"\bZielachse\b", "target axis"),
    (r"\bLinienkontrolle\b", "lineage control"),
    (r"\bHint\b", "background"),
    (r"\bbasi\b", "baseline-stratified"),
    (r"\bnaiv(e|er|en|em)?\b", "undifferentiated"),
    (r"\bnaive\b", "undifferentiated"),
]

# Columns that are pure machine keys and are NOT translated as values --
# otherwise the join between panel CSV and number check breaks.
NICHT_UEBERSETZEN = {"gse", "gen", "punkt", "code", "lauf", "matching"}


_DEUTSCH = re.compile(
    r"[äöüÄÖÜß]|"
    r"(kein|keine|ohne|nicht|nur|gegen|und|der|die|das|Arm|Zelle|Achse|"
    r"Spender|Probe|Grenze|Anteil|Datensatz|Fenster|Zerlegung|Modul|"
    r"Vorregistrierung|Protokoll|Laesion|Entitaet|bereits|schon|beide|"
    r"vorab|eigene|Punkt|Freigabe|Bericht|Schritt|Nachtrag)")


def _wert(v: object) -> object:
    if not isinstance(v, str):
        return v
    s = WERTE.get(v)
    if s is not None:
        return s
    # The preregistration titles carry an em dash (U+2014) in some places and
    # an en dash (U+2013) in others. Normalise them for the lookup key.
    s = PROSA.get(v) or PROSA.get(v.replace("—", "–"))
    if s is not None:
        return s
    # The pattern rules apply to text that is already English. On German
    # prose they would produce half-sentences ("kein undifferentiated Arm"),
    # so such strings are left standing and are caught by rest_deutsch().
    if _DEUTSCH.search(v):
        return v
    for muster, ersatz in MUSTER:
        v = re.sub(muster, ersatz, v)
    return v


def _ist_text(reihe: pd.Series) -> bool:
    """A text column? pandas 3 gives dtype 'str', pandas 2 gives 'object'."""
    return pd.api.types.is_string_dtype(reihe) or reihe.dtype == object


# ----------------------------------------------------------------- columns
# The column names of the panel and table files. Internally the pipeline
# computes with short names; what is delivered is English, and the mapping
# lives at this one place, exactly as it does for the display values. `None`
# means the column is dropped on write, because it was only a working note.
SPALTEN: dict[str, str | None] = {
    'OR': 'odds_ratio',
    'OR_gematcht': 'odds_ratio_matched',
    'OR_mde80': 'odds_ratio_detection_limit',
    'OR_roh': 'odds_ratio_raw',
    'S1_beobachtet': 'S1_observed',
    'S1_mde80': 'S1_detection_limit',
    'S1_null_mittel': 'S1_null_mean',
    'S2_beobachtet': 'S2_observed',
    'S2_mde80': 'S2_detection_limit',
    'S2_null_mittel': 'S2_null_mean',
    'S3a_beobachtet': 'S3a_observed',
    'S3a_mde80': 'S3a_detection_limit',
    'S3a_null_mittel': 'S3a_null_mean',
    'S3b_beobachtet': 'S3b_observed',
    'S3b_mde80': 'S3b_detection_limit',
    'S3b_null_mittel': 'S3b_null_mean',
    'U': 'mann_whitney_u',
    'above_mde80': 'above_detection_limit',
    'abrufdatum': 'retrieval_date',
    'achse': 'axis',
    'achse_lang': 'axis_label',
    'achse_name': 'axis_name',
    'amplitude_beobachtet': 'amplitude_observed',
    'amplitude_null_mittel': 'amplitude_null_mean',
    'ander_kontrast': 'other_axis_contrast',
    'ander_status': 'other_axis_status',
    'ander_z': 'other_axis_z',
    'anteil_abs_z_ueber_2': 'fraction_abs_z_above_2',
    'anteil_bestaetigt': 'fraction_confirmed',
    'anteil_entfernt': 'fraction_removed',
    'anteil_z_ueber_2': 'fraction_z_above_2',
    'art': 'type',
    'artefakt': 'artefact',
    'befund': 'finding',
    'begruendung': 'reason',
    'beispiel': 'example',
    'beleg': 'evidence',
    'beobachtet': 'observed',
    'bestanden': 'passed',
    'bestandteile': 'set_members',
    'datei': 'file',
    'datensatz': 'dataset',
    'delta_beobachtet': 'delta_observed',
    'delta_hyper_minus_mes': 'delta_hypertrophic_minus_mesenchymal',
    'diagnose': 'diagnosis',
    'diff_arm': 'differentiated_arm',
    'dwt_modulwert': 'dwt_module_value',
    'e2': 'criterion_e2',
    'ebene': 'level',
    'effekt': 'effect',
    'eichbar': 'calibratable',
    'eichstatus': 'calibration_status',
    'eichung': 'calibration',
    'eichung_bestanden': 'calibration_passed',
    'eichung_bestanden_abgelegt': 'calibration_passed_filed',
    'eichung_kontrast': 'calibration_contrast',
    'eichung_mde80': 'calibration_detection_limit',
    'eichung_status_abgelegt': 'calibration_status_filed',
    'eichung_z': 'calibration_z',
    'eichung_z_abgelegt': 'calibration_z_filed',
    'eichungstyp': 'calibration_type',
    'eigen_kontrast': 'own_axis_contrast',
    'eigen_status': 'own_axis_status',
    'eigen_z': 'own_axis_z',
    'eingriff': 'intervention',
    'einheit': 'unit',
    'entitaet': 'entity',
    'erstautor': 'first_author',
    'band': 'volume',
    'est_max': 'estimate_max',
    'est_min': 'estimate_min',
    'faktor': 'factor',
    'fall': 'case',
    'fassung': 'version',
    'fenster': 'window',
    'geeicht': 'calibrated',
    'gegen': 'compared_to',
    'gen': 'gene',
    'genkarte': 'gene_map',
    'geo_titel': 'geo_title',
    'gewebe': 'tissue',
    'go': 'go_term',
    'grenze': 'detection_limit',
    'grenze_max': 'detection_limit_max',
    'grenze_min': 'detection_limit_min',
    'groesse': 'quantity',
    'groesse_seite': 'set_size_per_side',
    'gse': 'accession',
    'hat_grenze': 'has_detection_limit',
    'herkunft': 'origin',
    'hi': 'ci_high',
    'hypothese': 'hypothesis',
    'i_naiv_faellt': 'step_i_leaves_undifferentiated',
    'ii_linie_erreicht': 'step_ii_lineage_reached',
    'iii_modul_ueber_grenze': 'step_iii_module_above_limit',
    'jahr': 'year',
    'kategorie': 'category',
    'kategorie_de': None,
    'kennzahl': 'statistic_id',
    'ki_schliesst_1_aus': 'ci_excludes_1',
    'klasse': 'class',
    'kohorte': 'cohort',
    'kohorte_lang': 'cohort_label',
    'kompartiment': 'compartment',
    'konkordant': 'concordant',
    'konkordanz': 'concordance',
    'konkordanz_mde80': 'concordance_detection_limit',
    'konkordanz_null': 'concordance_null',
    'konkordanz_null_sd': 'concordance_null_sd',
    'konkordanz_p': 'concordance_p',
    'konkordanz_z': 'concordance_z',
    'kontrast': 'contrast',
    'kontrast_hi': 'contrast_ci_high',
    'kontrast_invitro': 'contrast_in_vitro_reachable',
    'kontrast_lo': 'contrast_ci_low',
    'kontrast_median': 'contrast_median',
    'kontrast_modul': 'contrast_module',
    'kontrast_positivkontrolle': 'contrast_positive_control',
    'kontrast_vorreg': 'contrast_preregistered',
    'kurve': 'curve',
    'laesion': 'lesion',
    'lauf': 'run',
    'lo': 'ci_low',
    'mde80': 'detection_limit',
    'mde80_delta': 'detection_limit_delta',
    'mde80_median': 'detection_limit_median',
    'mde80_modul': 'detection_limit_module',
    'mde80_positivkontrolle': 'detection_limit_positive_control',
    'mde80_rho': 'detection_limit_rho',
    'med': 'median',
    'median_bestanden': 'median_passed',
    'median_bi': 'median_biallelic',
    'median_d_basis': 'median_baseline_shift',
    'median_durchgefallen': 'median_failed',
    'median_mono': 'median_monoallelic',
    'median_pub': 'median_publications',
    'median_pub_hintergrund': 'median_publications_background',
    'mittel': 'mean',
    'mittlerer_z': 'mean_z',
    'modul_konkordanz': 'module_concordance',
    'modul_mde80': 'module_detection_limit',
    'modul_n': 'module_n',
    'modul_z': 'module_z',
    'modus': 'mode_of_inheritance',
    'n1': 'n_group_1',
    'n2': 'n_group_2',
    'n_a': 'n_group_a',
    'n_anderer_fall': 'n_other_case',
    'n_b': 'n_group_b',
    'n_bestaetigt': 'n_confirmed',
    'n_bestanden': 'n_passed',
    'n_bi': 'n_biallelic',
    'n_datensaetze': 'n_datasets',
    'n_datensaetze_chondro': 'n_datasets_chondrogenic',
    'n_datensaetze_osteo': 'n_datasets_osteogenic',
    'n_durchgefallen': 'n_failed',
    'n_entfernt': 'n_removed',
    'n_gene': 'n_genes',
    'n_gene_messbar': 'n_genes_measurable',
    'n_hintergrund': 'n_background',
    'n_invitro_a': 'n_in_vitro_reachable',
    'n_klasse_genom': 'n_class_genome',
    'n_klasse_panel': 'n_class_panel',
    'n_kontrollen': 'n_controls',
    'n_modul': 'n_module_genes',
    'n_modul_ueber_grenze': 'n_module_above_limit',
    'n_module': 'n_module_genes',
    'n_modulgene': 'n_module_genes',
    'n_modulgene_gemeinsam': 'n_module_genes_shared',
    'n_mono': 'n_monoallelic',
    'n_paare': 'n_pairs',
    'n_panel_in_hg': 'n_panel_in_background',
    'n_patienten': 'n_patients',
    'n_proben': 'n_samples',
    'n_proben_gds': 'n_samples_gds',
    'n_punkte': 'n_datasets',
    'n_satz': 'n_gene_set',
    'n_satz_gesamt': 'n_gene_set_total',
    'n_satz_im_hg': 'n_gene_set_in_background',
    'n_schnitt': 'n_intersection',
    'n_spezimen': 'n_specimens',
    'n_studien': 'n_studies',
    'n_vorreg_a': 'n_preregistered',
    'n_widerlegt': 'n_refuted',
    'n_wiederholungen': 'n_repeats',
    'n_zellen': 'n_cells',
    'n_zellen_min': 'n_cells_min',
    'n_zerlegbar_i': 'n_decomposable',
    'naiv_kontrast': 'undifferentiated_contrast',
    'naiv_marker_achsendiff': 'undifferentiated_axis_difference',
    'naiv_marker_eigen': 'undifferentiated_own_markers',
    'naiv_marker_fremd': 'undifferentiated_other_markers',
    'naiv_mde80': 'undifferentiated_detection_limit',
    'naiv_modulwert_roh': 'undifferentiated_module_raw',
    'naiv_modulwert_z': 'undifferentiated_module_z',
    'naiv_status': 'undifferentiated_status',
    'naiv_z': 'undifferentiated_z',
    'naiver_arm': 'undifferentiated_arm',
    'null': 'null_model',
    'null_intakt': 'null_intact',
    'null_mittel': 'null_mean',
    'nzieh': 'n_draws',
    'nziehungen': 'n_draws',
    'ohne': 'left_out',
    'p': 'p_value',
    'p_emp': 'p_empirical',
    'p_invitro': 'p_in_vitro_reachable',
    'p_modul': 'p_module',
    'p_perm_oben': 'p_permutation_upper',
    'p_perm_two_sided': 'p_permutation_two_sided',
    'p_perm_upper': 'p_permutation_upper',
    'p_positivkontrolle': 'p_positive_control',
    'p_roh_fisher': 'p_raw_fisher',
    'p_vorreg': 'p_preregistered',
    'pool': 'genes_pooled',
    'probe': 'sample',
    'punkt': 'dataset_id',
    'quelle': 'source',
    'rang': 'rank',
    'rang_abs_delta': 'rank_abs_delta',
    'rang_null': 'rank_null',
    'rang_null_sd': 'rank_null_sd',
    'rang_p': 'rank_p',
    'rang_z': 'rank_z',
    'rausch_mittel': 'noise_mean',
    'rausch_q05': 'noise_q05',
    'rausch_q95': 'noise_q95',
    'rechenbar': 'computable',
    'repliziert': 'replicated',
    'rho_modulgene': 'rho_module_genes',
    'ri': 'direction_ri',
    'richtung': 'direction',
    'richtung_beobachtet': 'direction_observed',
    'richtung_null_mittel': 'direction_null_mean',
    'richtung_null_sd': 'direction_null_sd',
    'richtung_p': 'direction_p',
    'richtung_z': 'direction_z',
    'rolle': 'role',
    'runde': 'round',
    'satz': 'gene_set',
    'satz_a': 'marker_set',
    'satz_ander': 'other_marker_set',
    'satz_eigen': 'own_marker_set',
    'satzlab': 'gene_set_label',
    'schwelle': 'threshold',
    'schwelle_erreicht': 'threshold_reached',
    'seite': 'side',
    'seiten': 'pages',
    'share_match': 'directional_share',
    'spender': 'donor',
    'spezimen': 'specimen',
    'stadium': 'stage',
    'statistik': 'statistic',
    'statistik_z': 'statistic_z',
    'studie': 'study',
    'stufe': 'step',
    'symbol': 'gene_symbol',
    'teil': 'part',
    'teilmenge': 'subset',
    'titel': 'title',
    'treffer': 'hits',
    'ueber_bonferroni': 'above_bonferroni',
    'ueber_grenze': 'above_detection_limit',
    'ueber_mde80': 'above_detection_limit',
    'urteil': 'verdict',
    'urteil_invitro': 'verdict_in_vitro_reachable',
    'urteil_vorreg': 'verdict_preregistered',
    'urteilklasse': 'verdict_class',
    'variante': 'variant',
    'wert': 'value',
    'wilson_hi': 'wilson_ci_high',
    'wilson_lo': 'wilson_ci_low',
    'z_invitro': 'z_in_vitro_reachable',
    'z_korr': 'z_corrected',
    'z_modul': 'z_module',
    'z_positivkontrolle': 'z_positive_control',
    'z_ueber_2': 'fraction_z_above_2',
    'z_vorreg': 'z_preregistered',
    'zeitschrift': 'journal',
    'zelle': 'cell',
    'zerlegbar': 'decomposable',
    'zielgroesse': 'target_quantity',
    'zonenrang': 'zone_rank',
}


def spalten_englisch(df: pd.DataFrame) -> pd.DataFrame:
    """Translate the column headers and drop the working columns."""
    weg = [k for k in df.columns if SPALTEN.get(k, "") is None]
    d = df.drop(columns=weg) if weg else df
    return d.rename(columns={k: v for k, v in SPALTEN.items()
                             if v is not None and k in d.columns})


def englisch(df: pd.DataFrame) -> pd.DataFrame:
    """Put every display value of a panel data frame into English."""
    d = df.copy()
    for k in d.columns:
        if k in NICHT_UEBERSETZEN or not _ist_text(d[k]):
            continue
        # Joined diagnosis axes ("A|B") are resolved token by token FIRST;
        # after that no pattern can cut halfway into an identifier.
        if k == "achse" and d[k].astype(str).str.contains("|", regex=False).any():
            d[k] = d[k].map(diagnoseachse)
        d[k] = d[k].map(_wert).astype(object)
        # The two large name lists act on their own columns only, so that a
        # coincidentally identical free text is not translated with them.
        if k in ("datensatz", "dataset", "einheit"):
            d[k] = d[k].map(lambda v: DATENSATZ.get(v, v))
        if k in ("zelle", "ohne", "spender"):
            d[k] = d[k].map(lambda v: ZELLE.get(v, v))
        if k in ("quelle", "source"):
            d[k] = d[k].map(lambda v: QUELLE.get(v, v))
        if k in ("satz", "gene_set", "satzlab", "panel", "class", "klasse"):
            d[k] = d[k].map(lambda v: SATZ.get(v, v))
    # The accession column is not translated, but it is resolved.
    for k in ("gse", "accession", "studie"):
        if k in d.columns and _ist_text(d[k]):
            d[k] = d[k].map(lambda v: AKZESSION.get(v, v))
    return spalten_englisch(d)


def rest_deutsch(df: pd.DataFrame) -> list[str]:
    """What still looks German after translation -- for the check."""
    verdacht = re.compile(
        r"[äöüÄÖÜß]|"
        r"\b(naiv|Eichung|Zelle|Probe|Grenze|Anteil|Datensatz|Spender|Achse|"
        r"Fenster|Zerlegung|keine|kein|ohne|nicht|Linie|gegen|und|der|die|das)\b")
    treffer: set[str] = set()
    for k in df.columns:
        if not _ist_text(df[k]):
            continue
        for v in df[k].dropna().unique():
            if isinstance(v, str) and verdacht.search(v):
                treffer.add(f"{k}: {v}")
    return sorted(treffer)

# -*- coding: utf-8 -*-
"""
11_cohort_prescreen.py -- mechanical prescreening of the search hits from 50_.

Applies the exclusion codes preregistered in `PRAEREG_M_B.md` §2, as far as
they can be decided from the GEO metadata alone. Everything that survives
this stage goes to MANUAL SCREENING -- there each accession number is
checked individually, not adopted from memory.

The prescreening is deliberately generous: in doubt a series stays in.

Output:
  derived_data/M_patienten/sichtung_mechanisch.csv   all 1424 series with verdict
  derived_data/M_patienten/kandidaten.csv            the survivors
"""
from __future__ import annotations

import pathlib
import re
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
from _module import ERGEBNISSE  # noqa: E402

AUS = ERGEBNISSE / "M_patienten"

# Minimum sample count: 5 patients (E2) + 2 controls (E3).
MIN_PROBEN = 7

# Entities naming a skeletal disease (E1/M4).
SKELETT = re.compile(
    r"osteogenesis imperfecta|\bOI\b|mucopolysaccharid|\bMPS\b|Maroteaux|Hurler|Hunter|"
    r"Morquio|Sanfilippo|pseudoachondroplasia|epiphyseal dysplasia|achondroplasia|"
    r"hypochondroplasia|thanatophoric|dyschondrosteosis|Leri-Weill|SHOX|"
    r"cleidocranial|campomelic|fibrodysplasia ossificans|\bFOP\b|"
    r"skeletal dysplasia|chondrodysplasia|osteochondrodysplasia|short stature|dwarfism|"
    r"growth plate|osteoporosis|osteopetrosis|osteoarthritis|chondrocyte|cartilage|"
    r"osteoblast|bone marrow stromal|mesenchymal stem|MSC\b|fibroblast|iPSC|"
    r"skeletal|bone\b|osteogenic|chondrogenic", re.I)

# Formats without a gene-level matrix (A3) or wrong measurement layer.
KEINE_EXPRESSION = re.compile(
    r"Genome binding|Methylation|Genome variation|SNP|Non-coding RNA profiling by array$|"
    r"Protein profiling|Expression profiling by RT-PCR|Third-party|Genome tiling", re.I)

# Immortalized lines (M2, hard).
LINIEN = re.compile(
    r"\bHEK ?293|\bHeLa\b|\bU2OS\b|\bSaos-?2\b|\bMG-?63\b|\bMC3T3|\bATDC5\b|\bC2C12\b|"
    r"\bHT-?1080\b|\bK562\b|\bA549\b|\bHCT ?116\b|\bMCF-?7\b|\bSW1353\b|\bTC28|"
    r"\bhTERT\b|immortali[sz]ed|cell line model|\bLCL\b|lymphoblastoid", re.I)

# Perturbation formats without patient assignment (M3).
PERTURBATION = re.compile(
    r"knock-?down|knock-?out|siRNA|shRNA|CRISPR|overexpress|transfect|"
    r"\bKD\b|\bKO\b|treatment of .* cells with|stimulat", re.I)

NICHT_HUMAN = re.compile(r"Mus musculus|Rattus|Danio|Gallus|Bos taurus|Sus scrofa", re.I)


def urteile(r: pd.Series) -> tuple[str, str, str]:
    """-> (urteil, code, begruendung)"""
    txt = f"{r.titel} {r.zusammenfassung}"
    n = pd.to_numeric(r.n_proben, errors="coerce")
    typ = str(r.typ)

    if NICHT_HUMAN.search(str(r.taxon)) and "Homo sapiens" not in str(r.taxon):
        return "AUS", "A7", "nicht human"
    if not re.search(r"Expression profiling", typ, re.I) or KEINE_EXPRESSION.fullmatch(typ):
        if not re.search(r"Expression profiling by (array|high throughput sequencing)", typ, re.I):
            return "AUS", "A3", f"keine Genebene-Expressionsmatrix ({typ[:60]})"
    if pd.isna(n) or n < MIN_PROBEN:
        return "AUS", "M1", f"n_proben = {n} < {MIN_PROBEN} (E2/E3 nicht erfuellbar)"
    if not SKELETT.search(txt):
        return "AUS", "M4", "keine skelettale Entitaet im Titel/Abstract"
    if LINIEN.search(txt):
        return "AUS", "M2", "immortalisierte Linie benannt"
    if re.search(r"single[- ]cell|scRNA|snRNA|10x Genomics|Chromium", txt, re.I):
        return "AUS", "A9", "Einzelzellformat -- Pseudobulk je Patient nicht aus Metadaten belegbar"
    if PERTURBATION.search(txt) and not re.search(
            r"patient|donor|subject|biops|cohort|individual", txt, re.I):
        return "AUS", "M3", "Perturbationsformat ohne Patientenzuordnung"
    return "HAND", "-", "Handsichtung"


def main() -> None:
    T = pd.read_csv(AUS / "treffer_roh.csv")
    urt = T.apply(urteile, axis=1, result_type="expand")
    T[["urteil", "code", "begruendung"]] = urt
    T.to_csv(AUS / "sichtung_mechanisch.csv", index=False)

    print("=" * 78)
    print("Mechanical prescreening -- %d series" % len(T))
    print("=" * 78)
    print(T.code.value_counts().to_string())
    K = T[T.urteil == "HAND"].sort_values(
        "n_proben", key=lambda s: -pd.to_numeric(s, errors="coerce"))
    K.to_csv(AUS / "kandidaten.csv", index=False)
    print("\n-> %d series into manual screening" % len(K))


if __name__ == "__main__":
    main()

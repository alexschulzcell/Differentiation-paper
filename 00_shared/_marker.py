# -*- coding: utf-8 -*-
"""
_marker.py -- canonical marker sets for the positive controls of the
orthogonal layers.

These sets are textbook knowledge, not results of this project. They serve
exclusively to CALIBRATE a measurement layer: if a layer does not show the
separation osteogenic vs adipogenic (or osteogenic vs myogenic) at these
genes, its null finding on the 173-gene module says nothing about the
module.

The sets are fixed before any analysis and are not readjusted.
Disjointness has been verified.
"""

# Osteoblastic differentiation: transcription factors, matrix proteins,
# mineralization machinery.
OSTEOGEN = [
    "RUNX2", "SP7", "DLX5", "MSX2", "SATB2", "ATF4",
    "ALPL", "BGLAP", "IBSP", "SPP1", "COL1A1", "COL1A2",
    "POSTN", "MEPE", "PHEX", "DMP1", "SOST", "PTH1R",
]

# Adipogenic differentiation: nuclear regulators and lipid machinery.
ADIPOGEN = [
    "PPARG", "CEBPA", "CEBPB", "CEBPD", "SREBF1", "NR1H3",
    "FABP4", "ADIPOQ", "LEP", "LPL", "PLIN1", "PLIN4",
    "CIDEC", "CD36", "GPD1", "AQP7", "LIPE", "PNPLA2",
]

# Myogenic differentiation.
MYOGEN = [
    "MYOD1", "MYOG", "MYF5", "MYF6", "PAX7", "DES",
    "MYH3", "MYH8", "TNNT1", "TNNT2", "ACTA1", "CKM",
    "MYL4", "TTN", "DMD", "CAV3",
]

# Chondrogenic differentiation.
CHONDROGEN = [
    "SOX9", "SOX5", "SOX6", "ACAN", "COL2A1", "COL9A1",
    "COL11A1", "COMP", "HAPLN1", "MATN3", "PRG4", "WWP2",
]

# Stem-cell/naive markers (should go DOWN along EVERY axis).
NAIV = [
    "THY1", "ENG", "NT5E", "NGFR", "LIF", "KITLG",
    "CXCL12", "MKI67", "CCNB1", "TOP2A",
]

_alle = OSTEOGEN + ADIPOGEN + MYOGEN + CHONDROGEN + NAIV
assert len(_alle) == len(set(_alle)), "marker sets must be disjoint."

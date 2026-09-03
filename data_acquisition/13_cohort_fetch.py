# -*- coding: utf-8 -*-
"""
13_cohort_fetch.py -- download of the seven phase-B cohorts.

Runs AFTER screening (50_, 50b_, 50c_) and AFTER `PRAEREG_M_B.md`. Per
series, loads the series matrix (arrays: the normalized value deposited by
the authors, §3) and -- for the sequencing cohorts -- the supplementary
file with the gene-level matrix.

Storage: data_raw/<GSE>/
"""
from __future__ import annotations

import pathlib
import sys
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
from _module import DATEN  # noqa: E402

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

# GSE -> supplementary files to download (empty = series matrix only)
KOHORTEN = {
    "GSE186141": ["GSE186141_FPKM9.6Col1.vs.2Ctrl.xlsx"],
    "GSE292600": ["GSE292600_raw_counts.txt.gz"],
    "GSE160207": ["GSE160207_EE_OI_RNAseq_counts.txt.gz"],
    "GSE228522": ["GSE228522_matrix_deseq2rlogv1.txt.gz"],
    "GSE77758": [],
    "GSE22855": [],
    "GSE58435": [],
}


def praefix(gse: str) -> str:
    return "GSE" + (gse[3:-3] or "") + "nnn"


def hole(url: str, ziel: pathlib.Path) -> bool:
    if ziel.exists() and ziel.stat().st_size > 1000:
        print(f"   exists:    {ziel.name}")
        return True
    try:
        with urllib.request.urlopen(url, timeout=600) as r, open(ziel, "wb") as f:
            f.write(r.read())
        print(f"   fetched:   {ziel.name}  ({ziel.stat().st_size/1e6:.1f} MB)")
        return True
    except Exception as e:                       # noqa: BLE001
        print(f"   MISSING:   {ziel.name}  ({e})")
        if ziel.exists():
            ziel.unlink()
        return False


def main() -> None:
    for gse, zusatz in KOHORTEN.items():
        ordner = DATEN / gse
        ordner.mkdir(parents=True, exist_ok=True)
        print(f"\n=== {gse} ===")
        hole(f"{FTP}/{praefix(gse)}/{gse}/matrix/{gse}_series_matrix.txt.gz",
             ordner / f"{gse}_series_matrix.txt.gz")
        for datei in zusatz:
            hole(f"{FTP}/{praefix(gse)}/{gse}/suppl/{datei}", ordner / datei)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
11_external_differentiation_systems.py -- Discovery -> Lock -> Validation on INDEPENDENT
differentiation datasets (Figure 2I). NEW in this version of the paper.

The 173-gene programme was discovered and frozen from the 18 perturbation
datasets. Here it is scored, without any refitting, on independent public GEO
differentiation series that took no part in its derivation and are not
perturbation experiments (so they cannot be among the 18). For each dataset the
per-gene log-fold-change (differentiated - undifferentiated) is precomputed and
stored; this script tests the concordance of sign(logFC) with the programme's
frozen directions `ri` against a size-matched background null (the project's
standard concordance test).

Datasets (all independent of the 18; provenance in the panel file):
  GSE37558  primary bone-marrow MSC -> osteoblast (osteogenic), Illumina array
  GSE37558  vascular smooth-muscle cell -> calcifying vascular cell (a
            non-classical, vascular-calcification lineage)
  GSE283759 primary bone-marrow MSC, before -> after adipogenesis (adipogenic)
  GSE214987 iPSC-derived MSC (OA/AC), day 0 -> day 21 chondrogenic

The per-gene log-fold-changes are archived in
derived_data/reference_tables/external_differentiation_logfc.csv.gz (built once
from the public processed matrices; the accessions and grouping are in the
module docstring and the panel file). This script recomputes only the
concordance test, so it needs no network access.

Outputs  results/external_validation.csv and log; figures/data/F2I_external_validation.csv
Runtime  under a minute
"""
from __future__ import annotations
import os
import pathlib
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import _hardening as H  # noqa: E402

_env = os.environ.get("PAPER_V2_ROOT")
ROOT = (pathlib.Path(_env) if _env else pathlib.Path(__file__).resolve().parents[1])
RES = ROOT / "results"
DATA = ROOT / "figures" / "data"
TAB = ROOT / "derived_data" / "reference_tables"

# presentation order (top to bottom in the forest) and labels
LABELS = {
    "GSE37558_osteogenic": ("osteogenic (MSC → osteoblast, GSE37558)", "osteogenic"),
    "GSE283759_adipogenic": ("adipogenic (hMSC, GSE283759)", "adipogenic"),
    "GSE37558_vascular": ("vascular calcification (VSMC → CVC, GSE37558)", "vascular calcification"),
    "GSE214987_chondrogenic": ("chondrogenic (iPSC-derived MSC, GSE214987)", "chondrogenic"),
}


def main():
    mod = pd.read_csv(TAB / "S5_konvergente_gene.csv", dtype={"gen": str})
    ri_sym = mod.set_index("symbol")["ri"]
    ri_ens = mod.set_index("gen")["ri"]
    lf = pd.read_csv(TAB / "external_differentiation_logfc.csv.gz", dtype={"gene": str})

    rows, log = [], []
    for ds, (label, lineage) in LABELS.items():
        d = lf[lf.dataset == ds]
        id_type = d.id_type.iloc[0]
        acc = d.accession.iloc[0]
        vec = d.set_index("gene")["logfc"]
        ri = ri_ens if id_type == "ensembl" else ri_sym
        r = H.concordance(vec.reindex(ri.index).dropna(), ri, background=vec,
                          nziehungen=10000)
        rows.append({
            "dataset": label, "accession": acc, "lineage": lineage,
            "id_type": id_type, "n_tested": r.get("n"),
            "concordance": r.get("concordance"), "z": r.get("z"),
            "mde80": r.get("mde80"), "mde80_z": 2.8,
            "above_mde80": r.get("above_mde80"), "p": r.get("p"),
        })
        log.append(f"  {label:52s} conc {r['concordance']:.3f}  z {r['z']:+.2f}  "
                   f"above_MDE80 {r['above_mde80']}  (n={r['n']})")

    T = pd.DataFrame(rows)
    n_above = int(T.above_mde80.sum())
    header = (f"external validation: {n_above} of {len(T)} independent "
              f"differentiation systems run the locked programme above their limit")
    T.to_csv(RES / "external_validation.csv", index=False)
    T.to_csv(DATA / "F2I_external_validation.csv", index=False)
    (RES / "external_validation_log.txt").write_text(
        header + "\n" + "\n".join(log) + "\n", encoding="utf-8")
    try:
        print(header + "\n" + "\n".join(log))
    except UnicodeEncodeError:
        print((header + "\n" + "\n".join(log)).encode("ascii", "replace").decode())


if __name__ == "__main__":
    main()

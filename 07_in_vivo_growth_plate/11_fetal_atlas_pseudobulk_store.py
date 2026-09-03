# -*- coding: utf-8 -*-
"""
11_fetal_atlas_pseudobulk_store.py -- store the pseudobulk matrix of the limb atlas once.

Purpose  `07_in_vivo_growth_plate/10_fetal_limb_atlas_pseudobulk_build.py` builds
         a pseudobulk mean per (zone, sample) over a fixed block of columns
         (the marker and panel genes needed, plus 4 000 random background
         genes) and then discards it. Every further question put to the same
         matrix -- the gene decomposition of Figure 3C, for instance -- would
         have to load the 7.6 GB atlas again. This script stores the matrix
         **once**.

It recomputes NOTHING and changes nothing in the selection: the same seed, the
same random draw of the background, the same minimum cell count and the same
symbol resolution as in `10_fetal_limb_atlas_pseudobulk_build.py`. That it is the same
matrix is demonstrated rather than asserted: `code/27_...` recomputes the
per-sample module contrast from this file and compares it against
`derived_data/followup/ws4_modulwert_je_probe.csv`.

Inputs   the human fetal limb atlas `.h5ad` (path override: LIMB_ATLAS)
Outputs  results/invivo_pseudobulk.csv.gz  (rows: zone, sample, n_cells;
         columns: gene symbols)
Runtime  a few minutes, about 3 GB of memory
"""
from __future__ import annotations

import os
import pathlib
import sys

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
sys.path.insert(0, str(WURZEL / "00_shared"))
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import MODUL  # noqa: E402

REF = pathlib.Path(os.environ.get(
    "LIMB_ATLAS",
    str(WURZEL.parent / "Referenzdaten" /
        "Limb_Nature2023_s41586-023-06806-x")))
H5AD = REF / "221114LimbCellranger3annotated.minimal.h5ad"
RES = WURZEL / "results"
RES.mkdir(parents=True, exist_ok=True)
AUS = RES / "invivo_pseudobulk.csv.gz"

# --- everything from here on is word for word as in 10_fetal_limb_atlas_pseudobulk_build.py
MIN_ZELLEN = 5
SEED = 20260822
ZONEN = ["MesCond", "ChondroProg", "RestingChon", "ProlifChon",
         "PrehyperChon", "HyperChon"]


def lade_obs(f):
    def cat(name):
        g = f["obs"][name]
        cats = np.array([c.decode() if isinstance(c, bytes) else c
                         for c in g["categories"][:]])
        codes = g["codes"][:]
        return pd.Series(np.where(codes >= 0, cats[np.clip(codes, 0, None)],
                                  "NA"))
    return pd.DataFrame({"celltype": cat("celltype"),
                         "adj_sample": cat("adj_sample"),
                         "adj_stage": cat("adj_stage")})


def main() -> None:
    if AUS.exists() and not os.environ.get("INVIVO_PSEUDOBULK_NEU"):
        print("%s already exists -- nothing done." % AUS)
        print("To rebuild: INVIVO_PSEUDOBULK_NEU=1")
        return

    print("Loading h5ad ...", flush=True)
    f = h5py.File(H5AD, "r")
    obs = lade_obs(f)
    var_index = np.array([g.decode() if isinstance(g, bytes) else g
                          for g in f["var"]["_index"][:]])
    symbol_zu_idx = {}
    for i, s in enumerate(var_index):
        symbol_zu_idx.setdefault(s, i)

    modul_alle = MODUL["symbol"].tolist()
    panels = pd.read_csv(WURZEL / "derived_data" / "M_humangenetik" /
                         "panels.csv")
    pa309 = panels.loc[panels.panel == "PA309",
                       "symbol"].dropna().unique().tolist()
    noso = panels.loc[panels.panel == "NOSO",
                      "symbol"].dropna().unique().tolist()
    benoetigt = sorted(set(modul_alle) | set(OSTEOGEN) | set(ADIPOGEN)
                       | set(MYOGEN) | set(CHONDROGEN) | set(NAIV)
                       | set(pa309) | set(noso))
    benoetigt = [g for g in benoetigt if g in symbol_zu_idx]
    print("Required genes found in the atlas: %d" % len(benoetigt),
          flush=True)

    rng = np.random.default_rng(SEED)
    n_genes = len(var_index)
    hg_idx = rng.choice(n_genes, size=min(4000, n_genes), replace=False)
    hg_symbole = var_index[hg_idx]
    alle_symbole = sorted(set(benoetigt) | set(hg_symbole))
    alle_symbole = [g for g in alle_symbole if g in symbol_zu_idx]
    spalten_idx = np.array([symbol_zu_idx[g] for g in alle_symbole])
    print("Column block: %d genes x %d cells" % (len(spalten_idx), len(obs)),
          flush=True)

    indptr = f["X"]["indptr"][:]
    indices = f["X"]["indices"][:]
    data = f["X"]["data"][:]
    X = sparse.csr_matrix((data, indices, indptr), shape=(len(obs), n_genes))
    del data, indices, indptr
    Xs = X[:, spalten_idx].tocsc()
    del X
    f.close()
    print("Column block loaded.", flush=True)

    zeilen, werte = [], []
    for zone in ZONEN:
        maske_zone = (obs["celltype"].values == zone)
        for p in obs.loc[maske_zone, "adj_sample"].unique():
            m = maske_zone & (obs["adj_sample"].values == p)
            n = int(m.sum())
            if n < MIN_ZELLEN:
                continue
            zeilen.append({"zone": zone, "probe": p, "n_zellen": n})
            werte.append(np.asarray(Xs[m, :].mean(axis=0)).ravel())

    T = pd.concat([pd.DataFrame(zeilen),
                   pd.DataFrame(np.vstack(werte), columns=alle_symbole)],
                  axis=1)
    T.to_csv(AUS, index=False)
    print("%d (zone, sample) points x %d genes" % (len(T), len(alle_symbole)))
    print("-> %s" % AUS)


if __name__ == "__main__":
    main()

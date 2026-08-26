# -*- coding: utf-8 -*-
"""
01_internal_gene_map.py -- the internal RNA layer as a CONTINUOUS quantity
per gene.

So far the paper carries a yes/no list: 173 genes are convergent, all
others are not. For comparison with orthogonal measurement layers this is
the weakest possible form. The underlying matrix
(`_archiv/Sitzungen/20_Exploration/derived_data/20d_gene_*.csv`, 18 points)
holds per gene and dataset `dWT` (differentiation response), `iv` (lesion
interaction term) and `basis` (baseline level). Per gene this script forms:

    dwt_med     median of dWT across datasets
    dwt_v, n    sign majority and number of evaluable datasets
    dwt_kons    (2 * v/n - 1) * sign  -- the signed convergence,
                +1 = up in all datasets, -1 = down in all,
                0 = no common sign
    iv_med, iv_kons   the same for the lesion response
    basis_med   the mean baseline level (nuisance variable, carried along)

This makes the core claim of the paper formulable as a DOSE: does the
strength of a gene's internal convergence predict its movement on a foreign
measurement layer? The 173-gene set then becomes just the upper edge of the
same axis, no longer a category of its own.

Symbol mapping runs via Gencode v46 -- the same annotation the chromatin
layer uses. This keeps the mapping between layers consistent, which would
otherwise silently fail on symbol synonyms.
"""
from __future__ import annotations

import gzip
import pathlib
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import DATEN, ERGEBNISSE, MODUL, WURZEL  # noqa: E402

QUELLE = WURZEL / "_archiv" / "Sitzungen" / "20_Exploration" / "derived_data"
AUS = ERGEBNISSE / "R_intern"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def ensembl_symbol() -> pd.Series:
    """Ensembl gene ID -> gene symbol from Gencode v46 (no version suffix)."""
    zwischen = AUS / "ensembl_symbol_gencode46.csv"
    if zwischen.exists():
        t = pd.read_csv(zwischen)
        return pd.Series(t.symbol.values, index=t.ensembl.values)
    p = DATEN / "_referenz" / "gencode.v46.annotation.hg38.gtf.gz"
    ids, syms = [], []
    rid = re.compile(r'gene_id "([^".]+)')
    rnm = re.compile(r'gene_name "([^"]+)"')
    with gzip.open(p, "rt", encoding="utf-8") as f:
        for z in f:
            if z.startswith("#"):
                continue
            sp = z.split("\t")
            if len(sp) < 9 or sp[2] != "gene":
                continue
            a, b = rid.search(sp[8]), rnm.search(sp[8])
            if a and b:
                ids.append(a.group(1))
                syms.append(b.group(1))
    t = pd.DataFrame({"ensembl": ids, "symbol": syms}).drop_duplicates("ensembl")
    t.to_csv(zwischen, index=False)
    return pd.Series(t.symbol.values, index=t.ensembl.values)


def main() -> None:
    log("=" * 78)
    log("Internal RNA layer as a continuous gene map (18 datasets)")
    log("=" * 78)

    teile = sorted(QUELLE.glob("20d_gene_*.csv"))
    log("Part files: %d" % len(teile))
    G = pd.concat([pd.read_csv(f) for f in teile], ignore_index=True)
    log("Rows: %d | points: %d | genes: %d"
        % (len(G), G.punkt.nunique(), G.gen.nunique()))

    def karte(spalte: str) -> pd.DataFrame:
        P = G.pivot(index="gen", columns="punkt", values=spalte)
        n = P.notna().sum(axis=1)
        pos = (P > 0).sum(axis=1)
        neg = (P < 0).sum(axis=1)
        v = np.maximum(pos, neg)
        vz = np.where(pos >= neg, 1, -1)
        return pd.DataFrame({
            spalte + "_n": n,
            spalte + "_med": P.median(axis=1),
            spalte + "_v": v,
            spalte + "_vz": vz,
            spalte + "_kons": (2 * v / n - 1) * vz,
        })

    K = pd.concat([karte("dWT"), karte("iv")], axis=1)
    K["basis_med"] = G.pivot(index="gen", columns="punkt", values="basis").median(axis=1)
    K["pool_anteil"] = G.pivot(index="gen", columns="punkt", values="im_pool").mean(axis=1)

    MIN_N = 15
    log("\nGenes with >= %d evaluable datasets: %d of %d"
        % (MIN_N, int((K.dWT_n >= MIN_N).sum()), len(K)))
    K = K[K.dWT_n >= MIN_N].copy()

    es = ensembl_symbol()
    K["symbol"] = K.index.map(es)
    log("Symbol mapping via Gencode v46: %d of %d (%.1f %%)"
        % (K.symbol.notna().sum(), len(K), 100 * K.symbol.notna().mean()))
    K["im_module"] = K.index.isin(set(MODUL.ensembl))
    K["ri"] = K.index.map(dict(zip(MODUL.ensembl, MODUL.ri)))

    log("\n--- Check: does the map recover the module? -----------------------")
    log("Module genes in the map: %d of 173" % int(K.im_module.sum()))
    log("dWT_kons in module : median %+.3f (range %+.3f to %+.3f)"
        % (K.dWT_kons[K.im_module].abs().median(),
           K.dWT_kons[K.im_module].min(), K.dWT_kons[K.im_module].max()))
    log("dWT_kons elsewhere : median |%.3f|"
        % K.dWT_kons[~K.im_module].abs().median())
    stimmt = (np.sign(K.dWT_vz[K.im_module]) == K.ri[K.im_module]).mean()
    log("Direction `ri` agrees with dWT_vz: %.1f %% -- must be 100 %%"
        % (100 * stimmt))

    log("\n--- The central contrast, now continuous --------------------------")
    for nm, sp in (("Differentiation dWT", "dWT_kons"), ("Lesion response iv", "iv_kons")):
        a = K[sp].abs()
        log("%-22s |kons| median %.3f | share >= 0.8: %.4f | >= 0.9: %.4f"
            % (nm, a.median(), (a >= 0.8).mean(), (a >= 0.9).mean()))
    log("Same asymmetry as in Figure 4, here without any threshold.")

    K.reset_index().rename(columns={"gen": "ensembl"}).to_csv(
        AUS / "R_interne_genkarte.csv", index=False)
    (AUS / "R_interne_genkarte_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

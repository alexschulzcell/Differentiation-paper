# -*- coding: utf-8 -*-
"""
20_postnatal_growth_plate_annotation.py -- our own clustering of the postnatal
growth-plate series, with cluster-to-zone assignment made EXCLUSIVELY from the
stored literature sets.

Context
-------
The second, postnatal in vivo anchor needs a zone annotation that does not
come from us. The first route was negative: the annotated cell object is not
deposited anywhere (GEO holds only raw count matrices; cellxgene is empty; the
Zenodo record holds only the downstream notebook; the data files of the primary
publication contain no annotation). The second route therefore applies: our
own clustering, but the cluster-to-zone assignment comes exclusively from the
literature sets stored under data_raw/_referenz/wachstumsfuge_zonen/, under a
rule fixed BEFORE the run.

-----------------------------------------------------------------------------
THE DECISION RULE -- written here before the first number was computed, and
not changed afterwards.
-----------------------------------------------------------------------------
1) The score of a Leiden cluster C for zone Z is the mean, over the genes of
   the zone set Z, of their SCALED expression (sc.pp.scale, max 10, over all
   fresh cells).
2) Cluster C is assigned zone Z exactly when
      score(C,Z) = the maximum over the four zones   AND
      score(C,Z) - max(score of the other three) >= DELTA,  DELTA = 0.25.
3) Clusters that the rule does not assign unambiguously are called
   "not assignable" and drop out -- they are NOT corrected by hand. Several
   clusters may carry the same zone.
4) Sensitivity of the threshold (DESCRIPTIVE, decides nothing): DELTA = 0.15
   and 0.50 are logged alongside.
5) Disjointness check (BLOCKING): the genes that ANNOTATE the zones and the
   genes that CALIBRATE the level (chondrogenic against undifferentiated, from
   00_shared/_marker.py) must be disjoint. On a collision the
   CALIBRATION set is reduced, not the annotation set, and the reduced size is
   reported. By project rule 4 calibration sets are NOT filtered by
   measurability; only the collision enters here.
6) The section common to both atlases, fixed before the run: the fetal atlas
   runs MesCond > ChondroProg > Resting > Prolif > Prehyper > Hyper, the
   postnatal growth plate RZ > PZ > (PHZ) > HZ. What is compared is the common
   section Resting<->RZ, Prolif<->PZ, Prehyper<->PHZ, Hyper<->HZ; the direction
   is "immature to hypertrophic" in both. MesCond and ChondroProg have no
   postnatal counterpart and are not compared. If PHZ is missing from the
   points, the ranks of the others are unchanged (Spearman is invariant under
   a strictly monotone rank transformation).
7) Cells: PRIMARILY only the four freshly sequenced human samples. The
   cultured arms (vehicle and growth hormone) are a perturbation experiment
   and do not contribute to the native zone axis. Mouse samples are excluded
   by their genome reference (mm10).
   SENSITIVITY (descriptive, decides nothing; the rule was fixed before the
   second run): with the environment variable POSTNATALE_PROBEN=alle the same
   pipeline runs over ALL twelve human samples (fresh, vehicle and growth
   hormone), and the outputs carry the suffix "_alle". Treatment then stays in
   the data -- ComBat corrects donor means only -- and the verdict of 33_ is
   decoupled from the primary run.

What the script does (standard tools, documented): quality control
(min_genes 200, min_cells 3, at most 15 % mitochondrial reads, no upper bound,
no Scrublet, no regression), normalize_total 1e4, log1p, 2000 highly variable
genes (seurat), ComBat per donor, scale max 10, PCA 30, neighbours k = 15 on
20 principal components, Leiden resolution 1.0, seed 20260824.

Two CORRECTIONS were made before the decisive run, both logged; NO decision
was taken from the aborted first run, which wrote no file before it stopped:
1) PTHR1 (the literature name) is resolved to the official symbol **PTH1R**
   through the HGNC alias list -- pure nomenclature, no change to the content
   of the stored set (`saetze.csv` is unchanged).
2) Clustering WITH ComBat donor correction: the first attempt without batch
   correction produced donor-dominated megaclusters (one cluster held 3 796 of
   3 797 cells from a single donor) that swallowed the zone signal -- exactly
   the problem the primary publication addresses with Harmony or scVI. ComBat
   is the standard method built into scanpy for this purpose. Scores and
   pseudobulks stay on the uncorrected, log-normalised values (they are
   difference-based, and the donor is controlled through the permutation
   within the donor).

Inputs   the unpacked GEO archive of the postnatal series (four fresh samples)
         data_raw/_referenz/wachstumsfuge_zonen/saetze.csv
         00_shared/_marker.py, 00_shared/_module.py
Outputs  results/gse288028_zellannotation.csv.gz
         results/gse288028_cluster_scores.csv
         results/gse288028_pseudobulk.csv.gz
         results/gse288028_annotierung_log.txt
Runtime  a few minutes, about 2 GB of memory
"""
from __future__ import annotations

import os
import pathlib
import sys

import numpy as np
import pandas as pd
import scanpy as sc

_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
sys.path.insert(0, str(WURZEL / "00_shared"))
from _marker import CHONDROGEN, NAIV  # noqa: E402
from _module import MODUL  # noqa: E402

REF = WURZEL.parent / "Referenzdaten" / "GSE288028_ChuSciTranslMed2026" / "RAW"
SATZ_DATEI = WURZEL / "derived_data" / "reference_tables" / "growth_plate_zone_markers.csv"
RES = WURZEL / "results"
RES.mkdir(parents=True, exist_ok=True)

# --- fixed parameters (see the header) --------------------------------------
# All twelve human samples with the assignment taken from the GEO metadata
# (GSM title against file name); the treatment is part of the GSM title.
PROBEN = {
    "GSM9328218": ("GSM9328218_P30453_1001.h5", "P30453",
                   "Human growth plate, rep1"),
    "GSM9328219": ("GSM9328219_P30453_1002.h5", "P30453",
                   "Human growth plate, cultured vehicle, rep1"),
    "GSM9328220": ("GSM9328220_P30453_1003.h5", "P30453",
                   "Human growth plate, cultured gh, rep1"),
    "GSM9328221": ("GSM9328221_P31011_1001.h5", "P31011",
                   "Human growth plate, rep2"),
    "GSM9328222": ("GSM9328222_P31011_1002.h5", "P31011",
                   "Human growth plate, cultured vehicle, rep2"),
    "GSM9328223": ("GSM9328223_P31011_1003.h5", "P31011",
                   "Human growth plate, cultured gh, rep2"),
    "GSM9328224": ("GSM9328224_P25452_001.h5", "P25452",
                   "Human growth plate, rep3"),
    "GSM9328225": ("GSM9328225_P25452_004.h5", "P25452",
                   "Human growth plate, cultured vehicle, rep3 (file1)"),
    "GSM9328226": ("GSM9328226_P25452_005.h5", "P25452",
                   "Human growth plate, cultured vehicle, rep3 (file2)"),
    "GSM9328227": ("GSM9328227_P25452_007.h5", "P25452",
                   "Human growth plate, cultured gh, rep3 (file1)"),
    "GSM9328228": ("GSM9328228_P25452_008.h5", "P25452",
                   "Human growth plate, cultured gh, rep3 (file2)"),
    "GSM9328229": ("GSM9328229_P22202_1015.h5", "P22202",
                   "Human growth plate, rep4"),
}
MODUS = os.environ.get("POSTNATALE_PROBEN", "frisch")
assert MODUS in ("frisch", "alle")
if MODUS == "frisch":
    AUSWAHL = {g: v for g, v in PROBEN.items()
               if "cultured" not in v[2]}
else:
    AUSWAHL = dict(PROBEN)
SUFFIX = "" if MODUS == "frisch" else "_alle"
MIN_ZELLEN = 5           # as in WS4 and 07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py
ALIAS = {"PTHR1": "PTH1R"}   # HGNC official symbol; literature name Lee 1996
DELTA = 0.25             # assignment threshold (rule 2)
DELTA_SENS = [0.15, 0.50]
MIN_GENEN = 200
MIN_ZELLEN_GEN = 3
MAX_MT = 15.0
N_HVG = 2000
N_PCS = 30
N_NACHBARN = 15
N_PCS_KNN = 20
LEIDEN_RES = 1.0
SEED = 20260824          # clustering
SEED_HG = 20260822       # background draw, as in 07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py
N_HG = 4000
RANG = {"RZ": 1, "PZ": 2, "PHZ": 3, "HZ": 4}

LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def main() -> None:
    import anndata as ad

    saetze = pd.read_csv(SATZ_DATEI)
    zonen = list(RANG)

    def loese(g: str) -> str:
        return ALIAS.get(g, g)

    satz_gene = {z: [loese(g) for g in saetze.loc[saetze.zone == z, "gen"]]
                 for z in zonen}
    annot_gene = sorted({loese(g) for g in saetze.gen})
    log("Annotation sets from %s (symbols resolved via the HGNC alias list):"
        % SATZ_DATEI.name)
    for z in zonen:
        log("  %-3s (%d Gene): %s" % (z, len(satz_gene[z]), ", ".join(satz_gene[z])))

    # ---- disjointness check (blocking, rule 5) -----------------------------
    eich_union = sorted(set(CHONDROGEN) | set(NAIV))
    kollision = sorted(set(annot_gene) & set(eich_union))
    log("")
    log("Disjointness check (blocking):")
    log("  annotation genes: %d | calibration union (CHONDROGEN & NAIV): %d"
        % (len(annot_gene), len(eich_union)))
    log("  intersection: %s" % (kollision if kollision else "EMPTY"))
    naiv_red = [g for g in NAIV if g not in kollision]
    chondrogen_red = [g for g in CHONDROGEN if g not in kollision]
    log("  calibration set after collision reduction: CHONDROGEN %d genes "
        "(unchanged) | NAIV %d -> %d genes"
        % (len(chondrogen_red), len(NAIV), len(naiv_red)))
    if len(naiv_red) < 3 or len(chondrogen_red) < 3:
        raise SystemExit("Calibration sets too small after reduction -- "
                         "abort.")

    # ---- loading -----------------------------------------------------------
    adatas = []
    for gsm, (datei, spender, titel) in AUSWAHL.items():
        a = sc.read_10x_h5(REF / datei, gex_only=True)
        a.var_names_make_unique()
        a.obs_names_make_unique()
        a.obs["gsm"] = gsm
        a.obs["probe"] = datei[:-3]
        a.obs["spender"] = spender
        a.obs["geo_titel"] = titel
        genome = set(np.asarray(a.var["genome"]).tolist())
        if not genome <= {"GRCh38"}:
            raise SystemExit("%s: unexpected genome reference %s"
                             % (gsm, genome))
        log("%s (%s, donor %s): %d cells x %d genes [%s]"
            % (gsm, titel, spender, a.n_obs, a.n_vars, datei))
        adatas.append(a)
    A = ad.concat(adatas, join="outer", label="quelle")
    A.var_names_make_unique()
    A.obs_names_make_unique()
    log("total: %d cells x %d genes" % (A.n_obs, A.n_vars))

    # ---- QC ---------------------------------------------------------------
    A.var["mt"] = A.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    n0 = A.n_obs
    A = A[(A.obs["n_genes_by_counts"] >= MIN_GENEN)
          & (A.obs["pct_counts_mt"] <= MAX_MT)].copy()
    sc.pp.filter_genes(A, min_cells=MIN_ZELLEN_GEN)
    log("QC: %d -> %d cells | %d genes remain (min %d cells)"
        % (n0, A.n_obs, A.n_vars, MIN_ZELLEN_GEN))

    # ---- pipeline (see the header) -----------------------------------------
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=N_HVG, flavor="seurat")
    H = A[:, A.var.highly_variable].copy()
    sc.pp.combat(H, key="spender")
    H.X = np.array(H.X)          # ComBat returns a read-only array
    sc.pp.scale(H, max_value=10)
    sc.tl.pca(H, n_comps=N_PCS, random_state=SEED)
    sc.pp.neighbors(H, n_neighbors=N_NACHBARN, n_pcs=N_PCS_KNN,
                    random_state=SEED)
    sc.tl.leiden(H, resolution=LEIDEN_RES, key_added="leiden",
                 random_state=SEED, flavor="igraph", n_iterations=2,
                 directed=False)
    A.obs["leiden"] = H.obs["leiden"].reindex(A.obs_names).astype(str)
    log("Leiden (%.1f): %d clusters, sizes %s"
        % (LEIDEN_RES, A.obs["leiden"].nunique(),
           dict(A.obs["leiden"].value_counts().sort_index())))

    # ---- scores per cluster (rule 1) ---------------------------------------
    # Implementation: z-scaling PER marker gene over ALL fresh cells (on the
    # log-normalised values), then the cell mean per cluster, then the mean
    # over the genes of the set. The presence check runs against the full gene
    # complement, not against the list of highly variable genes.
    symbol_idx_full = {s: i for i, s in enumerate(A.var_names)}
    alle_marker = sorted(set(annot_gene) & set(symbol_idx_full))
    subX = A[:, alle_marker].X
    if hasattr(subX, "toarray"):
        subX = subX.toarray()
    M = pd.DataFrame(np.asarray(subX), index=A.obs_names, columns=alle_marker)
    M = (M - M.mean(axis=0)) / M.std(axis=0, ddof=0).replace(0.0, np.nan)
    zeilen = []
    leiden_v = A.obs["leiden"].to_numpy()
    for c in sorted(A.obs["leiden"].unique(), key=int):
        maske_c = (leiden_v == c)
        zeile = {"cluster": c, "n_zellen": int(maske_c.sum())}
        cm = M[maske_c].mean(axis=0)
        for z in zonen:
            gene_da = [g for g in satz_gene[z] if g in symbol_idx_full]
            gene_skalierbar = [g for g in gene_da
                               if g in cm.index and np.isfinite(cm[g])]
            zeile[z] = float(cm[gene_skalierbar].mean()) \
                if gene_skalierbar else np.nan
            zeile[z + "_gene_im_datensatz"] = len(gene_da)
        zeilen.append(zeile)
    S = pd.DataFrame(zeilen).set_index("cluster")
    log("")
    log("Marker presence in the full data set:")
    for z in zonen:
        fehlen = [g for g in satz_gene[z] if g not in symbol_idx_full]
        log("  %-3s: %d of %d genes present%s"
            % (z, len(satz_gene[z]) - len(fehlen), len(satz_gene[z]),
               "" if not fehlen else " -- missing: " + ", ".join(fehlen)))
    log("Cluster scores (z-scaled unit):")
    log(S[[z for z in zonen]].round(4).to_string())
    S.to_csv(RES / ("gse288028_cluster_scores%s.csv" % SUFFIX))

    # ---- assignment rule (rules 2-4) ---------------------------------------
    def zuordnung(delta: float) -> dict:
        out = {}
        for c, r in S.iterrows():
            werte = {z: r[z] for z in zonen if np.isfinite(r[z])}
            if len(werte) < 4 or any(r[z + "_gene_im_datensatz"] == 0
                                     for z in zonen):
                out[c] = "nicht zuordenbar"
                continue
            best = max(werte, key=werte.get)
            rest = [v for z, v in werte.items() if z != best]
            out[c] = best if (werte[best] - max(rest)) >= delta else \
                "nicht zuordenbar"
        return out

    zuw = zuordnung(DELTA)
    A.obs["zone"] = A.obs["leiden"].map(zuw).astype("string")
    log("")
    log("Assignment at DELTA=%.2f: %s" % (DELTA, zuw))
    for d in DELTA_SENS:
        log("Sensitivity DELTA=%.2f: %s" % (d, zuordnung(d)))
    kt = pd.crosstab(A.obs["leiden"], A.obs["spender"])
    log("")
    log("Cluster x donor (diagnostics, no batch correction):")
    log(kt.to_string())
    kz = pd.crosstab(A.obs["zone"], A.obs["spender"])
    log("")
    log("Zone x donor (cells):")
    log(kz.to_string())

    # ---- storing -----------------------------------------------------------
    A.obs[["gsm", "probe", "spender", "geo_titel", "leiden", "zone"]] \
        .to_csv(RES / ("gse288028_zellannotation%s.csv.gz" % SUFFIX))

    # ---- pseudobulk matrix (analogous to 07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py) ----
    modul_alle = MODUL["symbol"].tolist()
    benoetigt = sorted(set(modul_alle) | set(chondrogen_red) | set(naiv_red))
    symbol_idx = {s: i for i, s in enumerate(A.var_names)}
    fehlend = [g for g in benoetigt if g not in symbol_idx]
    benoetigt_da = [g for g in benoetigt if g in symbol_idx]
    log("")
    log("Required genes present in the data set: %d of %d%s"
        % (len(benoetigt_da), len(benoetigt),
           "" if not fehlend else " (missing: %s)" % ", ".join(fehlend[:10])))
    rng = np.random.default_rng(SEED_HG)
    hg_idx = rng.choice(A.n_vars, size=min(N_HG, A.n_vars), replace=False)
    hg_symbole = A.var_names[hg_idx]
    alle_symbole = sorted(set(benoetigt_da) | set(hg_symbole))
    spalten = [symbol_idx[g] for g in alle_symbole]

    X = A.X
    if hasattr(X, "tocsr"):
        X = X.tocsr()
    zone_v = np.asarray(A.obs["zone"].astype(object))
    spender_v = np.asarray(A.obs["spender"].astype(object))
    zeilen_pb, werte_pb = [], []
    for zone in zonen:
        maske_zone = (zone_v == zone)
        if not maske_zone.any():
            log("Zone %s: not represented at all -- dropped." % zone)
            continue
        for p in pd.unique(A.obs.loc[maske_zone, "spender"]):
            m = maske_zone & (spender_v == p)
            n = int(m.sum())
            if n < MIN_ZELLEN:
                log("Zone %s, donor %s: only %d cells (<%d) -- dropped."
                    % (zone, p, n, MIN_ZELLEN))
                continue
            sub = X[m][:, spalten]
            mw = np.asarray(sub.mean(axis=0)).ravel()
            zeilen_pb.append({"zone": zone, "spender": p, "n_zellen": n})
            werte_pb.append(mw)
    if not zeilen_pb:
        raise SystemExit("No (zone, donor) points -- abort.")
    T = pd.concat([pd.DataFrame(zeilen_pb),
                   pd.DataFrame(np.vstack(werte_pb), columns=alle_symbole)],
                  axis=1)
    T.to_csv(RES / ("gse288028_pseudobulk%s.csv.gz" % SUFFIX), index=False)
    log("Pseudobulk: %d points x %d genes -> gse288028_pseudobulk%s.csv.gz"
        % (len(T), len(alle_symbole), SUFFIX))

    (RES / ("gse288028_annotierung_log%s.txt" % SUFFIX)).write_text(
        chr(10).join(LOG) + chr(10), encoding="utf-8")


if __name__ == "__main__":
    main()

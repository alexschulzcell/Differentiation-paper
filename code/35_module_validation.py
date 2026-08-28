# -*- coding: utf-8 -*-
"""
35_module_validation.py -- held-out validation and robustness of the 173-gene
programme. NEW in this version of the paper.

The programme was discovered from all 18 datasets and, in v1, tested in the
same 18 (explicitly exploratory). This script adds three independent checks
that the programme is real, generalises, and is not an artefact of its
derivation or of a few genes:

  LEAVE-ONE-STUDY-OUT (Figure 2G). For each held-out GEO study, the convergence
     programme is re-derived from the REMAINING datasets only (same frozen
     rule, threshold scaled to the training-set size) and scored, untouched,
     on the held-out dataset(s). Datasets of one publication never sit in both
     training and test. Reported: held-out concordance z per dataset and how
     many of 18 exceed their own detection limit. Leave-one-dataset-out is run
     too, for completeness.

  MATCHED RANDOM NULLS (Figure 2A, second null). The cross-arm Spearman rho of
     the arm-pooled dWT vectors on the programme (the Fig. 2A statistic, 0.622)
     is compared with 10,000 draws of 173 genes matched to the programme's
     joint distribution of baseline expression rank, gene length and constraint
     (LOEUF) -- the "you just rediscovered well-studied ECM genes" objection
     made quantitative. Drawn on the same panel as the unmatched null.

  GENE ROBUSTNESS (Figure 2H). The same cross-arm rho recomputed after dropping
     the strongest / most-expressed / random gene fractions, plus a
     leave-one-gene-out jackknife over all 173 genes.

Inputs   derived_data/reference_tables/{20d_dWT_matrix.csv.gz,
         dataset_study_map.csv, S5_konvergente_gene.csv},
         derived_data/followup/ws1_genkarte_erweitert.csv
Outputs  results/module_validation_*.csv and log; figures/data/F2G_*,
         F2A_matched_nulls_*, F2H_dropout/jackknife panel files
Runtime  about two to three minutes
"""
from __future__ import annotations
import os
import pathlib
import numpy as np
import pandas as pd
from scipy import stats

import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import _hardening as H  # noqa: E402

_env = os.environ.get("PAPER_V2_ROOT")
ROOT = (pathlib.Path(_env) if _env else pathlib.Path(__file__).resolve().parents[1])
RES = ROOT / "results"
DATA = ROOT / "figures" / "data"
NDRAW = 10000
SEED = 20260827
LOG = []


def log(s=""):
    print(s)
    LOG.append(str(s))


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a or b) else np.nan


# ---------------------------------------------------------------- G ----------
def leave_out(wide, smap, group_col, tag):
    all_pts = list(wide.columns)
    frozen = set(pd.read_csv(H.paths()["s5"], dtype={"gen": str}).gen)
    if group_col == "gse":
        groups = smap.groupby("gse").punkt.apply(list).to_dict()
    else:
        groups = {int(p): [int(p)] for p in all_pts}
    rows = []
    for gname, held in sorted(groups.items(), key=lambda kv: min(kv[1])):
        train = [p for p in all_pts if p not in held]
        mod = H.derive_module(wide, train)
        ri = mod["ri"]
        rederived = set(mod.index)
        for p in held:
            dwt = wide[p].dropna()
            r = H.concordance(dwt.reindex(ri.index).dropna(), ri, background=dwt,
                              nziehungen=5000)
            ds = smap.set_index("punkt").loc[p]
            rows.append({
                "held_out_unit": gname if group_col == "gse" else ds["datensatz"],
                "punkt": p, "datensatz": ds["datensatz"], "gse": ds["gse"],
                "arm": ds["arm"], "n_train_datasets": len(train),
                "n_rederived_genes": len(rederived),
                "jaccard_with_frozen173": round(jaccard(rederived, frozen), 4),
                "n_tested": r.get("n"), "concordance": r.get("concordance"),
                "z": r.get("z"), "mde80_z": 2.8, "above_mde80": r.get("above_mde80"),
            })
    T = pd.DataFrame(rows)
    ok = T[T.z.notna()]
    log(f"  {tag}: {int(ok.above_mde80.sum())}/{len(ok)} held-out datasets above own MDE80; "
        f"median z {ok.z.median():+.2f} (range {ok.z.min():+.2f}..{ok.z.max():+.2f}); "
        f"re-derived genes median {int(ok.n_rederived_genes.median())}, "
        f"Jaccard vs 173 median {ok.jaccard_with_frozen173.median():.3f}")
    return T


# ------------------------------------------------------ arm vectors ----------
def arm_vectors(wide, smap):
    a = smap.set_index("punkt")
    osteo = [p for p in wide.columns if a.loc[p, "arm"] == "osteogen"]
    chond = [p for p in wide.columns if a.loc[p, "arm"] == "chondrogen"]
    O = wide[osteo].median(axis=1, skipna=True)
    C = wide[chond].median(axis=1, skipna=True)
    gem = O.dropna().index.intersection(C.dropna().index)
    return O.reindex(gem), C.reindex(gem), len(osteo), len(chond)


def rho_on(Ov, Cv, idx):
    return stats.spearmanr(Ov[idx], Cv[idx])[0]


# ---------------------------------------------------------------- H ----------
def matched_nulls(wide, smap, genkarte, module):
    O, C, n_ost, n_ch = arm_vectors(wide, smap)
    background = list(O.index)
    pos = {g: i for i, g in enumerate(background)}
    Ov, Cv = O.values, C.values
    rng = np.random.default_rng(SEED)

    mod_idx = [pos[g] for g in module if g in pos]
    k = len(mod_idx)
    obs = rho_on(Ov, Cv, mod_idx)

    # unmatched
    allidx = np.arange(len(background))
    un = np.array([rho_on(Ov, Cv, rng.choice(allidx, size=k, replace=False))
                   for _ in range(NDRAW)])
    z_un = (obs - un.mean()) / un.std(ddof=1)
    p_un = (1 + int((un >= obs).sum())) / (1 + NDRAW)

    # matched on expr rank + length + LOEUF, decile-of-5 joint cells
    gk = genkarte.set_index("ensembl").reindex(background)[
        ["expr_rank_med", "laenge", "loeuf"]].dropna()
    codes = {v: pd.qcut(gk[v].rank(method="first"), 5, labels=False)
             for v in ("expr_rank_med", "laenge", "loeuf")}
    cell = (codes["expr_rank_med"].astype(str) + "_" + codes["laenge"].astype(str)
            + "_" + codes["loeuf"].astype(str))
    cell = pd.Series(cell, index=gk.index)
    cell2pos = {c: np.array([pos[g] for g in gs])
                for c, gs in cell.groupby(cell).groups.items()}
    mod_cells = cell.reindex([g for g in module if g in cell.index]).dropna()
    counts = mod_cells.value_counts()
    ma = np.empty(NDRAW)
    for i in range(NDRAW):
        idx = np.concatenate([rng.choice(cell2pos[c], size=cnt,
                                         replace=len(cell2pos[c]) < cnt)
                              for c, cnt in counts.items()])
        ma[i] = rho_on(Ov, Cv, idx)
    obs_m = rho_on(Ov, Cv, [pos[g] for g in mod_cells.index])
    z_m = (obs_m - ma.mean()) / ma.std(ddof=1)
    p_m = (1 + int((ma >= obs_m).sum())) / (1 + NDRAW)

    summary = pd.DataFrame([
        {"null_type": "unmatched", "n_genes": k, "rho_observed": obs,
         "null_mean": un.mean(), "null_sd": un.std(ddof=1), "z_sd_units": z_un,
         "p_empirical": p_un, "n_background": len(background)},
        {"null_type": "matched", "n_genes": int(counts.sum()), "rho_observed": obs_m,
         "null_mean": ma.mean(), "null_sd": ma.std(ddof=1), "z_sd_units": z_m,
         "p_empirical": p_m, "n_background": len(gk)},
    ])
    draws = pd.concat([
        pd.DataFrame({"null_type": "unmatched", "rho": un}),
        pd.DataFrame({"null_type": "matched", "rho": ma}),
    ], ignore_index=True)
    log(f"  matched nulls: observed rho {obs_m:.3f} vs matched null "
        f"{ma.mean():.3f}+-{ma.std(ddof=1):.3f}  z {z_m:+.2f}  p {p_m:.2g} "
        f"(unmatched z {z_un:+.2f})")
    return summary, draws


# ---------------------------------------------------------------- I ----------
def robustness(wide, smap, genkarte, module):
    O, C, _, _ = arm_vectors(wide, smap)
    genes = [g for g in module.gen if g in O.index]
    Os, Cs = O[genes].values, C[genes].values
    n = len(genes)

    def rho(mask):
        return stats.spearmanr(Os[mask], Cs[mask])[0]

    full = rho(np.ones(n, bool))
    abs_dwt = np.abs(np.nanmedian(np.vstack([Os, Cs]), axis=0))
    expr = genkarte.set_index("ensembl").reindex(genes)["expr_rank_med"].values
    rng = np.random.default_rng(SEED)
    rows = [{"scheme": "full", "removed_frac": 0.0, "n_used": n, "rho": full, "rho_sd": np.nan}]
    for frac in (0.01, 0.05, 0.10, 0.20):
        nd = max(1, int(round(frac * n)))
        for key, order, name in (
            (abs_dwt, np.argsort(-abs_dwt), "drop_top_abs_dwt"),
            (expr, np.argsort(-np.nan_to_num(expr, nan=-1)), "drop_top_expression")):
            m = np.ones(n, bool); m[order[:nd]] = False
            rows.append({"scheme": name, "removed_frac": frac,
                         "n_used": int(m.sum()), "rho": rho(m), "rho_sd": np.nan})
        vals = []
        for _ in range(500):
            m = np.ones(n, bool); m[rng.choice(n, size=nd, replace=False)] = False
            vals.append(rho(m))
        rows.append({"scheme": "drop_random", "removed_frac": frac, "n_used": n - nd,
                     "rho": float(np.mean(vals)), "rho_sd": float(np.std(vals, ddof=1))})
    dropout = pd.DataFrame(rows)

    jk = np.empty(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        jk[i] = rho(m)
    jack = pd.DataFrame({
        "ensembl": genes,
        "symbol": module.set_index("gen").reindex(genes)["symbol"].values,
        "rho_without_gene": jk, "delta_from_full": jk - full}).sort_values("delta_from_full")
    log(f"  robustness: full rho {full:.3f}; top-20% |dWT| removed -> "
        f"{dropout.query('scheme==\"drop_top_abs_dwt\" and removed_frac==0.20').rho.iloc[0]:.3f}; "
        f"jackknife range {jk.min():.3f}..{jk.max():.3f} "
        f"(largest single-gene drop {full - jk.max():+.3f})")
    return dropout, jack, full


def main():
    wide = H.load_dwt_wide()
    smap = H.study_map()
    genkarte = pd.read_csv(H.paths()["genkarte"])
    module = pd.read_csv(H.paths()["s5"], dtype={"gen": str})

    # self-test: the full re-derivation must reproduce the frozen 173
    full_mod = H.derive_module(wide, list(wide.columns))
    assert set(full_mod.index) == set(module.gen), "re-derivation != frozen S5"
    log("self-test: full re-derivation reproduces the frozen 173 genes -- OK")

    log("\nG  leave-one-study-out / leave-one-dataset-out")
    Tstudy = leave_out(wide, smap, "gse", "leave-one-study-out")
    Tds = leave_out(wide, smap, "dataset", "leave-one-dataset-out")
    Tstudy.to_csv(RES / "module_validation_leave_one_study_out.csv", index=False)
    Tds.to_csv(RES / "module_validation_leave_one_dataset_out.csv", index=False)
    # panel file (primary = study level), sorted for the forest. The arm is
    # translated to English for the published panel file (project language rule).
    g = Tstudy.copy()
    g["arm"] = g["arm"].map({"osteogen": "osteogenic",
                             "chondrogen": "chondrogenic"}).fillna(g["arm"])
    g[["datensatz", "punkt", "arm", "n_tested", "concordance", "z", "mde80_z",
       "above_mde80", "n_rederived_genes", "jaccard_with_frozen173"]].to_csv(
        DATA / "F2G_leave_one_study_out.csv", index=False)

    log("\nH  matched random nulls")
    summ, draws = matched_nulls(wide, smap, genkarte, list(module.gen))
    summ.to_csv(RES / "module_validation_matched_nulls.csv", index=False)
    # the matched null joins the cross-arm concordance panel (Figure 2A)
    summ.to_csv(DATA / "F2A_matched_nulls_summary.csv", index=False)
    draws.to_csv(DATA / "F2A_matched_nulls_draws.csv", index=False)

    log("\nH  gene robustness")
    dropout, jack, full = robustness(wide, smap, genkarte, module)
    dropout.to_csv(RES / "module_validation_dropout.csv", index=False)
    jack.to_csv(RES / "module_validation_jackknife.csv", index=False)
    dropout.to_csv(DATA / "F2H_dropout.csv", index=False)
    jack.to_csv(DATA / "F2H_jackknife.csv", index=False)

    (RES / "module_validation_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    log("\nwritten: results/module_validation_*.csv, figures/data/F2G_*, F2A_matched_nulls_*, F2H_dropout/jackknife")


if __name__ == "__main__":
    main()

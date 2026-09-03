# -*- coding: utf-8 -*-
"""Follow-up Phase 4 data provenance: write derived_data/manuscript/f6_followup_*.csv
from the frozen follow-up analysis artifacts (26_Orthogonal_S12/Ergebnisse).
No re-analysis; pure export + derived descriptive aggregates (explicitly
exploratory)."""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd

# --- path parameters (2026-08-23) -----------------------------------------
# Two roots, kept separate because they are not the same tree:
#   * the repository itself, which holds derived_data/ and results/ --
#     found from this file's location, or from PAPER_V2_ROOT if set;
#   * the author's session folders, which hold the per-run intermediates
#     of the original exploration and are not part of this repository --
#     SCHERENPAPER_SITZUNGEN if set. Steps that need them are marked
#     "needs raw data" in the stage README.
import os, pathlib
_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
SITZUNGEN = pathlib.Path(os.environ.get(
    "SCHERENPAPER_SITZUNGEN", str(WURZEL / "_archiv" / "Sitzungen")))
ERG = SITZUNGEN / "26_Orthogonal_S12" / "derived_data"
DAT = WURZEL / "derived_data" / "manuscript"

ACC = ["GSE247491", "GSE247528"]
LINE = {  # already frozen in the follow-up Outputs
    "GSE247491": "Chondrogenesis (SERPINA3 cohort)",
    "GSE247528": "Osteogenesis (SERPINA3 cohort)",
}


def main() -> int:
    DAT.mkdir(parents=True, exist_ok=True)
    s = pd.read_csv(ERG / "followup_summary.csv", dtype=str)

    # atlas rows (study-level primary endpoints)
    atlas = []
    for acc in ACC:
        d = s[(s.linie == acc) & (s.endpoint == "atlas_day7")].iloc[0]
        d3 = s[(s.linie == acc) & (s.endpoint == "atlas_day3")].iloc[0]
        atlas.append(dict(
            dataset=acc, label=LINE[acc], n_units=3,
            share=float(d["share_match"]), n_matching=int(float(d["n_matching"])),
            n_measurable=int(float(d["n_measurable"])),
            wilson_lo=float(d["wilson_lo"]), wilson_hi=float(d["wilson_hi"]),
            z=float(d["z"]), p_perm=float(d["p_perm_two_sided"]),
            null_mean=float(d["perm_mean"]), null_sd=float(d["perm_sd"]),
            mde80=float(d["mde80"]), above_mde80=(d["above_mde80"] == "True"),
            loo_min=float(d["loo_min"]), loo_max=float(d["loo_max"]),
            not_single_unit=(d["not_single_unit"] == "True"),
            share_d3=float(d3["share_match"]), mde80_d3=float(d3["mde80"]),
            above_mde80_d3=(d3["above_mde80"] == "True"),
        ))
    atlas_df = pd.DataFrame(atlas)
    atlas_df.to_csv(DAT / "f6_followup_atlas.csv", index=False)

    # iv rows (2x2 convergent)
    iv = []
    for acc in ACC:
        d = s[(s.linie == acc) & (s.endpoint == "iv_day7")].iloc[0]
        d3 = s[(s.linie == acc) & (s.endpoint == "iv_day3")].iloc[0]
        iv.append(dict(
            dataset=acc, label=LINE[acc], n_units=3,
            unanim_observed=int(float(d["unanimous_observed"])),
            unanim_null_mean=float(d["unanimous_null_mean"]),
            unanim_null_sd=float(d["unanimous_null_sd"]),
            unanim_z=float(d["unanimous_z"]), unanim_p=float(d["unanimous_p"]),
            unanim_n_measurable=int(float(d["n_measurable"])),
            major_observed=int(float(d["majority_observed"])),
            major_null_mean=float(d["majority_null_mean"]),
            major_null_sd=float(d["majority_null_sd"]),
            major_p=float(d["majority_p"]),
            unanim_observed_d3=int(float(d3["unanimous_observed"])),
            unanim_null_mean_d3=float(d3["unanimous_null_mean"]),
            unanim_p_d3=float(d3["unanimous_p"]),
        ))
    iv_df = pd.DataFrame(iv)
    iv_df.to_csv(DAT / "f6_followup_iv.csv", index=False)

    # donor module score secondary
    score_rows = []
    for acc in ACC:
        m = pd.read_csv(ERG / f"followup_{acc}_module_score.csv", dtype=str).iloc[0]
        score_rows.append(dict(
            dataset=acc, label=LINE[acc], n_donors=int(m["n_donors"]),
            t_stat=float(m["t_stat"]), t_p=float(m["t_p"]),
            wilcoxon_p=float(m["wilcoxon_p"]), mean_diff=float(m["mean_diff"]),
            seed=int(m["seed"]), score_per_donor=m["score_per_donor"],
        ))
    pd.DataFrame(score_rows).to_csv(DAT / "f6_followup_module_score.csv", index=False)

    # descriptive aggregates (not a registered endpoint; for the integration map)
    pooled_rows = []
    for acc in ACC:
        g = pd.read_csv(ERG / f"followup_{acc}_gene_level.csv")
        pooled_rows.append(g.assign(dataset=acc))
    allg = pd.concat(pooled_rows, ignore_index=True)
    agg = dict(
        n_studies=len(ACC),
        n_units=int(allg.donor.nunique()),
        n_observations=int(len(allg)),
        iv_unanim_combined=int(iv_df.unanim_observed.sum()),
        iv_unanim_null_combined=float(iv_df.unanim_null_mean.sum()),
    )
    pg = allg.groupby(["dataset", "donor"]).apply(
        lambda x: float(np.sign(x.dWT).eq(x.ri).mean()), include_groups=False)
    agg["unit_shares"] = ";" .join(
        f"{a}:{d}:{v:.3f}" for (a, d), v in pg.items())
    pd.DataFrame([agg]).to_csv(DAT / "f6_followup_pooled.csv", index=False)
    print(pd.DataFrame([agg]).to_string(index=False))
    print("wrote f6_followup_*.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
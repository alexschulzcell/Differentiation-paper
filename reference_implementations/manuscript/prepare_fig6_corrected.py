# -*- coding: utf-8 -*-
"""Prepare provenance-checked data for the corrected Fig. 6."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

# --- path parameters (2026-08-23) -----------------------------------------
# Previously a hard-coded path pointed to ".../DFG Antrag/
# Scherenpaper_Folgeprojekt" -- a Windows junction, not clonable. Two
# different trees sat in this one variable:
#   * the session folders 25_Orthogonal_S11 / 26_Orthogonal_S12 -- they live
#     under `Paper v2/_archiv/Sitzungen/`. The audit of 2026-08-22 had listed
#     them as "missing"; they were just stored elsewhere.
#   * `derived_data/manuscript`, `derived_data/reference_tables` -- these live under `Paper v2` itself.
# Both are now separate and overridable via environment variables.
import os, pathlib
_env = os.environ.get("PAPER_V2_ROOT")
WURZEL = (pathlib.Path(_env) if _env
          else pathlib.Path(__file__).resolve().parents[1])
SITZUNGEN = pathlib.Path(os.environ.get(
    "SCHERENPAPER_SITZUNGEN", str(WURZEL / "_archiv" / "Sitzungen")))
DAT = WURZEL / "derived_data" / "manuscript"
ERG = SITZUNGEN / "26_Orthogonal_S12" / "derived_data"
EXT = SITZUNGEN / "26_Orthogonal_S12" / "data_raw" / "_followup_extension"


def atlas_rows() -> pd.DataFrame:
    s12 = pd.read_csv(DAT / "f6_s12_fixed173_by_study.csv")
    s12["source"] = "S12 corrected fixed-173"
    s12["cohort"] = s12["dataset"]
    s12["label"] = "S12 " + s12["dataset"]
    s12["p"] = s12["p_perm_two_sided"]
    s12["note"] = np.where(s12["dataset"] == "GSE37521", "sensitivity: n=2", "independent S12 study")

    follow = pd.read_csv(DAT / "f6_followup_atlas.csv")
    follow["source"] = "SERPINA3 follow-up"
    follow["cohort"] = "SERPINA3 shared donor cohort"
    follow["label"] = follow["dataset"].map({
        "GSE247491": "SERPINA3 chondro (GSE247491)",
        "GSE247528": "SERPINA3 osteo (GSE247528)",
    })
    follow["p"] = follow["p_perm"]
    follow["note"] = "one shared donor cohort"
    follow = follow.rename(columns={"n_units": "n_biological_units", "share": "share_match"})

    ext = pd.read_csv(ERG / "followup_extension_atlas_summary.csv")
    ext["source"] = "search extension"
    ext["cohort"] = ext["dataset"]
    ext["label"] = ext["dataset"].map({
        "GSE12266": "hMSC mineralization (GSE12266)",
        "GSE18043": "BM-MSC osteo (GSE18043)",
        "GSE63754": "adipose-MSC osteo (GSE63754)",
    })
    ext["p"] = ext["p_perm_two_sided"]
    ext["note"] = "independent extension cohort"
    ext["wilson_lo"] = np.nan
    ext["wilson_hi"] = np.nan

    cols = ["dataset", "label", "source", "cohort", "n_biological_units",
            "n_measurable", "n_matching", "share_match", "wilson_lo", "wilson_hi",
            "null_mean", "null_sd", "z", "p", "mde80", "above_mde80", "note"]
    s12 = s12.rename(columns={"n_units": "n_biological_units"})
    # S12 fixed-173 data have no Wilson interval; leave it missing rather than
    # constructing an independent-gene interval for the corrected figure.
    s12["wilson_lo"] = np.nan
    s12["wilson_hi"] = np.nan
    return pd.concat([s12[cols], follow[cols], ext[cols]], ignore_index=True)


def study_level_synthesis() -> tuple[pd.DataFrame, pd.DataFrame]:
    s12 = pd.read_csv(DAT / "f6_s12_fixed173_by_study.csv")
    s12_null = pd.read_csv(DAT / "f6_s12_fixed173_study_null.csv")
    ext = pd.read_csv(ERG / "followup_extension_atlas_summary.csv")
    ext_null = {
        acc: pd.read_csv(ERG / f"followup_extension_{acc}_null.csv")["share_null"].to_numpy()
        for acc in ext.dataset
    }

    records = []
    null_z = {}
    for row in s12.to_dict("records"):
        acc = row["dataset"]
        arr = s12_null[s12_null.dataset == acc]["share_null"].to_numpy()
        null_z[acc] = (arr - row["null_mean"]) / row["null_sd"]
        records.append({"dataset": acc, "n_units": row["n_units"],
                        "observed_z": row["z"], "source": "S12 corrected fixed-173"})
    for row in ext.to_dict("records"):
        acc = row["dataset"]
        null_z[acc] = (ext_null[acc] - row["null_mean"]) / row["null_sd"]
        records.append({"dataset": acc, "n_units": row["n_biological_units"],
                        "observed_z": row["z"], "source": "search extension"})

    by = pd.DataFrame(records)
    primary_ids = by.loc[by.n_units >= 3, "dataset"].tolist()
    all_ids = by["dataset"].tolist()
    primary_obs = float(by[by.dataset.isin(primary_ids)].observed_z.mean())
    sensitivity_obs = float(by.observed_z.mean())
    primary_null = np.mean(np.vstack([null_z[a] for a in primary_ids]), axis=0)
    sensitivity_null = np.mean(np.vstack([null_z[a] for a in all_ids]), axis=0)
    summaries = pd.DataFrame([
        {"analysis": "primary n>=3 studies", "n_studies": len(primary_ids),
         "omitted_sensitivity": "GSE37521", "observed_stat": primary_obs,
         "null_mean": float(primary_null.mean()), "null_sd": float(primary_null.std(ddof=1)),
         "z": (primary_obs - primary_null.mean()) / primary_null.std(ddof=1),
         "p_upper": float(np.mean(primary_null >= primary_obs)), "seed": 20260829},
        {"analysis": "sensitivity all S12+extension", "n_studies": len(all_ids),
         "omitted_sensitivity": "none", "observed_stat": sensitivity_obs,
         "null_mean": float(sensitivity_null.mean()), "null_sd": float(sensitivity_null.std(ddof=1)),
         "z": (sensitivity_obs - sensitivity_null.mean()) / sensitivity_null.std(ddof=1),
         "p_upper": float(np.mean(sensitivity_null >= sensitivity_obs)), "seed": 20260829},
    ])
    null = pd.DataFrame({"primary_stat_null": primary_null,
                         "sensitivity_stat_null": sensitivity_null})
    by.to_csv(DAT / "f6_study_level_by_study.csv", index=False)
    summaries.to_csv(DAT / "f6_study_level_summary.csv", index=False)
    null.to_csv(DAT / "f6_study_level_null.csv", index=False)
    return by, summaries


atlas = atlas_rows()
atlas.to_csv(DAT / "f6_final_atlas.csv", index=False)
_, synthesis = study_level_synthesis()
print(atlas.to_string(index=False))
print(synthesis.to_string(index=False))

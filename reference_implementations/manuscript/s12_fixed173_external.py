# -*- coding: utf-8 -*-
"""Recompute the S12 Route-A exploratory aggregate with a fixed 173-gene denominator."""
from __future__ import annotations

import importlib.util
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
SCRIPT = WURZEL / "manuscript" / "reference_implementations" / "analyse_explorativ_aggregat.py"
OUT = WURZEL / "derived_data" / "manuscript"
SEED = 20260830
N_PERM = 10_000

spec = importlib.util.spec_from_file_location("s12_loader", SCRIPT)
loader = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(loader)

tau = pd.read_csv(loader.TAU, dtype={"gen": str, "symbol": str, "ri": int})
specs = [
    ("GSE185951", loader.diff_gse185951, "ensg"),
    ("GSE161176", loader.diff_gse161176, "symbol"),
    ("GSE113253", loader.diff_gse113253, "symbol"),
    ("GSE210984", loader.diff_gse210984, "ensg"),
    ("GSE202080", loader.diff_gse202080, "symbol"),
    ("GSE37521", loader.diff_gse37521, "symbol"),
]

mats = {}
per_study = []
rng = np.random.default_rng(SEED)
study_null_rows = []
for name, fn, key_type in specs:
    diff, _ = fn()
    diff = diff.groupby(level=0, sort=False).mean(numeric_only=True)
    keys = tau["gen"] if key_type == "ensg" else tau["symbol"]
    module_keys = list(keys)
    mapping = tau.set_index("gen" if key_type == "ensg" else "symbol")["ri"]
    complete = diff.reindex(module_keys).notna().all(axis=1)
    fixed = diff.reindex(module_keys).fillna(0.0)
    mat = fixed.to_numpy(dtype=float)
    expected = mapping.reindex(module_keys).to_numpy(dtype=int)
    obs = int((np.sign(mat.mean(axis=1)) == expected).sum())
    share = obs / len(module_keys)
    mats[name] = (mat, expected)
    study_null = np.empty(N_PERM, dtype=float)
    for i in range(N_PERM):
        flips = rng.choice(np.array([-1.0, 1.0]), size=mat.shape[1])
        study_null[i] = np.mean(np.sign(mat @ flips / mat.shape[1]) == expected)
    study_mean = float(study_null.mean())
    study_sd = float(study_null.std(ddof=1))
    per_study.append({"dataset": name, "n_units": mat.shape[1],
                      "n_module": len(module_keys),
                      "n_measurable": int(complete.sum()), "n_matching": obs,
                      "share_match": share, "null_mean": study_mean,
                      "null_sd": study_sd, "z": (share - study_mean) / study_sd,
                      "p_perm_two_sided": float(np.mean(np.abs(study_null - study_mean) >= abs(share - study_mean))),
                      "p_perm_upper": float(np.mean(study_null >= share)),
                      "mde80": study_mean + 2.80 * study_sd,
                      "above_mde80": share > study_mean + 2.80 * study_sd})
    study_null_rows.extend({"dataset": name, "share_null": v} for v in study_null)

denom = sum(mat.shape[0] for mat, _ in mats.values())
obs_matches = sum(row["n_matching"] for row in per_study)
obs_share = obs_matches / denom

null = np.empty(N_PERM, dtype=float)
for i in range(N_PERM):
    matches = 0
    for mat, expected in mats.values():
        flips = rng.choice(np.array([-1.0, 1.0]), size=mat.shape[1])
        matches += int((np.sign(mat @ flips / mat.shape[1]) == expected).sum())
    null[i] = matches / denom
null_mean = float(null.mean())
null_sd = float(null.std(ddof=1))
z = (obs_share - null_mean) / null_sd
summary = dict(n_studies=len(mats), n_observations=denom, n_matching=obs_matches,
               pooled_share=obs_share, null_mean=null_mean, null_sd=null_sd,
               z=z, p_perm_two_sided=float(np.mean(np.abs(null - null_mean) >= abs(obs_share - null_mean))),
               p_perm_upper=float(np.mean(null >= obs_share)),
               mde80=null_mean + 2.80 * null_sd, seed=SEED, perm_rounds=N_PERM)
pd.DataFrame(per_study).to_csv(OUT / "f6_s12_fixed173_by_study.csv", index=False)
pd.DataFrame(study_null_rows).to_csv(OUT / "f6_s12_fixed173_study_null.csv", index=False)
pd.DataFrame([summary]).to_csv(OUT / "f6_s12_fixed173_summary.csv", index=False)
pd.DataFrame({"pooled_share_null": null}).to_csv(OUT / "f6_s12_fixed173_null.csv", index=False)
print(pd.DataFrame(per_study).to_string(index=False))
print(pd.DataFrame([summary]).to_string(index=False))

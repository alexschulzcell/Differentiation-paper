# -*- coding: utf-8 -*-
"""MDE80 reachability per S12 dataset/route.

Question: does the pre-fixed threshold (null + 2.80*null SD) lie above the
maximum reachable value for this statistic (share -> 1.0; |rho| -> 1.0)?
Descriptive statistics only, no analysis. Sources: s12_route_*_summary.csv.
"""
from pathlib import Path
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
OUT = WURZEL / "derived_data" / "manuscript"
RES = SITZUNGEN / "26_Orthogonal_S12" / "derived_data"

rows = []

a = pd.read_csv(RES / "s12_route_a_summary.csv")
for r in a.itertuples():
    rows.append(dict(route="A", dataset=r.dataset, effect_type="share",
                     effect=r.share_match, null_mean=r.perm_mean, null_sd=r.perm_sd,
                     mde80=r.mde80, max_reachable=1.0,
                     reachable_lt_max=bool(r.mde80 < 1.0),
                     signed_effect=r.share_match, p_perm=r.p_perm_two_sided))

b = pd.read_csv(RES / "s12_route_b_summary.csv")
for r in b.itertuples():
    rows.append(dict(route="B", dataset=r.dataset, effect_type="share",
                     effect=r.share_match, null_mean=r.perm_mean, null_sd=r.perm_sd,
                     mde80=r.mde80, max_reachable=1.0,
                     reachable_lt_max=bool(r.mde80 < 1.0),
                     signed_effect=r.share_match, p_perm=r.p_perm_two_sided))

c = pd.read_csv(RES / "s12_route_c_summary.csv")
for r in c.itertuples():
    rows.append(dict(route="C", dataset=r.dataset, effect_type="abs_rho",
                     effect=r.abs_rho, null_mean=r.perm_mean, null_sd=r.perm_sd,
                     mde80=r.mde80, max_reachable=1.0,
                     reachable_lt_max=bool(r.mde80 < 1.0),
                     signed_effect=r.rho, p_perm=r.p_perm_two_sided))

d = pd.read_csv(RES / "s12_route_d_summary.csv")
for r in d.itertuples():
    rows.append(dict(route="D", dataset=r.dataset, effect_type="share",
                     effect=r.share_match, null_mean=r.perm_mean, null_sd=r.perm_sd,
                     mde80=r.mde80, max_reachable=1.0,
                     reachable_lt_max=bool(r.mde80 < 1.0),
                     signed_effect=r.share_match, p_perm=r.p_perm_two_sided))

df = pd.DataFrame(rows)
df.to_csv(OUT / "f6_mde80_erreichbarkeit.csv", index=False)
n_unreach = int((~df.reachable_lt_max).sum())
print(df.to_string(index=False))
print(f"\nunerreichbare Schwellen (MDE80 > 1.0): {n_unreach} von {len(df)}")

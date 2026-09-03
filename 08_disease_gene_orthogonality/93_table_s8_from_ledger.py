# -*- coding: utf-8 -*-
"""Build Supplementary Table S8 mechanically from the S12 candidate ledger."""
from pathlib import Path
import csv
import pandas as pd

# --- path parameters (2026-08-23) ------------------------------------------
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
S12 = SITZUNGEN / "26_Orthogonal_S12"
OUT = WURZEL / "derived_data" / "reference_tables" / "S8_orthogonale_triangulation.csv"

candidate_path = S12 / "data_raw" / "s12_kandidaten.csv"
with candidate_path.open(encoding="utf-8", newline="") as fh:
    reader = csv.reader(fh)
    header = next(reader)
    records = []
    statuses = {"ausgewertet", "ausgewertet_deskriptiv", "ausgeschlossen"}
    for fields in reader:
        status_idx = next((i for i, value in enumerate(fields) if value in statuses), None)
        if status_idx is None or status_idx < 5:
            continue
        art_idx = status_idx - 2
        unit_idx = status_idx - 1
        records.append([
            fields[0], fields[1], ", ".join(fields[2:art_idx]), fields[art_idx],
            fields[unit_idx], fields[status_idx], ", ".join(fields[status_idx + 1:]),
        ])
candidate = pd.DataFrame(records, columns=header)
summary_files = {
    "A": "s12_route_a_summary.csv",
    "B": "s12_route_b_summary.csv",
    "C": "s12_route_c_summary.csv",
    "D": "s12_route_d_summary.csv",
}

rows = []
for route, filename in summary_files.items():
    path = S12 / "derived_data" / filename
    if not path.exists():
        continue
    summary = pd.read_csv(path)
    for row in summary.to_dict("records"):
        if route == "C":
            effect = row["abs_rho"]
            effect_type = "abs_rho"
        else:
            effect = row["share_match"]
            effect_type = "share"
        rows.append({
            "route": route,
            "gse": row["dataset"],
            "n_units": row["n_biological_units"],
            "effect_type": effect_type,
            "effect": effect,
            "null_mean": row["perm_mean"],
            "null_sd": row["perm_sd"],
            "mde80": row["mde80"],
            "p_perm": row["p_perm_two_sided"],
            "z": row["z"],
            "reachable_lt_1.0": row["mde80"] < 1.0,
            "final_verdict": "nicht konfirmatorisch wegen n"
            if row["n_biological_units"] < (5 if route == "C" else 3)
            else "unter Schwelle",
        })

stats = pd.DataFrame(rows)
if stats.empty:
    raise RuntimeError("No S12 summary CSVs found")

candidate = candidate.rename(columns={"titel": "kandidat"})
out = candidate[["route", "gse", "kandidat", "art", "einheiten_je_zeitpunkt", "status"]].copy()
out = out.rename(columns={"einheiten_je_zeitpunkt": "candidate_units"})
out = out.merge(stats, on=["route", "gse"], how="left")
out.loc[out["status"] == "ausgeschlossen", "final_verdict"] = "ausgeschlossen"
out.loc[
    (out["status"] == "ausgewertet_deskriptiv") & out["final_verdict"].isna(),
    "final_verdict",
] = "nicht konfirmatorisch wegen n"
out = out[[
    "route", "gse", "kandidat", "art", "candidate_units", "status", "n_units",
    "effect_type", "effect", "null_mean", "null_sd", "mde80", "p_perm", "z",
    "reachable_lt_1.0", "final_verdict",
]]
out.to_csv(OUT, index=False)
print(out.to_string(index=False))
print(f"Wrote {OUT}")

# -*- coding: utf-8 -*-
"""
reproduce.py -- run the whole analysis from this repository, in order.

Everything this runner touches reads only files that are committed here. No
raw data, no network. It re-derives every number of the paper, redraws every
figure, and then verifies the manuscript against what it just computed.

    python reproduce.py                # the whole chain (~10 min)
    python reproduce.py --list         # show the steps, run nothing
    python reproduce.py --only checks  # one group only
    python reproduce.py --from figures # from that group onward
    python reproduce.py --package      # also rebuild submission/ (needs pandoc + Word)

Groups, in dependency order:

    analysis   the calibration, the decomposition, held-out and external
               validation, and the in vivo growth plate
    figures    one CSV per panel, then F1-F6, S1-S9 and the graphical abstract
    checks     every number of the manuscript against its panel file, plus
               the reference and language checks

The last group is the self-test: if `checks` passes, the manuscript is
consistent with the numbers this runner just produced.

Not run here, because they need input this repository does not carry:
  * everything marked "needs raw data" in the stage READMEs -- the roughly
    98 GB of public raw data under `data_raw/` (see 00_setup.md)
  * `06_orthogonal_layers/60_gene_sets_build.R` and `61_gene_set_enrichment.R`
    -- they need MSigDB and `org.Hs.eg.db`; their outputs are frozen under
    `results/` and the figure steps read those
  * `07_in_vivo_growth_plate/11_fetal_atlas_pseudobulk_store.py` -- needs the
    7.6 GB limb atlas; its pseudobulk is frozen under `derived_data/`
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent

# (group, script, what it produces)
STEPS: list[tuple[str, str, str]] = [
    ("analysis", "03_lineage_calibration/10_calibration_18_datasets.py",
     "the lineage calibration: 2 of 18 data sets pass, 7 of 14 donor cells"),
    ("analysis", "03_lineage_calibration/11_calibration_sensitivity.py",
     "does the calibration verdict depend on its thresholds?"),
    ("analysis", "03_lineage_calibration/12_calibration_gene_space.py",
     "which gene space belongs under the calibration"),
    ("analysis", "04_programme_decomposition/10_decomposition_18_datasets.py",
     "the three-way decomposition: 8 confirmations, 2 other, 0 refutations"),
    ("analysis", "05_programme_validation/10_heldout_and_robustness.py",
     "leave-one-study-out re-derivation, matched nulls, gene dropout"),
    ("analysis", "05_programme_validation/11_external_differentiation_systems.py",
     "the locked programme on four independent cohorts: 3 of 4 above limit"),
    ("analysis", "07_in_vivo_growth_plate/12_fetal_donor_trend_test.py",
     "the fetal growth-plate trend, stratified by donor"),
    ("analysis", "07_in_vivo_growth_plate/13_fetal_gene_decomposition.py",
     "is the in vivo trend broadly carried, or by few genes?"),
    ("analysis", "07_in_vivo_growth_plate/14_hypertrophic_zone_sensitivity.py",
     "does the trend hinge on the terminal hypertrophic point?"),
    ("analysis", "07_in_vivo_growth_plate/21_postnatal_growth_plate_test.py",
     "the postnatal anchor -- fails its own calibration, carries no verdict"),

    ("figures", "09_figures/10_panel_data_main.py",
     "one CSV per main-figure panel, under figures/data/"),
    ("figures", "09_figures/11_panel_data_supplement.py",
     "the supplement panels and Tables S1-S14"),
    ("figures", "09_figures/12_panel_data_second_cohort.py",
     "the panels of supplementary figure S9"),
    ("figures", "09_figures/20_figures_main.R",
     "F1-F6 as PDF and PNG at 600 dpi"),
    ("figures", "09_figures/21_figures_supplement.R",
     "S1-S9 as PDF and PNG at 600 dpi"),
    ("figures", "09_figures/30_graphical_abstract.py",
     "the graphical abstract at 300 dpi"),

    ("checks", "10_manuscript_checks/10_check_numbers.py",
     "every load-bearing number of the manuscript against its panel file"),
    ("checks", "10_manuscript_checks/11_check_references.py",
     "references and citations, in both directions"),
    ("checks", "10_manuscript_checks/12_check_language.py",
     "the language rules of the material"),
]

PACKAGE = ("10_manuscript_checks/21_build_submission.py",
           "the submission package under submission/")

GROUPS = ["analysis", "figures", "checks"]


def runner(script: str) -> list[str]:
    return ([sys.executable] if script.endswith(".py") else ["Rscript"]) + \
        [str(ROOT / script)]


def run(script: str, what: str, n: int, total: int) -> bool:
    print(f"\n[{n}/{total}] {script}\n        {what}", flush=True)
    t = time.time()
    r = subprocess.run(runner(script), cwd=ROOT,
                       capture_output=True, text=True)
    if r.returncode == 0:
        print(f"        ok ({time.time() - t:.0f} s)", flush=True)
        return True
    print(f"        FAILED after {time.time() - t:.0f} s", flush=True)
    tail = (r.stdout + r.stderr).strip().splitlines()[-25:]
    for ln in tail:
        print("        | " + ln, flush=True)
    return False


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Reproduce the analysis from this repository alone.")
    ap.add_argument("--list", action="store_true",
                    help="show the steps and exit")
    ap.add_argument("--only", choices=GROUPS, help="run one group only")
    ap.add_argument("--from", dest="start", choices=GROUPS,
                    help="run from this group onward")
    ap.add_argument("--package", action="store_true",
                    help="also rebuild submission/ (needs pandoc and Word)")
    a = ap.parse_args()

    steps = list(STEPS)
    if a.only:
        steps = [s for s in steps if s[0] == a.only]
    elif a.start:
        i = GROUPS.index(a.start)
        keep = set(GROUPS[i:])
        steps = [s for s in steps if s[0] in keep]

    if a.package and (not a.only or a.only == "checks"):
        steps.append(("checks", *PACKAGE))

    if a.list:
        group = None
        for g, s, what in steps:
            if g != group:
                print(f"\n{g}")
                group = g
            print(f"  {s}\n      {what}")
        return 0

    print("=" * 78)
    print("Reproducing from this repository alone -- no raw data, no network.")
    print("=" * 78)

    failed = []
    for i, (_, s, what) in enumerate(steps, 1):
        if not run(s, what, i, len(steps)):
            failed.append(s)

    print("\n" + "=" * 78)
    if failed:
        print(f"{len(failed)} of {len(steps)} steps FAILED:")
        for s in failed:
            print("  " + s)
        print("=" * 78)
        return 1
    print(f"All {len(steps)} steps passed.")
    if any(g == "checks" for g, _, _ in steps):
        print("The manuscript is consistent with the numbers just computed.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

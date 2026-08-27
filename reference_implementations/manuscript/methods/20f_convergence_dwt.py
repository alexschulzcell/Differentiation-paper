# -*- coding: utf-8 -*-
"""Reproduce the exploratory 173-gene dWT module.

This is the public, standalone equivalent of section (3) in the archived
``20i_dexamethasone.R``.  It does not define a new module: it reconstructs the
existing exploratory selection and checks it against the frozen S5 table.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = ROOT / "derived_data" / "reference_tables" / "20d_dWT_matrix.csv.gz"
DEFAULT_MODULE = ROOT / "derived_data" / "reference_tables" / "S5_konvergente_gene.csv"


def build_module(matrix_path: Path) -> tuple[pd.DataFrame, int]:
    """Apply the archived dWT convergence rule to the long matrix."""
    long = pd.read_csv(matrix_path, compression="gzip", dtype={"gen": str})
    required = {"gen", "punkt", "dWT"}
    missing = required.difference(long.columns)
    if missing:
        raise ValueError(f"matrix is missing columns: {sorted(missing)}")
    if long[["gen", "punkt"]].duplicated().any():
        raise ValueError("matrix contains duplicate gene/point rows")

    wide = long.pivot(index="gen", columns="punkt", values="dWT")
    wide = wide.reindex(columns=sorted(wide.columns))
    if list(wide.columns) != list(range(1, 19)):
        raise ValueError(f"expected points 1..18, found {list(wide.columns)}")

    # This order is material: the archived analysis first defines the
    # 14-of-18 universe and only then computes the per-point medians.
    available = wide.notna().sum(axis=1)
    universe = wide.loc[available >= 14]
    medians = universe.median(axis=0)
    signs = np.sign(universe.subtract(medians, axis="columns")).fillna(0)

    positive = (signs > 0).sum(axis=1)
    n = (signs != 0).sum(axis=1)
    v = np.maximum(positive, n - positive)
    ri = np.where(positive >= n - positive, 1, -1)

    result = pd.DataFrame(
        {
            "gen": universe.index.astype(str),
            "n": n.to_numpy(dtype=int),
            "v": v.to_numpy(dtype=int),
            "ri": ri.astype(int),
            "med": universe.median(axis=1).to_numpy(dtype=float),
        }
    )
    result = result.loc[result["n"] > 0]
    result = result.loc[result["v"] / result["n"] >= 0.90]
    return result.sort_values("gen").reset_index(drop=True), len(universe)


def verify(result: pd.DataFrame, module_path: Path, universe_size: int) -> None:
    """Verify the derived genes and directions against the frozen S5 table."""
    frozen = pd.read_csv(module_path, dtype={"gen": str, "ri": int})
    required = {"gen", "n", "v", "ri"}
    missing = required.difference(frozen.columns)
    if missing:
        raise ValueError(f"frozen module is missing columns: {sorted(missing)}")

    derived = result.set_index("gen").sort_index()
    expected = frozen.set_index("gen").sort_index()
    if set(derived.index) != set(expected.index):
        raise AssertionError(
            "derived gene set differs from S5: "
            f"missing={len(set(expected.index) - set(derived.index))}, "
            f"extra={len(set(derived.index) - set(expected.index))}"
        )

    for column in ("n", "v", "ri"):
        if not derived[column].equals(expected[column].astype(derived[column].dtype)):
            raise AssertionError(f"derived column {column!r} differs from S5")

    if len(derived) != 173:
        raise AssertionError(f"expected 173 genes, derived {len(derived)}")
    if (derived["ri"] > 0).sum() != 129 or (derived["ri"] < 0).sum() != 44:
        raise AssertionError("expected direction split +129/-44")

    print(f"PASS: {len(derived)} genes reproduced from {universe_size} evaluable genes")
    print("PASS: gene identities, n, v and ri match S5_konvergente_gene.csv")
    print("PASS: direction split is 129 up / 44 down")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--module", type=Path, default=DEFAULT_MODULE)
    parser.add_argument(
        "--write",
        type=Path,
        help="optionally write the reconstructed table to this path",
    )
    args = parser.parse_args()

    result, universe_size = build_module(args.matrix)
    verify(result, args.module, universe_size)
    if args.write:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.write, index=False)
        print(f"wrote: {args.write}")


if __name__ == "__main__":
    main()

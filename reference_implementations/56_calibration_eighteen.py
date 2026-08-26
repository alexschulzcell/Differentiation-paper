# -*- coding: utf-8 -*-
"""
56_calibration_eighteen.py -- the built-in calibration on the EIGHTEEN
perturbation datasets.

The same rule as in phase M-D, only one resolution level higher: does the
**own** differentiation contrast of a dataset find the canonical lineage
markers of **its own** axis?

    _module.kontrast(dWT, satz_a, all other marker sets)
    satz_a = OSTEOGEN or CHONDROGEN depending on the arm, directed, p < 0.05

The rule has been fixed since `PRAEREG_M_D.md` §6 and is NOT changed here.
Only the object is new: not the 14 donor-resolved cells, but the 18 datasets
from which the module originates.

**No second implementation**: `dWT` per dataset comes unchanged from the
frozen tables `20d_gene_*.csv` (computed there with `kern()`), the marker
sets from `_marker.py`, the contrast from `_module.kontrast`.

-----------------------------------------------------------------------------
ADDENDUM 2026-08-24 -- the gene space, and the passing criterion
-----------------------------------------------------------------------------
Two things are corrected; both are justified in KONSISTENZ_PROTOKOLL §16 and
re-computed in `code/29_calibration_gene_space.py`.

1. **The marker set is NO LONGER filtered.** The retired version of this
   script was computed with the internal gene map, which only carries genes
   with `dWT_n >= 15` of 18 datasets. This filter is right for defining the
   module and wrong for a positive control per experiment: it preferentially
   removes the TERMINAL markers (SP7, BGLAP, IBSP, DMP1, SOST, MEPE; ACAN,
   COL2A1, COMP, PRG4; FABP4, ADIPOQ), because these are off in the naive
   arm and only become measurable where differentiation actually runs --
   while NAIV remains complete at 10/10. The filter thus acts one-sidedly
   against the half of the calibration that measures ARRIVAL, and the
   datasets do measure these genes (point 15 measures 18 of 18 osteogenic
   markers, the internal map lets 12 through). The decomposition in `29_`
   shows that the marker set alone moves the z values (median |dz| 0.378),
   the background does not (0.017).

2. **`bestanden` now follows the preregistered criterion z >= 2**
   (`PRAEREG_M_D.md` §6). The script previously checked `p < 0.05`. Both
   criteria give the same number here; the captions have always said z >= 2,
   and now the code says it too.

The old output is kept as `eichung_achtzehn_MIN_N15_ALT.csv`.

**Honesty note that belongs in every analysis:** the 18 datasets and the 14
phase-D cells are **not independent**. All six study units of phase D are
also inside the 18; the difference is resolution (study level vs donor
level), not material. Both counts are therefore reported separately.

Output: derived_data/M_kalibrierung/eichung_achtzehn.csv
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import ERGEBNISSE, kontrast  # noqa: E402

WURZEL = pathlib.Path(__file__).resolve().parents[1]
ANTRAG = WURZEL.parent
AUS = ERGEBNISSE / "M_kalibrierung"
AUS.mkdir(parents=True, exist_ok=True)
GENE20D = (ANTRAG / "backups" / "_backup_2026-08-19_vor_paperaufbau" /
           "20_Exploration" / "derived_data")

MARKERSAETZE = {"OSTEOGEN": OSTEOGEN, "ADIPOGEN": ADIPOGEN,
                "MYOGEN": MYOGEN, "CHONDROGEN": CHONDROGEN, "NAIV": NAIV}
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def symbolkarte() -> dict:
    return gencode_karte()

def gencode_karte() -> dict:
    """Ensembl -> symbol from the Gencode reference -- the same map that
    `54b_cells.py` uses in phase M-D. The project's internal gene map covers
    only about 11 500 genes and leaves only 5 of the 12 chondrogenic markers;
    with it, calibration would be underdetermined on the chondrogenic axis.
    This choice is a matter of coverage, not of the result, and is therefore
    recorded here."""
    import gzip
    import re
    from _module import DATEN as _D
    for p in sorted((_D / "_referenz").glob("*.gtf*")):
        m = {}
        op = gzip.open if p.suffix == ".gz" else open
        with op(p, "rt", encoding="utf-8", errors="replace") as f:
            for ln in f:
                if ln[0] == "#" or "	gene	" not in ln:
                    continue
                g = re.search(r'gene_id "([^".]+)', ln)
                s = re.search(r'gene_name "([^"]+)', ln)
                if g and s:
                    m[g.group(1)] = s.group(1)
        if m:
            return m
    raise RuntimeError("no Gencode reference found")



def main() -> None:
    log("=" * 78)
    log("The built-in calibration on the 18 perturbation datasets")
    log("Rule unchanged from PRAEREG_M_D.md §6")
    log("=" * 78)

    karte = symbolkarte()
    KO = pd.read_csv(WURZEL / "derived_data" / "manuscript" / "f1_kohorte.csv")
    arm = dict(zip(KO.punkt, KO.arm))
    name = dict(zip(KO.punkt, KO.datensatz))

    zeilen = []
    for f in sorted(GENE20D.glob("20d_gene_*.csv")):
        G = pd.read_csv(f)
        p = int(G.punkt.iloc[0])
        G["symbol"] = [karte.get(str(g).split(".")[0]) for g in G.gen]
        G = G[G.symbol.notna() & G.dWT.notna()]
        dwt = G.groupby("symbol").dWT.median()

        a = arm.get(p, "")
        satz_a = "OSTEOGEN" if a == "osteogen" else "CHONDROGEN"
        setzt_a = MARKERSAETZE[satz_a]
        setzt_b = [g for n, s in MARKERSAETZE.items() if n != satz_a for g in s]
        r = kontrast(dwt, setzt_a, setzt_b)
        # Preregistered criterion: z >= 2 against the own null
        # (PRAEREG_M_D.md §6). See addendum in the header.
        best = bool(r.get("status") == "ok" and r.get("z", 0) >= 2
                    and r["kontrast"] > 0)
        r.update(punkt=p, datensatz=name.get(p, ""), arm=a, satz_a=satz_a,
                 n_gene_messbar=len(dwt), bestanden=best)
        zeilen.append(r)
        log("  %2d  %-28s %-11s %-11s contrast %+7.3f | z %+6.2f | p %8.4g "
            "-> %s"
            % (p, str(name.get(p, ""))[:28], a, satz_a.lower(),
               r.get("kontrast", np.nan), r.get("z", np.nan),
               r.get("p", np.nan),
               "PASSED" if best else "FAILED"))

    T = pd.DataFrame(zeilen)
    T.to_csv(AUS / "eichung_achtzehn.csv", index=False)

    n_ok = int(T.bestanden.sum())
    log("")
    log("-" * 78)
    log("STUDY LEVEL (18 perturbation datasets): %d of %d passed, "
        "%d failed (%.0f %%)"
        % (n_ok, len(T), len(T) - n_ok, 100 * (len(T) - n_ok) / len(T)))
    for a in sorted(T.arm.dropna().unique()):
        s = T[T.arm == a]
        log("   of which %-11s %d of %d passed" % (a, int(s.bestanden.sum()),
                                                    len(s)))

    # ---- joint view with phase M-D, explicitly NOT as independent ----------
    D = pd.read_csv(ERGEBNISSE / "M_donoren" / "eichung.csv")
    n_d_ok = int(D.bestanden.sum())
    log("")
    log("DONOR LEVEL (14 phase-D cells): %d of %d passed, %d "
        "failed (%.0f %%)"
        % (n_d_ok, len(D), len(D) - n_d_ok, 100 * (len(D) - n_d_ok) / len(D)))
    log("")
    log("TOTAL: %d of %d calibration attempts failed (%.0f %%)."
        % ((len(T) - n_ok) + (len(D) - n_d_ok), len(T) + len(D),
           100 * ((len(T) - n_ok) + (len(D) - n_d_ok)) / (len(T) + len(D))))
    log("WARNING: these are NOT %d independent experiments. All six"
        % (len(T) + len(D)))
    log("study units of phase D are also inside the 18; the two counts")
    log("differ in RESOLUTION, not in material.")
    log("The robust headline count is the study level: %d of %d."
        % (len(T) - n_ok, len(T)))
    log("-" * 78)

    (AUS / "56_log.txt").write_text("\n".join(LOG), encoding="utf-8")
    log("-> %s" % (AUS / "eichung_achtzehn.csv"))


if __name__ == "__main__":
    main()

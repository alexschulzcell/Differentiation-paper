# -*- coding: utf-8 -*-
"""
10_dnam_GSE33896.py -- orthogonal layer A: DNA methylation.

Dataset
-------
GSE33896, Illumina HumanMethylation27 (GPL8490), 27 578 CpG.
Donor-matched series from human adipose stem cells (hASC):

    hASC (naive)  ->  in vitro induced OSTEOCYTES   donor 1, 2, 3
    hASC (naive)  ->  in vitro induced MYOCYTES     donor 1, 2, 3

Both axes start in the same naive cell of the same donor. The myogenic
axis is therefore the **lineage control** for the osteogenic one: it asks
whether any methylation signature of the module characterizes osteogenic
differentiation or mere differentiation as such.

Cell lines (MG-63, TE 32.T, RD) and the non-donor-matched primary samples
are excluded -- cell lines are an exclusion criterion throughout the
project.

Hypothesis and sign
-------------------
Promoter CpG methylation and transcription act in opposite directions. A
gene with `ri = +1` (upward in the module) should therefore be methylated
LOWER in the differentiated state. The expected sign of the beta difference
is thus `-ri`. This sign rule is fixed before any computation and is not
reversed.

The unit of analysis is the DONOR (three paired donors), not the sample and
not the probe.
"""
from __future__ import annotations

import gzip
import io
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from _module import DATEN, ERGEBNISSE, MODUL, konkordanz, wilson  # noqa: E402

AUS = ERGEBNISSE / "A_dnam"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------- manifest
def lade_manifest() -> pd.DataFrame:
    p = DATEN / "GPL8490" / "GPL8490_HumanMethylation27_270596_v.1.2.csv.gz"
    roh = gzip.open(p, "rt", encoding="latin-1").read()
    kopf = roh.index("IlmnID,Name,")
    tab = pd.read_csv(io.StringIO(roh[kopf:]), low_memory=False)
    tab = tab[["Name", "Symbol", "Distance_to_TSS", "CPG_ISLAND", "Chr", "MapInfo"]]
    tab = tab.rename(columns={"Name": "sonde", "Symbol": "symbol",
                              "Distance_to_TSS": "abstand_tss",
                              "CPG_ISLAND": "cpg_insel"})
    tab["symbol"] = tab["symbol"].astype(str).str.strip()
    tab = tab[tab.symbol.notna() & (tab.symbol != "") & (tab.symbol != "nan")]
    return tab


# ------------------------------------------------------------ beta matrix
def lade_beta() -> tuple[pd.DataFrame, pd.DataFrame]:
    p = DATEN / "GSE33896" / "GSE33896_series_matrix.txt.gz"
    zeilen = gzip.open(p, "rt", encoding="utf-8", errors="replace").read().split("\n")
    gsm = titel = quelle = None
    start = None
    for i, z in enumerate(zeilen):
        if z.startswith("!Sample_geo_accession"):
            gsm = [x.strip('"') for x in z.split("\t")[1:]]
        elif z.startswith("!Sample_source_name_ch1"):
            quelle = [x.strip('"') for x in z.split("\t")[1:]]
        elif z.startswith("!Sample_title"):
            titel = [x.strip('"') for x in z.split("\t")[1:]]
        elif z.startswith("!series_matrix_table_begin"):
            start = i + 1
        elif z.startswith("!series_matrix_table_end"):
            ende = i
            break
    tab = pd.read_csv(io.StringIO("\n".join(zeilen[start:ende])), sep="\t",
                      index_col=0, low_memory=False)
    tab.index.name = "sonde"

    meta = pd.DataFrame({"gsm": gsm, "titel": titel, "quelle": quelle})
    def donor(t: str):
        m = re.search(r"donor (\d+)", t)
        return int(m.group(1)) if m else None
    meta["donor"] = meta.titel.map(donor)
    def arm(q: str, t: str):
        if q.startswith("hASC"):
            return "naiv"
        if q.startswith("DIF.O"):
            return "osteo"
        if q.startswith("DIF.M"):
            return "myo"
        return "ausgeschlossen"
    meta["arm"] = [arm(q, t) for q, t in zip(meta.quelle, meta.titel)]
    return tab, meta


def main() -> None:
    log("=" * 78)
    log("Layer A -- DNA methylation, GSE33896 (27K, hASC -> osteocyte / myocyte)")
    log("=" * 78)

    man = lade_manifest()
    beta, meta = lade_beta()
    log("Probes in the series matrix: %d | samples: %d" % beta.shape)

    behalte = meta[meta.arm.isin(["naiv", "osteo", "myo"]) & meta.donor.isin([1, 2, 3])]
    log("\nSamples after excluding cell lines and unpaired primary samples:")
    for _, r in behalte.iterrows():
        log("  %-12s donor %d  %-6s  %s" % (r.gsm, r.donor, r.arm, r.quelle))
    ausg = meta[~meta.gsm.isin(behalte.gsm)]
    log("  excluded (%d): %s" % (len(ausg), ", ".join(ausg.quelle)))

    beta = beta[behalte.gsm.tolist()].apply(pd.to_numeric, errors="coerce")
    log("\nBeta values: %d probes x %d samples, %.2f %% missing"
        % (beta.shape[0], beta.shape[1], 100 * beta.isna().mean().mean()))

    # ---------------------------------------------------- probes -> genes
    man = man[man.sonde.isin(beta.index)]
    log("Probes with gene annotation: %d over %d symbols"
        % (len(man), man.symbol.nunique()))

    # Promoter-proximal probes: the 27K array is a promoter array, but its
    # distances scatter. Fixed BEFORE analysis: |distance to TSS| <= 1500 bp.
    FENSTER = 1500
    man["abstand_tss"] = pd.to_numeric(man.abstand_tss, errors="coerce")
    manp = man[man.abstand_tss.abs() <= FENSTER]
    log("of these promoter-proximal (|d(TSS)| <= %d bp): %d probes over %d symbols"
        % (FENSTER, len(manp), manp.symbol.nunique()))

    # ------------------------------------------------------- coverage
    sym_module = set(MODUL.symbol)
    abgedeckt = sorted(sym_module & set(manp.symbol))
    log("\n--- Coverage of the fixed 173-gene module on the 27K array -------")
    log("covered: %d of 173 (%.1f %%)" % (len(abgedeckt), 100 * len(abgedeckt) / 173))
    if len(abgedeckt) / 173 < 0.50:
        log("WARNING: below 50 %% -- the layer carries the module test only partially.")

    # ------------------------------------------- donor-paired difference
    def achse(arm: str) -> tuple[pd.Series, pd.DataFrame]:
        """Gene-weighted beta difference `differentiated minus naive`, donor-paired."""
        proben = {}
        for d in (1, 2, 3):
            n = behalte[(behalte.donor == d) & (behalte.arm == "naiv")].gsm.tolist()
            t = behalte[(behalte.donor == d) & (behalte.arm == arm)].gsm.tolist()
            assert len(n) == 1 and len(t) == 1, (d, arm, n, t)
            proben[d] = (n[0], t[0])
        # per donor the probe difference, then average over donors
        dif = pd.DataFrame({d: beta[t] - beta[n] for d, (n, t) in proben.items()})
        dif.index = beta.index
        sonde = dif.mean(axis=1)
        vz_konsistenz = (np.sign(dif).abs().sum(axis=1) == 3) & \
                        (np.sign(dif).sum(axis=1).abs() == 3)
        # probe -> gene: mean over the promoter-proximal probes of a gene
        z = manp.assign(delta=manp.sonde.map(sonde),
                        konsistent=manp.sonde.map(vz_konsistenz))
        gen = z.groupby("symbol").agg(delta=("delta", "mean"),
                                      n_sonden=("delta", "size"),
                                      donor_konsistent=("konsistent", "mean"))
        return gen.delta.dropna(), gen

    ergebnisse = []
    genlisten = {}
    for arm, name in (("osteo", "osteogen (Zielachse)"),
                      ("myo", "myogen (Linienkontrolle)")):
        delta, gentab = achse(arm)
        genlisten[arm] = gentab
        log("\n--- Axis: %s -------------------------------------------" % name)
        log("Genes with promoter-proximal measurement: %d" % len(delta))
        log("global drift (median delta beta over all genes): %+.4f" % delta.median())

        d_mod = delta.reindex(abgedeckt).dropna()
        erw = pd.Series({s: -int(MODUL.loc[MODUL.symbol == s, "ri"].iloc[0])
                         for s in d_mod.index})
        res = konkordanz(d_mod, erw, hintergrund=delta)
        k = int((np.sign(d_mod.values) == erw.reindex(d_mod.index).values).sum())
        lo, hi = wilson(k, len(d_mod))
        log("Module genes measurable: %d | concordant: %d (%.3f, Wilson %.3f-%.3f)"
            % (len(d_mod), k, res["konkordanz"], lo, hi))
        log("Null (background draw, same sign set): %.3f +- %.3f"
            % (res["konkordanz_null"], res["konkordanz_null_sd"]))
        log("z = %+.2f | p = %.4g | MDE80 threshold = %.3f"
            % (res["konkordanz_z"], res["konkordanz_p"], res["konkordanz_mde80"]))
        log("continuous (signed rank): %+.4f, z %+.2f, p %.4g"
            % (res["rang"], res["rang_z"], res["rang_p"]))
        log("threshold reached: %s"
            % ("YES" if res["konkordanz"] >= res["konkordanz_mde80"] else "no"))
        ergebnisse.append(dict(ebene="DNAm", datensatz="GSE33896", achse=arm,
                               achse_name=name, k=k, wilson_lo=lo, wilson_hi=hi,
                               **{a: b for a, b in res.items() if a != "status"}))

    E = pd.DataFrame(ergebnisse)
    E.to_csv(AUS / "A_dnam_GSE33896_modultest.csv", index=False)

    # gene table for the figure
    G = pd.DataFrame({"symbol": abgedeckt})
    G["ri"] = G.symbol.map(dict(zip(MODUL.symbol, MODUL.ri)))
    for arm in ("osteo", "myo"):
        G["delta_" + arm] = G.symbol.map(genlisten[arm].delta)
        G["n_sonden_" + arm] = G.symbol.map(genlisten[arm].n_sonden)
        G["donorkons_" + arm] = G.symbol.map(genlisten[arm].donor_konsistent)
    G["erwartet"] = -G.ri
    G.to_csv(AUS / "A_dnam_GSE33896_modulgene.csv", index=False)

    # background distribution for the figure
    for arm in ("osteo", "myo"):
        genlisten[arm].reset_index().to_csv(
            AUS / ("A_dnam_GSE33896_alle_gene_%s.csv" % arm), index=False)

    # --------------------------------------- direct axis comparison
    log("\n--- Osteogenic vs myogenic, same donors, same genes --------------")
    a = G["delta_osteo"] * G.erwartet
    b = G["delta_myo"] * G.erwartet
    ok = a.notna() & b.notna()
    w = stats.wilcoxon(a[ok], b[ok])
    log("mean expectation-directed delta beta: osteo %+.4f vs myo %+.4f"
        % (a[ok].mean(), b[ok].mean()))
    log("Wilcoxon paired over %d module genes: W %.0f, p %.4g"
        % (ok.sum(), w.statistic, w.pvalue))
    log("\nReading note: the myogenic axis is not a negative control in the")
    log("strict sense -- it shares with the osteogenic axis the exit from the")
    log("naive state. It separates only between 'differentiation as such'")
    log("and 'osteogenic differentiation'.")

    (AUS / "A_dnam_GSE33896_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

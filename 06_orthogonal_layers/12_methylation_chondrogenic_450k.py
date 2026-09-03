# -*- coding: utf-8 -*-
"""
12_methylation_chondrogenic_450k.py -- layer A, second dataset: 450K, CHONDROGENESIS.

GSE33896 covers the osteogenic half of the paper's axis on the thin 27K
array. GSE129266 adds the chondrogenic half on the 450K array (485 553
probes), human bone-marrow MSCs, day 0 vs day 14.

Sample situation (from the series matrix):
    donor 1: day 0 -> 4 samples, day 14 -> 3 samples
    donor 2: day 0 -> 2 samples, day 14 -> 1 sample
    donor 3: only day 0        donor 4: only day 14
Donors 1 and 2 are thus evaluable as pairs. Donors 3 and 4 are excluded --
an unpaired arm would be a pure donor difference. The unit of analysis is
the donor.

Beta from the non-normalized signals:
    beta = methylated / (methylated + unmethylated + 100)
The offset 100 is the Illumina convention and stabilizes probes with low
total intensity.

Probe -> gene via the Illumina manifest (`UCSC_RefGene_Name` /
`UCSC_RefGene_Group`); only promoter-proximal groups count
(TSS1500, TSS200, 5'UTR, 1stExon). Sign rule as on 27K: expected is `-ri`.
"""
from __future__ import annotations

import gzip
import io
import pathlib
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
from _marker import CHONDROGEN, NAIV  # noqa: E402
from _module import DATEN, ERGEBNISSE, MODUL, konkordanz, kontrast, wilson  # noqa: E402

AUS = ERGEBNISSE / "A_dnam"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []
PROMOTORGRUPPEN = {"TSS1500", "TSS200", "5'UTR", "1stExon"}


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


def sonde_gen() -> pd.DataFrame:
    """Probe -> gene symbol, only promoter-proximal assignments."""
    zwischen = AUS / "_450k_sonde_gen.csv"
    if zwischen.exists():
        return pd.read_csv(zwischen)
    p = DATEN / "GPL13534" / "manifest450k.csv.gz"
    zeilen = []
    with gzip.open(p, "rt", encoding="latin-1") as f:
        for z in f:
            if z.startswith("IlmnID,Name,"):
                kopf = z.rstrip("\n").split(",")
                break
        i_name = kopf.index("Name")
        i_gen = kopf.index("UCSC_RefGene_Name")
        i_grp = kopf.index("UCSC_RefGene_Group")
        for z in f:
            sp = z.rstrip("\n").split(",")
            if len(sp) <= i_grp or not sp[i_gen]:
                continue
            gene = sp[i_gen].split(";")
            grup = sp[i_grp].split(";")
            gefunden = {g for g, gr in zip(gene, grup) if gr in PROMOTORGRUPPEN}
            for g in gefunden:
                zeilen.append((sp[i_name], g))
    t = pd.DataFrame(zeilen, columns=["sonde", "symbol"]).drop_duplicates()
    t.to_csv(zwischen, index=False)
    return t


def lade_metadaten() -> pd.DataFrame:
    p = DATEN / "GSE129266" / "GSE129266_series_matrix.txt.gz"
    txt = gzip.open(p, "rt", encoding="utf-8", errors="replace").read().split("\n")
    felder = {}
    for z in txt:
        if z.startswith("!Sample_title"):
            felder["titel"] = [x.strip('"') for x in z.split("\t")[1:]]
        elif z.startswith("!Sample_geo_accession"):
            felder["gsm"] = [x.strip('"') for x in z.split("\t")[1:]]
        elif z.startswith("!Sample_characteristics_ch1"):
            w = [x.strip('"') for x in z.split("\t")[1:]]
            if w and w[0].startswith("donor:"):
                felder["donor"] = [int(x.split(":")[1]) for x in w]
            elif w and w[0].startswith("treatment:"):
                felder["arm"] = [x.split(": ")[1] for x in w]
    M = pd.DataFrame(felder)
    M["spalte"] = ["SAMPLE %d" % (i + 1) for i in range(len(M))]
    return M


def main() -> None:
    log("=" * 78)
    log("Layer A (450K) -- GSE129266, MSC chondrogenesis day 0 vs day 14")
    log("=" * 78)

    M = lade_metadaten()
    log("\nSamples according to the series matrix:")
    for _, r in M.iterrows():
        log("  %-10s %-22s donor %d  %s" % (r.gsm, r.spalte, r.donor, r.arm))

    gepaart = [d for d in sorted(M.donor.unique())
               if M[M.donor == d].arm.nunique() == 2]
    log("\nDonors with both arms: %s -- only these are analyzed."
        % ", ".join(map(str, gepaart)))
    log("excluded: donors %s (only one arm)"
        % ", ".join(str(d) for d in sorted(set(M.donor) - set(gepaart))))
    M = M[M.donor.isin(gepaart)]

    # ------------------------------------------------------------ beta
    p = DATEN / "GSE129266" / "GSE129266_non-normalized.txt.gz"
    T = pd.read_csv(p, sep="\t", index_col=0, low_memory=False)
    log("\nSignal table: %d probes x %d columns" % T.shape)
    B = {}
    for _, r in M.iterrows():
        u = pd.to_numeric(T["%s Unmethylated Signal" % r.spalte], errors="coerce")
        m = pd.to_numeric(T["%s Methylated signal" % r.spalte], errors="coerce")
        B[r.gsm] = m / (m + u + 100)
    B = pd.DataFrame(B)
    log("Beta matrix: %d probes x %d samples | median %.3f"
        % (B.shape[0], B.shape[1], np.nanmedian(B.values)))

    # ---------------------------------------------- donor-paired difference
    teil = []
    for d in gepaart:
        s0 = M[(M.donor == d) & (M.arm == "undifferentiated")].gsm.tolist()
        s1 = M[(M.donor == d) & (M.arm == "chondrogenic")].gsm.tolist()
        log("  donor %d: %d day-0 samples, %d day-14 samples" % (d, len(s0), len(s1)))
        teil.append(B[s1].mean(axis=1) - B[s0].mean(axis=1))
    DIF = pd.concat(teil, axis=1)
    DIF.columns = ["Donor%d" % d for d in gepaart]
    sonde_delta = DIF.mean(axis=1)
    donor_kons = (np.sign(DIF).abs().sum(axis=1) == len(gepaart)) & \
                 (np.sign(DIF).sum(axis=1).abs() == len(gepaart))
    log("\nProbes with the same sign in both donors: %d of %d (%.1f %%)"
        % (donor_kons.sum(), len(donor_kons), 100 * donor_kons.mean()))
    log("global drift (median delta beta): %+.4f" % sonde_delta.median())

    # ------------------------------------------------------ probe -> gene
    SG = sonde_gen()
    SG = SG[SG.sonde.isin(sonde_delta.index)]
    log("\npromoter-proximal probe-gene assignments: %d over %d symbols"
        % (len(SG), SG.symbol.nunique()))
    SG = SG.assign(delta=SG.sonde.map(sonde_delta),
                   kons=SG.sonde.map(donor_kons))
    GEN = SG.groupby("symbol").agg(delta=("delta", "mean"),
                                   n_sonden=("delta", "size"),
                                   donorkons=("kons", "mean")).dropna(subset=["delta"])
    log("Genes with promoter-proximal measurement: %d" % len(GEN))
    GEN.to_csv(AUS / "A_dnam450_GSE129266_gene.csv")

    # ------------------------------------------------------- calibration
    log("\n--- Calibration: chondrogenic minus naive markers ----------------")
    log("Expectation: chondrogenic markers lose promoter methylation,")
    log("i.e. NEGATIVE delta beta -- the contrast must come out negative.")
    r = kontrast(GEN.delta, CHONDROGEN, NAIV)
    if r.get("status") == "ok":
        # expected direction is negative -> flip the sign for the threshold
        ok = (-r["kontrast"]) >= (-r["null_mittel"] + 2.8 * r["null_sd"])
        log("contrast %+.4f | null %+.4f +- %.4f | z %+.2f | p %.4g | %s"
            % (r["kontrast"], r["null_mittel"], r["null_sd"], r["z"], r["p"],
               "PASSED" if ok else "failed"))
    else:
        ok = False
        log("calibration not computable: %s" % r.get("status"))
    pd.DataFrame([dict(datensatz="GSE129266", bestanden=ok,
                       **{k: v for k, v in r.items() if k != "status"})]).to_csv(
        AUS / "A_dnam450_GSE129266_eichung.csv", index=False)

    # ------------------------------------------------------- module test
    log("\n--- Module test: the fixed 173 genes, expected sign -ri ----------")
    sym_ri = dict(zip(MODUL.symbol, MODUL.ri))
    hg = GEN.delta.dropna()
    mod = hg[hg.index.isin(sym_ri)]
    erw = pd.Series({s: -sym_ri[s] for s in mod.index})
    log("Module genes measurable on 450K: %d of 173 (%.1f %%)"
        % (len(mod), 100 * len(mod) / 173))
    res = konkordanz(mod, erw, hintergrund=hg)
    k = int((np.sign(mod.values) == erw.reindex(mod.index).values).sum())
    lo, hi = wilson(k, len(mod))
    log("concordant %d -> C %.3f [%.3f-%.3f] | null %.3f+-%.3f | MDE80 %.3f"
        % (k, res["konkordanz"], lo, hi, res["konkordanz_null"],
           res["konkordanz_null_sd"], res["konkordanz_mde80"]))
    log("z %+.2f | p %.4g | continuous %+.4f (z %+.2f, p %.4g)"
        % (res["konkordanz_z"], res["konkordanz_p"], res["rang"],
           res["rang_z"], res["rang_p"]))
    log("threshold reached: %s"
        % ("YES" if res["konkordanz"] >= res["konkordanz_mde80"] else "no"))
    pd.DataFrame([dict(ebene="DNAm450", datensatz="GSE129266", achse="chondrogen",
                       geeicht=ok, k=k, wilson_lo=lo, wilson_hi=hi,
                       **{x: y for x, y in res.items() if x != "status"})]).to_csv(
        AUS / "A_dnam450_GSE129266_modultest.csv", index=False)

    (AUS / "A_dnam450_GSE129266_log.txt").write_text("\n".join(LOG) + "\n",
                                                     encoding="utf-8")
    print("\nwritten to", AUS)


if __name__ == "__main__":
    main()

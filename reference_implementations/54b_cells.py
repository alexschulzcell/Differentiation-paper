# -*- coding: utf-8 -*-
"""
54b_cells.py -- phase M-D: build the cells, check for LiCl, calibrate.

Preregistration: `preregistrations/PRAEREG_M_D.md` (2026-08-22, before the
first download) including addendum 1 (before the first statistic).

A **cell** is a donor x differentiation axis x study with a complete 2 x 2.
Per cell are formed (§3 of the preregistration, identical to `03_metric.R`,
function `kern`):

    dWT_p = z(control_p, diff) - z(control_p, naive)
    iv_p  = [z(lesion_p, diff) - z(lesion_p, naive)] - dWT_p

The z-standardization runs per gene over the included samples of THE SAME
study. If an arm has multiple samples (replicates or several
differentiation time points), they are averaged -- the project convention
from `03a_data.R`: no time point is selected.

This script computes NO statistics of the ladder. It builds the vectors,
records the LiCl check and runs the built-in calibration (`_module.kontrast`
on `dWT` against the lineage markers of its own axis).

Output: derived_data/M_donoren/{zellen.pkl, zellen_sichtung.csv,
         eichung.csv, licl_pruefung.csv, gene_level.csv.gz, 54b_log.txt}
"""
from __future__ import annotations

import gzip
import pathlib
import pickle
import re
import sys
import tarfile

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import DATEN, ERGEBNISSE, MODUL, kontrast  # noqa: E402

AUS = ERGEBNISSE / "M_donoren"
AUS.mkdir(parents=True, exist_ok=True)
WURZEL = pathlib.Path(__file__).resolve().parents[1]
ANTRAG = WURZEL.parent
ALT19 = ANTRAG / "_backup_2026-08-19_vor_paperaufbau" / "03_Metrik_Elf_Punkte" / "data_raw"

MARKERSAETZE = {"OSTEOGEN": OSTEOGEN, "ADIPOGEN": ADIPOGEN,
                "MYOGEN": MYOGEN, "CHONDROGEN": CHONDROGEN, "NAIV": NAIV}
LOG: list[str] = []


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def ens_symbol_karte() -> dict:
    """Ensembl -> symbol from the project's Gencode reference."""
    for p in sorted((DATEN / "_referenz").glob("*.gtf*")):
        m = {}
        op = gzip.open if p.suffix == ".gz" else open
        with op(p, "rt", encoding="utf-8", errors="replace") as f:
            for ln in f:
                if ln[0] == "#" or "\tgene\t" not in ln:
                    continue
                g = re.search(r'gene_id "([^".]+)', ln)
                s = re.search(r'gene_name "([^"]+)', ln)
                if g and s:
                    m[g.group(1)] = s.group(1)
        if m:
            return m
    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    return dict(zip(K.ensembl, K.symbol))


KARTE = ens_symbol_karte()


def auf_symbole(X: pd.DataFrame, ist_ensembl: bool = True) -> pd.DataFrame:
    if ist_ensembl:
        sym = pd.Series([KARTE.get(str(i).split(".")[0]) for i in X.index],
                        index=X.index)
    else:
        sym = pd.Series(X.index, index=X.index)
    X = X.copy()
    X["symbol"] = sym.values
    X = X[X.symbol.notna() & (X.symbol.astype(str) != "")]
    return X.groupby("symbol").median()


def zmat(X: pd.DataFrame) -> pd.DataFrame:
    """z per gene over the included samples of the same study."""
    sd = X.std(axis=1, ddof=1).replace(0, np.nan)
    return X.sub(X.mean(axis=1), axis=0).div(sd, axis=0).dropna(how="all")


def zelle(Z: pd.DataFrame, k_naiv, k_diff, l_naiv, l_diff) -> tuple:
    """dWT, iv and baseline (naive control arm) of a cell."""
    m = lambda c: Z[list(c)].mean(axis=1)              # noqa: E731
    dwt = m(k_diff) - m(k_naiv)
    iv = (m(l_diff) - m(l_naiv)) - dwt
    return dwt, iv, m(k_naiv)


def licl_pruefung(name: str, proben: list) -> dict:
    """Preregistration §10: ALWAYS check, ALWAYS record."""
    treffer = [p for p in proben
               if re.search(r"licl|lithium", str(p), flags=re.I)]
    return {"studie": name, "n_proben": len(proben),
            "licl_treffer": len(treffer),
            "treffer": "; ".join(map(str, treffer)),
            "urteil": "L1 -- ausgeschlossen" if treffer else "keine LiCl-Probe"}


# ---------------------------------------------------------------------------
# The studies.  Each returns:
#   Z       DataFrame symbol x sample, z-standardized
#   zellen  list of dicts with donor, axis, arm assignment, e2
#   proben  list of all sample names (for the LiCl check)
# ---------------------------------------------------------------------------
def lade_GSE218101():
    """MPS VI (ARSB), chondrogenic. 4 patient lines, each empty vector (EV,
    lesion) against gene-corrected (GE, control); iPS = naive, D14 = diff.
    The two healthy lines exist only at day 0 and do not belong in the
    isogenic 2x2 -- they are not used (as in `03a_data.R`)."""
    D = pd.read_csv(ALT19 / "GSE218101_CPM.csv.gz")
    D = D[~D.iloc[:, 0].duplicated()].set_index(D.columns[0])
    D = D.apply(pd.to_numeric, errors="coerce")
    proben = list(D.columns)
    sp = [c for c in proben if c.startswith("Line #")]
    X = np.log2(D[sp] + 1)
    Z = zmat(auf_symbole(X, ist_ensembl=False))
    zellen = []
    for k in (1, 2, 3, 4):
        f = lambda a, t: [c for c in sp if c.startswith(f"Line #{k}_{a}_{t}")]  # noqa: E731
        zellen.append(dict(spender=f"GSE218101_Line{k}", achse="chondrogen",
                           k_naiv=f("GE", "iPS"), k_diff=f("GE", "D14"),
                           l_naiv=f("EV", "iPS"), l_diff=f("EV", "D14"),
                           e2=True, laesion="ARSB-Patientenmutation, isogen korrigiert"))
    return Z, zellen, proben


def lade_GSE221128():
    """FOP (ACVR1), iMSC, chondrogenic. FOP (lesion) against resFOP (isogenic
    corrected, control), day 0 vs day 6. ex1-ex3 are experiments of the same
    line -- replicates, not a donor (rule 0.4)."""
    D = pd.read_csv(ALT19 / "GSE221128_TPM_allsamples.txt.gz", sep="\t")
    D = D.set_index("gene_id").drop(columns=["gene_name"])
    D = D.apply(pd.to_numeric, errors="coerce")
    proben = list(D.columns)
    Z = zmat(auf_symbole(np.log2(D + 1)))
    f = lambda k: [c for c in proben if k in c]        # noqa: E731
    zellen = [dict(spender="GSE221128_FOP", achse="chondrogen",
                   k_naiv=f("_R0_"), k_diff=f("_R6_"),
                   l_naiv=f("_F0_"), l_diff=f("_F6_"),
                   e2=True, laesion="ACVR1 R206H, isogen korrigiert (resFOP)")]
    return Z, zellen, proben


def lade_GSE245585():
    """RB1 (retinoblastoma), patient-derived MSC, osteogenic. WT1 (control)
    against MT1 (RB+/- patient line), day 0 vs day 7/14/21. By the project
    convention ALL differentiation time points enter as `diff`; no time
    point is selected."""
    ordner = ALT19 / "GSE245585"
    titel = {}
    with gzip.open(DATEN / "GSE245585" / "GSE245585_family.soft.gz", "rt",
                   errors="replace") as f:
        gsm = None
        for ln in f:
            if ln.startswith("^SAMPLE"):
                gsm = ln.split("=")[1].strip()
            elif ln.startswith("!Sample_title") and gsm:
                titel[gsm] = ln.split("=", 1)[1].strip()
    spalten = {}
    for p in sorted(ordner.glob("*.txt.gz")):
        g = p.name.split("_")[0]
        d = pd.read_csv(p, sep="\t", index_col=0)
        spalten[titel.get(g, p.name)] = d.iloc[:, 0]
    D = pd.DataFrame(spalten)
    proben = list(D.columns)
    Z = zmat(auf_symbole(np.log2(D.clip(lower=0) / D.sum() * 1e6 + 1)))
    f = lambda a, t: [c for c in proben if c.startswith(a) and f"_{t}_" in c]  # noqa: E731
    diff = lambda a: [c for c in proben                      # noqa: E731
                      if c.startswith(a) and not c.startswith(a + "_D0")]
    zellen = [dict(spender="GSE245585_WT1", achse="osteogen",
                   k_naiv=f("WT1", "D0"), k_diff=diff("WT1"),
                   l_naiv=f("MT1", "D0"), l_diff=diff("MT1"),
                   e2=True, laesion="RB1 +/- (Retinoblastom-Patientenlinie)")]
    return Z, zellen, proben


def _serpina3(acc: str, achse: str):
    """SERPINA3-siRNA in primary MSCs of three donors, day 0/3/7.

    The lesion is an siRNA knockdown in healthy cells -- E2 NOT fulfilled
    (addendum 1 b). The cells carry `dWT`; their `iv` is kept separately as
    an engineering response and never mixed with the lesion response.
    """
    t2g = pd.read_json(DATEN / "_protokolle_und_listen" /
                       "26_Orthogonal_S12___followup_t2g.json", typ="series")
    tar = DATEN / acc / f"{acc}_RAW.tar"
    spalten = {}
    with tarfile.open(tar) as t:
        for m in t.getmembers():
            if not m.name.endswith(".gz"):
                continue
            name = m.name.split("_", 1)[1].replace(".abundances.txt.gz", "")
            with gzip.open(t.extractfile(m), "rt", errors="replace") as f:
                d = pd.read_csv(f, sep="\t")
            d["gene"] = d.target_id.str.split(".").str[0].map(t2g)
            spalten[name] = d.dropna(subset=["gene"]).groupby("gene").tpm.sum()
    D = pd.DataFrame(spalten)
    proben = list(D.columns)
    Z = zmat(auf_symbole(np.log2(D + 1)))
    zellen = []
    for d in (1, 2, 3):
        f = lambda a, t: [c for c in proben                  # noqa: E731
                          if c == f"Day{t}_{a}_{d}"]
        dif = lambda a: [c for c in proben                   # noqa: E731
                         if c.endswith(f"_{a}_{d}") and not c.startswith("Day0")]
        zellen.append(dict(spender=f"SERPINA3_D{d}", achse=achse,
                           k_naiv=f("Control", 0), k_diff=dif("Control"),
                           l_naiv=f("KD", 0), l_diff=dif("KD"),
                           e2=False, laesion="SERPINA3-siRNA (Engineering, kein Patientendefekt)"))
    return Z, zellen, proben


def lade_GSE247491():
    return _serpina3("GSE247491", "chondrogen")


def lade_GSE247528():
    return _serpina3("GSE247528", "osteogen")


def lade_LAMA5_USC():
    """Own data: LAMA5 knockout in a healthy USC line, naive against
    chondrogenic and osteogenic.

    By rule 0.4, WT1-3 and KO9/46/75 are **clones of one line**, hence
    **one** donor. The series is an isogenic lesion series, not a donor
    series. E2 is not fulfilled (knockout in a healthy line), so the cells
    carry `dWT`; `iv` is kept separately as an engineering response.

    The matrix is the already z-standardized `Z` of the bulk run -- the
    quantity that `03_metric.R` also computes.
    """
    ordner = DATEN / "LAMA5_USC"
    Z = pd.read_csv(ordner / "LAMA5_USC_Z.csv", index_col=0)
    me = pd.read_csv(ordner / "LAMA5_USC_meta.csv")
    Z = auf_symbole(Z)
    proben = list(me["sample"])
    s = lambda g, c: list(me[(me.genotype == g) & (me.condition == c)]["sample"])  # noqa: E731
    zellen = []
    for arm, achse in (("chondro", "chondrogen"), ("osteo", "osteogen")):
        zellen.append(dict(spender="LAMA5_USC", achse=achse,
                           k_naiv=s("WT", "naiv"), k_diff=s("WT", arm),
                           l_naiv=s("KO", "naiv"), l_diff=s("KO", arm),
                           e2=False,
                           laesion="LAMA5-KO in gesunder Linie (Engineering)"))
    return Z, zellen, proben


STUDIEN = {
    "GSE218101": (lade_GSE218101, "MPS VI (ARSB), iPSC-Chondrogenese",
                  "Punkt der achtzehn"),
    "GSE221128": (lade_GSE221128, "FOP (ACVR1), iMSC-Chondrogenese",
                  "Punkt der achtzehn"),
    "GSE245585": (lade_GSE245585, "RB1, patienteneigene MSC, Osteogenese",
                  "Punkt der achtzehn"),
    "GSE247491": (lade_GSE247491, "SERPINA3-KD, MSC-Chondrogenese",
                  "Abb. S3C, studienweise"),
    "GSE247528": (lade_GSE247528, "SERPINA3-KD, MSC-Osteogenese",
                  "Abb. S3C, studienweise"),
    "LAMA5_USC": (lade_LAMA5_USC, "LAMA5-KO, eigene USC-Linie",
                  "zwei Punkte der achtzehn"),
}


# ---------------------------------------------------------------------------
def eichung(dwt: pd.Series, achse: str) -> dict:
    """The built-in positive control (§6): does the OWN differentiation
    contrast find the lineage markers of its own axis?"""
    satz_a = "OSTEOGEN" if achse == "osteogen" else "CHONDROGEN"
    a = MARKERSAETZE[satz_a]
    b = [g for n, s in MARKERSAETZE.items() if n != satz_a for g in s]
    r = kontrast(dwt.dropna(), a, b)
    r["satz_a"] = satz_a
    r["bestanden"] = bool(r.get("status") == "ok" and r["p"] < 0.05
                          and r["kontrast"] > 0)
    return r


def main() -> None:
    log("=" * 78)
    log("Phase M-D  --  cell construction, LiCl check, built-in calibration")
    log("Preregistration PRAEREG_M_D.md (2026-08-22) + addendum 1")
    log("=" * 78)

    daten, sicht, eich, licl = {}, [], [], []
    for gse, (lader, kurz, herkunft) in STUDIEN.items():
        log("\n" + "-" * 78)
        log("%s  --  %s" % (gse, kurz))
        try:
            Z, zellen, proben = lader()
        except Exception as e:                       # noqa: BLE001
            log("   LOAD ERROR: %s" % e)
            sicht.append(dict(studie=gse, zelle="-", urteil="AUS", code="A3",
                              begruendung=f"Matrix nicht lesbar: {e}"))
            continue

        l = licl_pruefung(gse, proben)
        licl.append(l)
        log("   LiCl check: %s (%d samples checked)" % (l["urteil"], l["n_proben"]))
        if l["licl_treffer"]:
            sicht.append(dict(studie=gse, zelle="-", urteil="AUS", code="L1",
                              begruendung=l["treffer"]))
            continue

        dwt_sp, iv_sp, basis_sp, meta = {}, {}, {}, []
        for z in zellen:
            fehlt = [k for k in ("k_naiv", "k_diff", "l_naiv", "l_diff")
                     if not z[k]]
            if fehlt:
                log("   %s: 2x2 incomplete (%s) -> excluded (A2)"
                    % (z["spender"], ", ".join(fehlt)))
                sicht.append(dict(studie=gse, zelle=z["spender"], urteil="AUS",
                                  code="A2", begruendung="Arm fehlt: " + ",".join(fehlt)))
                continue
            d, i, b = zelle(Z, z["k_naiv"], z["k_diff"], z["l_naiv"], z["l_diff"])
            name = f"{z['spender']}_{z['achse'][:6]}"
            nmod = len(set(MODUL.symbol) & set(d.dropna().index))
            if nmod < 60:
                log("   %s: only %d/173 module genes measurable -> excluded (E4)" % (name, nmod))
                sicht.append(dict(studie=gse, zelle=name, urteil="AUS",
                                  code="E4", begruendung=f"Modulgene {nmod}"))
                continue
            dwt_sp[name], iv_sp[name], basis_sp[name] = d, i, b
            ei = eichung(d, z["achse"])
            ei.update(studie=gse, zelle=name, spender=z["spender"],
                      achse=z["achse"], e2=z["e2"], laesion=z["laesion"],
                      n_modulegene=nmod, herkunft=herkunft,
                      n_proben=sum(len(z[k]) for k in
                                   ("k_naiv", "k_diff", "l_naiv", "l_diff")))
            eich.append(ei)
            log("   %-28s %d/173 module genes | calibration %s: contrast %+.3f, "
                "z %+5.2f, p %8.4g -> %s"
                % (name, nmod, ei["satz_a"], ei.get("kontrast", np.nan),
                   ei.get("z", np.nan), ei.get("p", np.nan),
                   "PASSED" if ei["bestanden"] else "FAILED"))
            meta.append(dict(studie=gse, zelle=name, spender=z["spender"],
                             achse=z["achse"], e2=z["e2"],
                             laesion=z["laesion"], herkunft=herkunft,
                             n_kontrollproben=len(z["k_naiv"]) + len(z["k_diff"]),
                             eichung_bestanden=ei["bestanden"]))
            sicht.append(dict(studie=gse, zelle=name,
                              urteil="EIN" if ei["bestanden"] else "AUS",
                              code="-" if ei["bestanden"] else "M6",
                              begruendung=(f"{z['laesion']}, {z['achse']}, "
                                           f"{nmod}/173 Modulgene"
                                           if ei["bestanden"]
                                           else "Eichung nicht bestanden")))
        if dwt_sp:
            daten[gse] = dict(dwt=pd.DataFrame(dwt_sp), iv=pd.DataFrame(iv_sp),
                              basis=pd.DataFrame(basis_sp).mean(axis=1),
                              meta=pd.DataFrame(meta))

    with open(AUS / "zellen.pkl", "wb") as f:
        pickle.dump(daten, f)
    pd.DataFrame(sicht).to_csv(AUS / "zellen_sichtung.csv", index=False)
    pd.DataFrame(eich).to_csv(AUS / "eichung.csv", index=False)
    pd.DataFrame(licl).to_csv(AUS / "licl_pruefung.csv", index=False)

    # module genes per cell, for verification -- like `followup_*_gene_level.csv`
    zeilen = []
    for gse, d in daten.items():
        for c in d["dwt"].columns:
            for s in MODUL.symbol:
                if s in d["dwt"].index:
                    zeilen.append(dict(studie=gse, zelle=c, symbol=s,
                                       dWT=d["dwt"].at[s, c],
                                       iv=d["iv"].at[s, c]))
    pd.DataFrame(zeilen).to_csv(AUS / "gene_level.csv.gz", index=False)

    E = pd.DataFrame(eich)
    log("\n" + "=" * 78)
    if len(E):
        ok = E[E.bestanden]
        log("Cells built: %d | calibration passed: %d | studies with "
            "passed calibration: %d" % (len(E), len(ok), ok.studie.nunique()))
        laesion = ok[ok.e2]
        log("Of these E2 fulfilled (true lesion): %d cells from %d studies"
            % (len(laesion), laesion.studie.nunique()))
        log("Minimum of the preregistration (>= 6 cells from >= 3 studies): %s"
            % ("reached" if len(laesion) >= 6 and laesion.studie.nunique() >= 3
               else "NOT reached -- abort criterion §9 applies"))
    log("=" * 78)
    (AUS / "54b_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()

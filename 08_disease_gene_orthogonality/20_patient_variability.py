# -*- coding: utf-8 -*-
"""
20_patient_variability.py -- phase M-B: does the program run in the same
direction across patients, while the lesion response is individual?

Preregistration: `preregistrations/PRAEREG_M_B.md` (2026-08-21, before the
first download). This script implements NOTHING anew: the direction test
comes from `_module.konkordanz`, the two-set contrast from
`_module.kontrast`, the marker sets from `_marker.py`, the module from
`S5_konvergente_gene.csv`.

Flow per cohort -- in this order, no other:

  1. load, map to symbols, average per patient
  2. CALIBRATION (addendum 1 of the preregistration): does the cohort find
     the textbook markers of its own tissue type?
  3. only after a passed calibration: the main computation

The quantity per gene is the **same-direction agreement between patients**

    w_g = mean over patients of sign(x_pg - median of the control group)

-- the per-patient counterpart of `delta_g` on the orthogonal layers. The
descriptive statistic of the preregistration is

    u_g = (1 + s_g * w_g) / 2 ,   U = mean over the set genes

with `s_g` the predicted direction. `U` and the concordance test are the
same quantity in two notations (U = 0.5 + 0.5 * mean signed w value); the
test is run with `_module.konkordanz`.

Output: derived_data/M_patienten/{streuung,streuung_null,eichung,
         kohorten_sichtung}.csv
"""
from __future__ import annotations

import gzip
import pathlib
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]
                       / "00_shared"))
from _marker import ADIPOGEN, CHONDROGEN, MYOGEN, NAIV, OSTEOGEN  # noqa: E402
from _module import (DATEN, ERGEBNISSE, MODUL, SEED, konkordanz,  # noqa: E402
                    kontrast, wilson)

AUS = ERGEBNISSE / "M_patienten"
AUS.mkdir(parents=True, exist_ok=True)
LOG: list[str] = []

MARKERSAETZE = {"OSTEOGEN": OSTEOGEN, "ADIPOGEN": ADIPOGEN,
                "MYOGEN": MYOGEN, "CHONDROGEN": CHONDROGEN, "NAIV": NAIV}


def log(s: str = "") -> None:
    print(s)
    LOG.append(s)


# ---------------------------------------------------------------------------
# The two fixed gene sets (preregistration §4)
# ---------------------------------------------------------------------------
def laesionssatz() -> pd.DataFrame:
    """The 173 genes of highest `iv` consistency, deterministically sorted."""
    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    K = K[K.symbol.notna() & K.iv_vz.notna()].copy()
    K["a"] = K.iv_kons.abs()
    K["b"] = K.iv_med.abs()
    K = K.sort_values(["a", "b", "ensembl"], ascending=[False, False, True])
    s = K.head(173)[["ensembl", "symbol", "iv_vz"]].rename(
        columns={"iv_vz": "ri"})
    s["ri"] = s.ri.astype(int)
    return s


# ---------------------------------------------------------------------------
# Cohort loaders.  Each returns (X, patient, gruppe):
#   X       DataFrame  symbol x sample, log scale
#   patient series      sample -> patient identifier
#   gruppe  series      sample -> "P" (patient) or "K" (control)
# ---------------------------------------------------------------------------
def _ens_symbol(idx: pd.Index) -> pd.Series:
    """Ensembl -> symbol via the project's internal gene map."""
    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    m = dict(zip(K.ensembl, K.symbol))
    return pd.Series([m.get(str(i).split(".")[0]) for i in idx], index=idx)


def _gencode_symbol(idx: pd.Index) -> pd.Series:
    """Ensembl -> symbol via the Gencode reference in data_raw/_referenz."""
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
            return pd.Series([m.get(str(i).split(".")[0]) for i in idx],
                             index=idx)
    return _ens_symbol(idx)


def _zaehlungen(df: pd.DataFrame) -> pd.DataFrame:
    """CPM, log2(x+1) -- preregistration §3."""
    s = df.sum(axis=0).replace(0, np.nan)
    return np.log2(df.div(s, axis=1) * 1e6 + 1)


def _serienmatrix(pfad: pathlib.Path) -> pd.DataFrame:
    with gzip.open(pfad, "rt", encoding="utf-8", errors="replace") as f:
        zeilen = f.read().splitlines()
    a = zeilen.index("!series_matrix_table_begin") + 1
    e = zeilen.index("!series_matrix_table_end")
    tab = [z.replace('"', "").split("\t") for z in zeilen[a:e]]
    d = pd.DataFrame(tab[1:], columns=tab[0]).set_index(tab[0][0])
    return d.apply(pd.to_numeric, errors="coerce")


def _gpl_karte(gpl: str, spalte: str) -> dict:
    """Probe -> symbol from the GEO platform annotation.

    `spalte` is the column name; on Affymetrix Gene-ST arrays the symbol
    sits in the `gene_assignment` field as the second field of a `//`
    triple.
    """
    kand = [DATEN / gpl / f"{gpl}.annot.txt", DATEN / gpl / f"{gpl}.annot.gz"]
    pfad = next(p for p in kand if p.exists())
    op = gzip.open if pfad.suffix == ".gz" else open
    kopf, karte, i = None, {}, None
    with op(pfad, "rt", encoding="utf-8", errors="replace") as f:
        for ln in f:
            if ln.startswith(("#", "^", "!")):
                continue
            t = ln.rstrip("\n").split("\t")
            if kopf is None:
                kopf = t
                if spalte not in kopf:
                    raise KeyError(f"{spalte} fehlt in {gpl}: {kopf[:12]}")
                i = kopf.index(spalte)
                continue
            if len(t) > i and t[i]:
                v = t[i]
                if spalte == "gene_assignment":
                    teile = [x.strip() for x in v.split("//")]
                    v = teile[1] if len(teile) > 1 else ""
                else:
                    v = v.split("///")[0].strip()
                if v and v != "---":
                    karte[t[0]] = v
    return karte


def _auf_symbole(d: pd.DataFrame, sym: pd.Series) -> pd.DataFrame:
    d = d.copy()
    d["symbol"] = sym.values
    d = d[d.symbol.notna() & (d.symbol != "")]
    return d.groupby("symbol").median()


def lade_GSE186141():
    x = pd.read_excel(DATEN / "GSE186141" / "GSE186141_FPKM9.6Col1.vs.2Ctrl.xlsx")
    proben = ["B14", "B15", "S501", "S144", "B23_L4", "S502", "S503", "S372"]
    X = x.set_index("gene_short_name")[proben]
    X = X.apply(pd.to_numeric, errors="coerce")
    X = np.log2(X.groupby(level=0).median() + 1)          # FPKM, log2
    grp = pd.Series({p: ("K" if p in ("B14", "B15") else "P") for p in proben})
    return X, pd.Series({p: p for p in proben}), grp


def lade_GSE292600():
    d = pd.read_csv(DATEN / "GSE292600" / "GSE292600_raw_counts.txt.gz",
                    sep="\t", index_col=0)
    X = _auf_symbole(_zaehlungen(d), _gencode_symbol(d.index))
    kontrollen = {"GD014E", "GD016E", "GD017E"}
    grp = pd.Series({c: ("K" if c in kontrollen else "P") for c in X.columns})
    return X, pd.Series({c: c for c in X.columns}), grp


def lade_GSE160207():
    d = pd.read_csv(DATEN / "GSE160207" / "GSE160207_EE_OI_RNAseq_counts.txt.gz",
                    sep="\t", index_col=0)
    sym = d.pop("symbol")
    d = d.drop(columns=[c for c in ("Chr",) if c in d])
    X = _auf_symbole(_zaehlungen(d.apply(pd.to_numeric, errors="coerce")), sym)
    grp = pd.Series({c: ("K" if c.startswith("C-") else "P") for c in X.columns})
    return X, pd.Series({c: c for c in X.columns}), grp


def lade_GSE228522():
    d = pd.read_csv(DATEN / "GSE228522" / "GSE228522_matrix_deseq2rlogv1.txt.gz",
                    sep="\t", index_col=0)
    # column keys jz2104_1..24 -> sample titles from the series matrix
    with gzip.open(DATEN / "GSE228522" / "GSE228522_series_matrix.txt.gz",
                   "rt", encoding="utf-8", errors="replace") as f:
        txt = f.read().splitlines()
    tit = [z for z in txt if z.startswith("!Sample_title")][0]
    bes = [z for z in txt if z.startswith("!Sample_description")][0]
    tit = [t.strip('"') for t in tit.split("\t")[1:]]
    bes = [t.strip('"') for t in bes.split("\t")[1:]]
    name = dict(zip(bes, tit))
    d = d.rename(columns=name)
    # only the untreated arm -- the Activin arm is an intervention, not a state
    d = d[[c for c in d.columns if "- Activin-A" in c]]
    X = _auf_symbole(np.log2(d.clip(lower=0) + 1), _gencode_symbol(d.index))
    grp = pd.Series({c: ("K" if c.startswith("Control") else "P")
                     for c in X.columns})
    return X, pd.Series({c: c for c in X.columns}), grp


def lade_GSE77758():
    d = _serienmatrix(DATEN / "GSE77758" /
                      "GSE77758-GPL6244_series_matrix.txt.gz")
    karte = _gpl_karte("GPL6244", "gene_assignment")
    X = _auf_symbole(d, pd.Series([karte.get(i) for i in d.index], index=d.index))
    # GSM2058569-73 patients, GSM2058574-78x controls (series matrix GPL6244)
    tit = _titel(DATEN / "GSE77758" / "GSE77758-GPL6244_series_matrix.txt.gz")
    X = X.rename(columns=tit)
    grp = pd.Series({c: ("K" if "Control" in c else "P") for c in X.columns})
    return X, pd.Series({c: c for c in X.columns}), grp


def lade_GSE22855():
    p = DATEN / "GSE22855" / "GSE22855_series_matrix.txt.gz"
    d = _serienmatrix(p)
    karte = _gpl_karte("GPL6884", "Symbol")
    X = _auf_symbole(d, pd.Series([karte.get(i) for i in d.index], index=d.index))
    X = X.rename(columns=_titel(p))
    kontrollen = {"L2687", "L2600", "L1142", "L1234", "L2392", "L2603"}
    grp = pd.Series({c: ("K" if c in kontrollen else "P") for c in X.columns})
    return X, pd.Series({c: c for c in X.columns}), grp


def lade_GSE58435():
    p = DATEN / "GSE58435" / "GSE58435_series_matrix.txt.gz"
    d = _serienmatrix(p)
    karte = _gpl_karte("GPL570", "Gene symbol", )
    X = _auf_symbole(d, pd.Series([karte.get(i) for i in d.index], index=d.index))
    X = X.rename(columns=_titel(p))
    grp = pd.Series({c: ("K" if c.startswith("Control") else "P")
                     for c in X.columns})
    return X, pd.Series({c: c for c in X.columns}), grp


def _titel(pfad: pathlib.Path) -> dict:
    with gzip.open(pfad, "rt", encoding="utf-8", errors="replace") as f:
        txt = f.read().splitlines()
    t = [z for z in txt if z.startswith("!Sample_title")][0]
    g = [z for z in txt if z.startswith("!Sample_geo_accession")][0]
    return dict(zip([x.strip('"') for x in g.split("\t")[1:]],
                    [x.strip('"') for x in t.split("\t")[1:]]))


KOHORTEN = {
    "GSE186141": dict(lader=lade_GSE186141, entitaet="Osteogenesis imperfecta",
                      gewebe="primaere Osteoblasten", satz_a="OSTEOGEN"),
    "GSE22855": dict(lader=lade_GSE22855, entitaet="Enchondromatose (Ollier)",
                     gewebe="Knorpel / Enchondrom", satz_a="CHONDROGEN"),
    "GSE292600": dict(lader=lade_GSE292600,
                      entitaet="akromele Dysplasie (ADAMTSL2/FBN1)",
                      gewebe="dermale Fibroblasten", satz_a="NAIV"),
    "GSE77758": dict(lader=lade_GSE77758, entitaet="EDS-HT / JHS",
                     gewebe="dermale Fibroblasten", satz_a="NAIV"),
    "GSE160207": dict(lader=lade_GSE160207, entitaet="Osteogenesis imperfecta",
                      gewebe="Vollblut", satz_a=None),
    "GSE228522": dict(lader=lade_GSE228522, entitaet="FOP (ACVR1)",
                      gewebe="CD14+-Monozyten", satz_a=None),
    "GSE58435": dict(lader=lade_GSE58435, entitaet="Turner-Syndrom",
                     gewebe="Fruchtwasser, zellfreie mRNA", satz_a=None),
}


# ---------------------------------------------------------------------------
def eichung(X: pd.DataFrame, grp: pd.Series, satz_a: str | None) -> dict:
    """Tissue identity control (preregistration, addendum 1)."""
    if satz_a is None:
        return {"status": "nicht eichbar", "bestanden": False}
    basis = X.loc[:, grp[grp == "K"].index].mean(axis=1)
    a = MARKERSAETZE[satz_a]
    b = [g for n, s in MARKERSAETZE.items() if n != satz_a for g in s]
    r = kontrast(basis, a, b)
    r["satz_a"] = satz_a
    r["bestanden"] = bool(r.get("status") == "ok" and r["p"] < 0.05
                          and r["kontrast"] > 0)
    return r


def gleichsinnigkeit(X: pd.DataFrame, grp: pd.Series) -> tuple[pd.Series, pd.Series]:
    """w_g and the baseline expression of the control group."""
    K = X.loc[:, grp[grp == "K"].index]
    P = X.loc[:, grp[grp == "P"].index]
    basis = K.median(axis=1)
    d = P.sub(basis, axis=0)
    w = np.sign(d).mean(axis=1)
    return w.dropna(), basis


def rechne(name: str, satz: pd.DataFrame, w: pd.Series,
           basis: pd.Series) -> dict:
    g = [s for s in satz.symbol if s in w.index]
    e = satz.set_index("symbol").ri.reindex(g)
    r = konkordanz(w.reindex(g), e, hintergrund=w,
                   schichtung=basis.reindex(w.index))
    r["satz"] = name
    if r.get("status") == "ok":
        u = (1 + e.values * w.reindex(g).values) / 2
        r["U"] = float(np.mean(u))
        lo, hi = wilson(int(round(r["konkordanz"] * r["n"])), r["n"])
        r["konkordanz_ki_lo"], r["konkordanz_ki_hi"] = lo, hi
        # two-set contrast: up genes against down genes of the same set
        auf = [s for s, v in zip(g, e.values) if v > 0]
        ab = [s for s, v in zip(g, e.values) if v < 0]
        k = kontrast(w, auf, ab)
        for sch, v in k.items():
            r[f"kontrast_{sch}"] = v
    return r


def main() -> None:
    log("=" * 78)
    log("Phase M-B  --  patient cohorts, program vs lesion response")
    log("Preregistration PRAEREG_M_B.md (2026-08-21), seed %d" % SEED)
    log("=" * 78)

    programm = MODUL[["ensembl", "symbol", "ri"]].copy()
    laesion = laesionssatz()
    ueber = set(programm.symbol) & set(laesion.symbol)
    log("\nProgram    %3d genes (%d up, %d down)"
        % (len(programm), (programm.ri > 0).sum(), (programm.ri < 0).sum()))
    log("Lesion     %3d genes (%d up, %d down)"
        % (len(laesion), (laesion.ri > 0).sum(), (laesion.ri < 0).sum()))
    log("Overlap of the two sets: %d genes" % len(ueber))
    laesion.to_csv(AUS / "laesionssatz_173.csv", index=False)

    eich_zeilen, erg_zeilen, sicht = [], [], []
    for gse, cfg in KOHORTEN.items():
        log("\n" + "-" * 78)
        log("%s  --  %s, %s" % (gse, cfg["entitaet"], cfg["gewebe"]))
        try:
            X, pat, grp = cfg["lader"]()
        except Exception as e:                      # noqa: BLE001
            log("   LOAD ERROR: %s" % e)
            sicht.append(dict(gse=gse, urteil="AUS", code="A3",
                              begruendung=f"Matrix nicht lesbar: {e}"))
            continue

        X = X.loc[X.var(axis=1) > 0]
        nP = int((grp == "P").sum())
        nK = int((grp == "K").sum())
        nmod = len(set(programm.symbol) & set(X.index))
        log("   %d genes, %d patients, %d controls, %d/173 module genes measurable"
            % (len(X), nP, nK, nmod))

        if nP < 5 or nK < 2 or nmod < 60:
            code = "M1" if nP < 5 else ("A2" if nK < 2 else "M5")
            log("   excluded (%s)" % code)
            sicht.append(dict(gse=gse, urteil="AUS", code=code,
                              begruendung=f"n_P={nP}, n_K={nK}, Modulgene={nmod}"))
            continue

        ei = eichung(X, grp, cfg["satz_a"])
        ei.update(gse=gse, entitaet=cfg["entitaet"], gewebe=cfg["gewebe"],
                  n_patienten=nP, n_kontrollen=nK, n_modulegene=nmod)
        eich_zeilen.append(ei)
        if ei["status"] == "nicht eichbar":
            log("   CALIBRATION: not calibratable -- no marker set of this tissue type")
        else:
            log("   CALIBRATION %s vs rest: contrast %+.4f, z %+.2f, p %.4g -> %s"
                % (ei["satz_a"], ei["kontrast"], ei["z"], ei["p"],
                   "PASSED" if ei["bestanden"] else "FAILED"))
        if not ei["bestanden"]:
            sicht.append(dict(gse=gse, urteil="AUS", code="M6",
                              begruendung=("nicht eichbar (kein Markersatz "
                                           "dieser Gewebeart)"
                                           if ei["status"] == "nicht eichbar"
                                           else "Positivkontrolle nicht bestanden")))
            log("   -> carries no finding; main computation is not reported")

        w, basis = gleichsinnigkeit(X, grp)
        for nm, satz in (("Programm", programm), ("Laesionsantwort", laesion)):
            r = rechne(nm, satz, w, basis)
            r.update(gse=gse, entitaet=cfg["entitaet"], gewebe=cfg["gewebe"],
                     n_patienten=nP, n_kontrollen=nK,
                     eichung_bestanden=ei["bestanden"])
            erg_zeilen.append(r)
            if r.get("status") != "ok":
                log("   %-16s %s" % (nm, r.get("status")))
                continue
            log("   %-16s n %3d | U %.3f | C %.3f (null %.3f) | z %+5.2f | "
                "p %8.4g | MDE80 %.3f | contrast z %+5.2f"
                % (nm, r["n"], r["U"], r["konkordanz"], r["konkordanz_null"],
                   r["konkordanz_z"], r["konkordanz_p"], r["konkordanz_mde80"],
                   r.get("kontrast_z", np.nan)))

        if ei["bestanden"]:
            sicht.append(dict(gse=gse, urteil="EIN", code="-",
                              begruendung=f"{cfg['entitaet']}, {cfg['gewebe']}, "
                                          f"n={nP} Patienten / {nK} Kontrollen"))

    E = pd.DataFrame(eich_zeilen)
    R = pd.DataFrame(erg_zeilen)
    S = pd.DataFrame(sicht)
    E.to_csv(AUS / "eichung.csv", index=False)
    R.to_csv(AUS / "streuung.csv", index=False)
    S.to_csv(AUS / "kohorten_sichtung.csv", index=False)

    # The per-cell null, stored separately (preregistration §9)
    nullspalten = ["gse", "satz", "n", "konkordanz_null", "konkordanz_null_sd",
                   "konkordanz_mde80", "rang_null", "rang_null_sd",
                   "kontrast_null_mittel", "kontrast_null_sd", "kontrast_mde80",
                   "eichung_bestanden"]
    R.reindex(columns=nullspalten).to_csv(AUS / "streuung_null.csv", index=False)

    # Study synthesis over the CALIBRATED cohorts -- descriptive, as Fig. 4B.
    # No inference test: that would require a joint cohort null model.
    G = R[R.eichung_bestanden & (R.status == "ok")]
    syn = (G.groupby("satz")
             .agg(n_kohorten=("gse", "nunique"),
                  mittlerer_z=("konkordanz_z", "mean"),
                  median_z=("konkordanz_z", "median"),
                  min_z=("konkordanz_z", "min"), max_z=("konkordanz_z", "max"),
                  mittleres_U=("U", "mean"),
                  ueber_mde80=("konkordanz", lambda s: int(
                      (s.values > G.loc[s.index, "konkordanz_mde80"].values).sum())))
             .reset_index())
    syn.to_csv(AUS / "synthese.csv", index=False)
    log("\nStudy synthesis over the calibrated cohorts (descriptive):")
    for _, r in syn.iterrows():
        log("   %-16s %d cohorts | mean z %+5.2f (%.2f .. %.2f) | "
            "mean U %.3f | above own MDE80: %d"
            % (r.satz, r.n_kohorten, r.mittlerer_z, r.min_z, r.max_z,
               r.mittleres_U, r.ueber_mde80))

    log("\n" + "=" * 78)
    n_ein = int((S.urteil == "EIN").sum()) if len(S) else 0
    log("GATE B: %d cohort(s) with passed calibration." % n_ein)
    log({0: "Layer does not carry -- no conversion from phase B (PLAN 5.4).",
         1: "Only one cohort -- finding is exploratory (rule AB1)."}.get(
        n_ein, "Two or more cohorts -- gate B passed."))
    log("=" * 78)
    (AUS / "51_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()

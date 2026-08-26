# -*- coding: utf-8 -*-
"""
52b_fetch_panels.py -- reference panels of phase M-A, with retrieval date.

Per `PRAEREG_M_A.md` §3. Loaded:

  HPO   `phenotype_to_genes.txt` (already propagated to all parent terms)
        -> HP:0002652 Skeletal dysplasia   = panel `NOSO`
        -> HP:0004322 Short stature        = panel `KLEIN`
  GWAS  GWAS Catalog, all associations; filtered on EFO_0004339
        (body height), genes from MAPPED_GENE              = panel `GWAS`
  PanelApp 309 / 1471, confidence level 3, from the existing JSON files of
        the old version                         = panels `PA309` / `PA1471`

Storage: data_raw/_panels/, output derived_data/M_humangenetik/panels.csv
"""
from __future__ import annotations

import json
import os
import pathlib
import sys
import urllib.request
from datetime import date

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from _module import DATEN, ERGEBNISSE  # noqa: E402

ROH = DATEN / "_panels"
ROH.mkdir(parents=True, exist_ok=True)
AUS = ERGEBNISSE / "M_humangenetik"
AUS.mkdir(parents=True, exist_ok=True)
HEUTE = date.today().isoformat()

HPO_URL = "https://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"
HPO_OBO = "https://purl.obolibrary.org/obo/hp.obo"
GWAS_URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
            "gwas-catalog-associations_ontology-annotated-full.zip")
# PanelApp JSON snapshots live in the repository under data_raw/_panels/panelapp;
# override with PAPER_V2_PANELAPP to read a different snapshot directory.
_env = os.environ.get("PAPER_V2_PANELAPP")
PANELAPP = pathlib.Path(_env) if _env else DATEN / "_panels" / "panelapp"


def hole(url: str, ziel: pathlib.Path) -> pathlib.Path:
    if ziel.exists() and ziel.stat().st_size > 10000:
        print(f"   vorhanden: {ziel.name}")
        return ziel
    print(f"   lade {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=900) as r, open(ziel, "wb") as f:
        f.write(r.read())
    print(f"   geladen: {ziel.name} ({ziel.stat().st_size/1e6:.1f} MB)")
    return ziel


def symbol_zu_ensembl() -> dict:
    K = pd.read_csv(ERGEBNISSE / "R_intern" / "R_interne_genkarte.csv")
    return dict(zip(K.symbol.astype(str), K.ensembl))


def nachkommen(obo: pathlib.Path, wurzel: str) -> set:
    """All descendants of an HPO term including the root, from `hp.obo`.

    `phenotype_to_genes.txt` is not propagated over the ontology; without
    this step `NOSO` would contain only the 188 genes attached literally to
    HP:0002652 instead of the genes of all dysplasia forms below it.
    """
    kind_von: dict[str, set] = {}
    akt = None
    for ln in obo.read_text("utf-8", "replace").splitlines():
        ln = ln.strip()
        if ln == "[Term]":
            akt = None
        elif ln.startswith("id: HP:"):
            akt = ln[4:]
        elif ln.startswith("is_a: HP:") and akt:
            e = ln[6:].split("!")[0].strip()
            kind_von.setdefault(e, set()).add(akt)
    aus, rand = {wurzel}, [wurzel]
    while rand:
        for kind in kind_von.get(rand.pop(), ()):
            if kind not in aus:
                aus.add(kind)
                rand.append(kind)
    return aus


def main() -> None:
    print("=" * 78)
    print("Referenzpanels M-A  --  Abrufdatum", HEUTE)
    print("=" * 78)
    s2e = symbol_zu_ensembl()
    zeilen = []

    print("\n[HPO]")
    p = hole(HPO_URL, ROH / "phenotype_to_genes.txt")
    obo = hole(HPO_OBO, ROH / "hp.obo")
    H = pd.read_csv(p, sep="\t", comment=None)
    H.columns = [c.strip("#").strip() for c in H.columns]
    spalte_hp = next(c for c in H.columns if c.lower().startswith("hpo"))
    spalte_gen = next(c for c in H.columns if "symbol" in c.lower())
    for kuerzel, hp in (("NOSO", "HP:0002652"), ("KLEIN", "HP:0004322")):
        terme = nachkommen(obo, hp)
        sym = sorted(set(H.loc[H[spalte_hp].isin(terme), spalte_gen].astype(str)))
        ens = sorted({s2e[s] for s in sym if s in s2e})
        print(f"   {kuerzel:6s} {hp} + {len(terme)-1} Subterme: "
              f"{len(sym)} Symbole, {len(ens)} im Rechnungspool")
        zeilen += [dict(panel=kuerzel,
                        quelle=f"HPO {hp} + Subterme ({len(terme)} Terme)",
                        abruf=HEUTE, symbol=s, ensembl=s2e.get(s)) for s in sym]

    print("\n[GWAS-Katalog]")
    p = hole(GWAS_URL, ROH / "gwas_catalog_full.zip")
    import zipfile
    with zipfile.ZipFile(p) as z:
        name = [n for n in z.namelist() if n.endswith((".tsv", ".txt"))][0]
        with z.open(name) as f:
            G = pd.read_csv(f, sep="\t", low_memory=False,
                            encoding="utf-8", encoding_errors="replace")
    mg = next(c for c in G.columns if c.strip() == "MAPPED_GENE")
    sel = G[G["MAPPED_TRAIT"].astype(str).str.strip().str.lower() == "body height"]
    uris = sorted(set(sel["MAPPED_TRAIT_URI"].astype(str)))
    print("   Merkmal 'body height' -> %s" % (uris[:3],))
    sym = set()
    for v in sel[mg].dropna().astype(str):
        for teil in v.replace(" - ", ",").split(","):
            t = teil.strip()
            if t and t != "NR":
                sym.add(t)
    sym = sorted(sym)
    ens = sorted({s2e[s] for s in sym if s in s2e})
    print(f"   GWAS  body height: {len(sel)} Assoziationen, "
          f"{len(sym)} Symbole, {len(ens)} im Rechnungspool")
    zeilen += [dict(panel="GWAS", quelle="GWAS-Katalog MAPPED_TRAIT=body height, MAPPED_GENE",
                    abruf=HEUTE, symbol=s, ensembl=s2e.get(s)) for s in sym]

    print("\n[PanelApp]")

    def panelapp(nr: int, nur_stufe3: bool) -> dict:
        d = json.load(open(PANELAPP / f"panel_{nr}.json", encoding="utf-8"))
        o = {}
        for g in d["genes"]:
            if nur_stufe3 and str(g.get("confidence_level")) != "3":
                continue
            e = g["gene_data"].get("ensembl_genes", {}).get("GRch38")
            if e:
                o[g["gene_data"].get("gene_symbol", "")] = \
                    list(e.values())[0]["ensembl_id"]
        return o

    for kuerzel, nr in (("PA309", 309), ("PA1471", 1471)):
        o = panelapp(nr, True)
        print(f"   {kuerzel:6s} PanelApp {nr}, Stufe 3: {len(o)} Gene "
              f"(Altfassung, unveraendert)")
        zeilen += [dict(panel=kuerzel, quelle=f"PanelApp {nr} Stufe 3",
                        abruf="2026-08-20 (bestehend)", symbol=s, ensembl=e)
                   for s, e in o.items()]

    # Nosology replacement (addendum 1 to PRAEREG_M_A): PanelApp all levels
    # combined with the HPO subtree -- the order of magnitude of the Nosology.
    for kuerzel, nr, hp in (("NOSO_BREIT", 309, "HP:0002652"),
                            ("KLEIN_BREIT", 1471, "HP:0004322")):
        o = panelapp(nr, False)
        terme = nachkommen(obo, hp)
        for s in set(H.loc[H[spalte_hp].isin(terme), spalte_gen].astype(str)):
            o.setdefault(s, s2e.get(s))
        o = {s: e for s, e in o.items() if e}
        ens = {e for e in o.values() if e in set(s2e.values())}
        print(f"   {kuerzel:11s} PanelApp {nr} alle Stufen + {hp}-Teilbaum: "
              f"{len(o)} Gene")
        zeilen += [dict(panel=kuerzel,
                        quelle=f"PanelApp {nr} alle Stufen + HPO {hp} Teilbaum",
                        abruf=HEUTE, symbol=s, ensembl=e) for s, e in o.items()]

    P = pd.DataFrame(zeilen).drop_duplicates(["panel", "symbol"])
    P.to_csv(AUS / "panels.csv", index=False)
    print("\n%d Panelzeilen -> %s" % (len(P), AUS / "panels.csv"))
    print(P.groupby("panel").ensembl.nunique().to_string())


if __name__ == "__main__":
    main()

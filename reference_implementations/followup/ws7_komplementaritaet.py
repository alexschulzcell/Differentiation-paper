"""WS7 -- are program and disease genes complementary layers? EXPLORATORY.

Program    = the 147/173 module genes (matrix output + cell-cycle exit)
Machinery  = the distal secretion machinery (S_DISTAL), on which the
             established human-genetic anchor rests (OR 3.70-5.75)

Claim tested: the program is DEPLETED for the secretion machinery while the
disease genes are ENRICHED for it -- two complementary layers. Draw
mechanism unchanged from reference_implementations/52_human_genetics_anchor.py.
"""
import sys, pathlib, importlib.util
import numpy as np, pandas as pd

W = pathlib.Path.cwd()
spec = importlib.util.spec_from_file_location("anker", W/"reference_implementations"/"52_human_genetics_anchor.py")
anker = importlib.util.module_from_spec(spec); sys.modules["anker"] = anker
spec.loader.exec_module(anker)

K = pd.read_csv(W/"derived_data"/"R_intern"/"R_interne_genkarte.csv")
HG = anker.raster(K)
go = pd.read_csv(W/"derived_data"/"M_humangenetik"/"go_saetze.csv")
panels = pd.read_csv(W/"derived_data"/"M_humangenetik"/"panels.csv")

PROGRAMM = set(K.loc[K.im_module, "ensembl"])
SAETZE = {s: set(g.ensembl) for s, g in go.groupby("satz")}
print("GO sets:", {k: len(v) for k, v in SAETZE.items()})

zeilen = []
for name, satz in SAETZE.items():
    r = anker.gematcht(PROGRAMM, satz, HG)
    r.update(seite="Programm", gegen=name); zeilen.append(r)
for pan in ["PA309", "NOSO", "NOSO_BREIT"]:
    pg = set(panels.loc[panels.panel == pan, "ensembl"].dropna())
    pg = {g for g in pg if g in HG.index}
    for name, satz in SAETZE.items():
        r = anker.gematcht(pg, satz, HG)
        r.update(seite=f"Krankheitsgene_{pan}", gegen=name); zeilen.append(r)

E = pd.DataFrame(zeilen)
E.to_csv(W/"Neu"/"derived_data"/"ws7_komplementaritaet.csv", index=False)
ok = E[E.status == "ok"]
print("\nTests:", len(ok), " Bonferroni alpha =", round(0.05/len(ok), 5))
print(ok[["seite","gegen","n_satz","n_panel","beobachtet","null_mittel",
          "z","p","OR_gematcht","OR_mde80"]].to_string(index=False))

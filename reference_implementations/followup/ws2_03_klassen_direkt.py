"""WS2 -- mechanism classes against the fixed program. EXPLORATORY.

Uses the draw mechanism from reference_implementations/52_human_genetics_anchor.py unchanged
(expression- and length-matched decile grid). No second implementation.

The question is no longer "do disease genes enrich in the program" (pooled
OR 1.00, failed), but: does the pool hide structure by disease mechanism?
"""
import sys, pathlib, importlib.util
import numpy as np, pandas as pd

W = pathlib.Path.cwd()
spec = importlib.util.spec_from_file_location("anker", W/"reference_implementations"/"52_human_genetics_anchor.py")
anker = importlib.util.module_from_spec(spec)
sys.modules["anker"] = anker
spec.loader.exec_module(anker)

K = pd.read_csv(W/"derived_data"/"R_intern"/"R_interne_genkarte.csv")
HG = anker.raster(K)
panels = pd.read_csv(W/"derived_data"/"M_humangenetik"/"panels.csv")
klassen = pd.read_csv(W/"Neu"/"derived_data"/"ws2_mechanismusklassen_go.csv")

PROGRAMM = set(K.loc[K.im_module, "ensembl"])
print("Program:", len(PROGRAMM), "background:", len(HG))

zeilen = []
for pan, gp in panels.groupby("panel"):
    pgene = set(gp.ensembl.dropna())
    for kl, gk in klassen.groupby("klasse"):
        # disease genes of THIS mechanism class
        teil = pgene & set(gk.ensembl)
        r = anker.gematcht(PROGRAMM, teil, HG)
        r.update(panel=pan, klasse=kl, n_klasse_panel=len(teil))
        zeilen.append(r)
    # reference: the whole panel (must reproduce OR ~1.00)
    r = anker.gematcht(PROGRAMM, pgene, HG); r.update(panel=pan, klasse="_GESAMT", n_klasse_panel=len(pgene))
    zeilen.append(r)

E = pd.DataFrame(zeilen)
E.to_csv(W/"Neu"/"derived_data"/"ws2_klassentest.csv", index=False)
ok = E[E.status == "ok"].copy()
ok["p_bonf"] = (ok.p * len(ok)).clip(upper=1.0)
print("\ntests computed:", len(ok), " Bonferroni threshold alpha =", round(0.05/len(ok), 5))
print(ok.sort_values("p")[["panel","klasse","n_satz","n_panel","beobachtet",
     "null_mittel","z","p","p_bonf","OR_gematcht","OR_mde80"]].head(25).to_string(index=False))

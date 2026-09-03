"""WS1 -- two layers: are disease genes constitutive infrastructure? EXPLORATORY.

The project's gene map carries only `basis_med` on the **z scale**. On it,
"constitutively highly expressed" is structurally unanswerable, and any
adjustment to it also adjusts part of the target quantity
(cor(Basis, dWT) = -0.566, `gefallene-hypothesen-guards.md`).

Everything here is therefore done on the **absolute** scale:
`expr_rank_med` is the expression rank averaged over the 18 points from
`DATEN[[i]]$expr` (rowMeans of the processed matrix). Matching runs over
deciles of this rank and of union exon length -- NOT over the z scale.

Checks
  P1  absolute expression of panel genes against a length-matched background
  P2  |dWT| of panel genes against an expression- AND length-matched background
      -> the core claim "constitutive, not dynamic"
  P3  gnomAD constraint (LOEUF), matched; plus dynamic vs constitutive
      panel genes
  P4  share of panel genes measurable at all in the 18 datasets
      (trivial alternative explanation, must be excluded)
  P5  mode of inheritance (PanelApp): monoallelic vs biallelic

Negative control throughout: a cell-cycle set. From
`Vorwissen/distale-sekretion-krankheitsgene.md` it is known that
well-studied genes generally enrich (there OR 2.1). Only **relative**
comparisons between sets are robust.
"""
import pathlib
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

W = pathlib.Path.cwd()
AUS = W / "derived_data" / "followup"
NDEZ, NZIEH, SEED = 10, 2000, 20260822
rng_global = np.random.default_rng(SEED)

G = pd.read_csv(AUS / "ws1_genkarte_erweitert.csv")
panels = pd.read_csv(W / "derived_data" / "M_humangenetik" / "panels.csv")
zellzyklus = pd.read_csv(AUS / "ws1_zellzyklus.csv")
moi = pd.read_csv(AUS / "ws1_panelapp_moi.csv")

PANELS = ["PA309", "NOSO", "NOSO_BREIT", "KLEIN", "GWAS"]
SETS = {p: set(panels.loc[panels.panel == p, "ensembl"].dropna()) for p in PANELS}
SETS["ZELLZYKLUS_NK"] = set(zellzyklus.ensembl.dropna())
PROGRAMM = set(G.loc[G.im_module, "ensembl"])
SETS["PROGRAMM"] = PROGRAMM


def raster(d: pd.DataFrame, spalten) -> pd.Series:
    """Matching cell from deciles of the given columns. Absolute scale."""
    z = pd.Series(0, index=d.index, dtype=int)
    for s in spalten:
        dez = pd.qcut(d[s].rank(method="first"), NDEZ, labels=False)
        z = z * NDEZ + dez
    return z


def gematcht_stetig(satz, HG, wert, seed=SEED, nzieh=NZIEH):
    """Median of a continuous quantity in the set against the matched null."""
    s = [g for g in satz if g in HG.index]
    if len(s) < 15:
        return {"status": "too small", "n_satz": len(s)}
    beob = float(HG.loc[s, wert].median())
    zellen = defaultdict(list)
    for i, z in enumerate(HG.zelle.values):
        zellen[z].append(i)
    zellen = {z: np.array(v) for z, v in zellen.items()}
    vals = HG[wert].values
    z_satz = HG.zelle.reindex(s).values
    rng = np.random.default_rng(seed)
    null = np.empty(nzieh)
    for i in range(nzieh):
        null[i] = np.median(vals[[rng.choice(zellen[z]) for z in z_satz]])
    mu, sd = float(null.mean()), float(null.std(ddof=1))
    zz = (beob - mu) / sd if sd > 0 else np.nan
    p_ein = (1 + int((null >= beob).sum())) / (1 + nzieh)
    p = float(min(1.0, 2 * min(p_ein, 1 - p_ein + 1 / (1 + nzieh))))
    return {"status": "ok", "n_satz": len(s), "beobachtet": beob,
            "null_mittel": mu, "null_sd": sd, "z": float(zz), "p": p,
            "mde80_wert": mu + 2.8 * sd, "mde80_delta": 2.8 * sd,
            "delta_beobachtet": beob - mu}


zeilen = []

# ---- P4 first: measurability. The trivial alternative explanation. --------
mess = []
for name, s in SETS.items():
    n_ges = len(s)
    n_mess = len(s & set(G.ensembl))
    mess.append({"satz": name, "n_panel_gesamt": n_ges, "n_in_genkarte": n_mess,
                 "anteil_messbar": n_mess / n_ges if n_ges else np.nan})
P4 = pd.DataFrame(mess)
P4.to_csv(AUS / "ws1_p4_messbarkeit.csv", index=False)
print("== P4  Measurability in the 18 datasets ==")
print(P4.to_string(index=False))

# ---- P1: absolute expression, length-matched only --------------------------
H1 = G.dropna(subset=["expr_rank_med", "laenge"]).copy()
H1["zelle"] = raster(H1, ["laenge"])
H1 = H1.set_index("ensembl")
print("\n== P1  absolute expression (rank), length-matched ==")
for name, s in SETS.items():
    r = gematcht_stetig(s, H1, "expr_rank_med")
    r.update(pruefung="P1_absolute_expression", satz=name); zeilen.append(r)
    if r["status"] == "ok":
        print(f"  {name:14s} median rank {r['beobachtet']:.3f} vs "
              f"{r['null_mittel']:.3f}  z {r['z']:+6.2f}  p {r['p']:.4f}")

# ---- P2: dynamics |dWT|, expression- AND length-matched --------------------
H2 = G.dropna(subset=["expr_rank_med", "laenge", "dWT_abs"]).copy()
H2["zelle"] = raster(H2, ["expr_rank_med", "laenge"])
H2 = H2.set_index("ensembl")
print("\n== P2  |dWT| at equal absolute expression and length ==")
for name, s in SETS.items():
    r = gematcht_stetig(s, H2, "dWT_abs")
    r.update(pruefung="P2_dynamik_dWT_abs", satz=name); zeilen.append(r)
    if r["status"] == "ok":
        print(f"  {name:14s} |dWT| {r['beobachtet']:.4f} vs "
              f"{r['null_mittel']:.4f}  z {r['z']:+6.2f}  p {r['p']:.4f}  "
              f"(MDE80 delta {r['mde80_delta']:.4f})")

# ---- P3: constraint --------------------------------------------------------
H3 = G.dropna(subset=["expr_rank_med", "laenge", "loeuf"]).copy()
H3["zelle"] = raster(H3, ["expr_rank_med", "laenge"])
H3 = H3.set_index("ensembl")
print("\n== P3  gnomAD LOEUF, expression- and length-matched "
      "(lower = more constrained) ==")
for name, s in SETS.items():
    r = gematcht_stetig(s, H3, "loeuf")
    r.update(pruefung="P3_loeuf", satz=name); zeilen.append(r)
    if r["status"] == "ok":
        print(f"  {name:14s} LOEUF {r['beobachtet']:.3f} vs "
              f"{r['null_mittel']:.3f}  z {r['z']:+6.2f}  p {r['p']:.4f}")

E = pd.DataFrame(zeilen)
E.to_csv(AUS / "ws1_p1_p3_gematcht.csv", index=False)

# ---- P3b: dynamic vs constitutive panel genes ------------------------------
print("\n== P3b  within panel genes: dynamic vs constitutive ==")
sub = []
for name in ["PA309", "NOSO_BREIT", "KLEIN"]:
    d = G[G.ensembl.isin(SETS[name])].dropna(subset=["dWT_abs", "loeuf"])
    if len(d) < 40:
        continue
    hi = d[d.dWT_abs >= d.dWT_abs.quantile(0.75)]
    lo = d[d.dWT_abs <= d.dWT_abs.quantile(0.25)]
    u, p = stats.mannwhitneyu(hi.loeuf, lo.loeuf)
    sub.append({"panel": name, "n_dynamisch": len(hi), "n_konstitutiv": len(lo),
                "loeuf_dynamisch": hi.loeuf.median(),
                "loeuf_konstitutiv": lo.loeuf.median(), "U": u, "p": p})
    print(f"  {name:12s} LOEUF dynamic {hi.loeuf.median():.3f} vs "
          f"constitutive {lo.loeuf.median():.3f}   p {p:.4f}")
pd.DataFrame(sub).to_csv(AUS / "ws1_p3b_dynamisch_gegen_konstitutiv.csv", index=False)

# ---- P5: mode of inheritance -----------------------------------------------
print("\n== P5  mode of inheritance (PanelApp), PA309 only ==")
m = moi[moi.panel == "PA309"].dropna(subset=["ensembl", "mode_of_inheritance"]).copy()
m["mono"] = m.mode_of_inheritance.str.contains("MONOALLELIC", na=False)
m["bi"] = m.mode_of_inheritance.str.contains("BIALLELIC", na=False)
m = m.merge(G[["ensembl", "dWT_abs", "loeuf", "expr_rank_med", "im_module"]],
            on="ensembl", how="inner")
rein_mono = m[m.mono & ~m.bi]
rein_bi = m[m.bi & ~m.mono]
p5 = []
for var in ["dWT_abs", "loeuf", "expr_rank_med"]:
    a = rein_mono[var].dropna(); b = rein_bi[var].dropna()
    if len(a) >= 10 and len(b) >= 10:
        u, p = stats.mannwhitneyu(a, b)
        p5.append({"variable": var, "n_mono": len(a), "n_bi": len(b),
                   "median_mono": a.median(), "median_bi": b.median(),
                   "U": u, "p": p})
        print(f"  {var:16s} mono {a.median():.3f} (n {len(a)}) vs "
              f"bi {b.median():.3f} (n {len(b)})   p {p:.4f}")
pd.DataFrame(p5).to_csv(AUS / "ws1_p5_vererbungsmodus.csv", index=False)

print("\nTests total (P1-P3):", int((E.status == "ok").sum()),
      " Bonferroni alpha =", round(0.05 / max(1, int((E.status == "ok").sum())), 5))

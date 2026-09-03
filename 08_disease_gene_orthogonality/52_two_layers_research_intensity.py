"""WS1 follow-up -- P1 and P3 with control for research intensity.
EXPLORATORY.

P1 (absolute expression) and P3 (gnomAD constraint) failed in
`51_two_layers_analysis.py` at the cell-cycle negative control: it enriched more
strongly than the disease panels in both. The suspicion -- known since
`Vorwissen/distale-sekretion-krankheitsgene.md` -- is that both quantities
in truth measure "well-studied gene".

This script tests the suspicion directly: publication count per gene from
NCBI `gene2pubmed` (retrieval 2026-08-22, tax_id 9606, mapped to Ensembl via
`gene2ensembl`) is added as a **third matching axis**.

The decision rule is fixed before the run:
  * If the cell-cycle negative control falls to the null with publication
    matching while a disease panel stays above it -> the finding of that
    panel is real and was previously only masked.
  * If both fall to the null -> P1/P3 measure research intensity, and the
    rows stay out of the paper.
  * If the negative control stays on top -> matching did not capture the
    confounder, and still nothing is claimed.
"""
import pathlib
from collections import defaultdict

import numpy as np
import pandas as pd

W = pathlib.Path.cwd()
AUS = W / "derived_data" / "followup"
NDEZ, NZIEH, SEED = 6, 2000, 20260822   # 6 deciles: 3 axes -> 216 cells

G = pd.read_csv(AUS / "ws1_genkarte_erweitert.csv")
pub = pd.read_csv(W / "data_raw" / "_referenz" / "ncbi_gene" / "gen_publikationen.csv")
panels = pd.read_csv(W / "derived_data" / "M_humangenetik" / "panels.csv")
zellzyklus = pd.read_csv(AUS / "ws1_zellzyklus.csv")

# one Ensembl gene can map to several Entrez IDs -- the publication count is
# then aggregated, otherwise duplicates arise
pub = (pub.groupby("ensembl", as_index=False).n_publikationen.max())
G = G.merge(pub[["ensembl", "n_publikationen"]], on="ensembl", how="left")
assert not G.ensembl.duplicated().any(), "duplicates after the merge"
G["log_pub"] = np.log1p(G.n_publikationen)
print("Gene map:", len(G), " of which with publication count:",
      int(G.n_publikationen.notna().sum()))

PANELS = ["PA309", "NOSO", "NOSO_BREIT", "KLEIN", "GWAS"]
SETS = {p: set(panels.loc[panels.panel == p, "ensembl"].dropna()) for p in PANELS}
SETS["ZELLZYKLUS_NK"] = set(zellzyklus.ensembl.dropna())
SETS["PROGRAMM"] = set(G.loc[G.im_module, "ensembl"])

# --- How strong is the confounder at all? ------------------------------------
print("\n== Research intensity per set (median publications) ==")
besch = []
for name, s in SETS.items():
    d = G[G.ensembl.isin(s)].n_publikationen.dropna()
    hg = G.n_publikationen.dropna()
    besch.append({"satz": name, "n": len(d), "median_pub": d.median(),
                  "median_pub_hintergrund": hg.median(),
                  "faktor": d.median() / hg.median() if hg.median() else np.nan})
    print(f"  {name:14s} n {len(d):5d}   median {d.median():7.1f}   "
          f"background {hg.median():5.1f}   factor {d.median()/hg.median():5.2f}x")
pd.DataFrame(besch).to_csv(AUS / "ws1_p6_untersuchungsintensitaet.csv", index=False)


def raster(d, spalten):
    z = pd.Series(0, index=d.index, dtype=int)
    for s in spalten:
        z = z * NDEZ + pd.qcut(d[s].rank(method="first"), NDEZ, labels=False)
    return z


def gematcht_stetig(satz, HG, wert, seed=SEED, nzieh=NZIEH):
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
            "mde80_delta": 2.8 * sd, "delta_beobachtet": beob - mu}


zeilen = []
LAEUFE = [
    # (name, target quantity, matching axes)
    ("P1_expression_ohne_pub", "expr_rank_med", ["laenge"]),
    ("P1_expression_mit_pub",  "expr_rank_med", ["laenge", "log_pub"]),
    ("P3_loeuf_ohne_pub", "loeuf", ["expr_rank_med", "laenge"]),
    ("P3_loeuf_mit_pub",  "loeuf", ["expr_rank_med", "laenge", "log_pub"]),
    ("P2_dynamik_mit_pub", "dWT_abs", ["expr_rank_med", "laenge", "log_pub"]),
]
for lauf, ziel, achsen in LAEUFE:
    H = G.dropna(subset=[ziel, "log_pub"] + achsen).copy()
    H["zelle"] = raster(H, achsen)
    H = H.set_index("ensembl")
    print(f"\n== {lauf}   matching: {' x '.join(achsen)}   (n {len(H)}) ==")
    for name, s in SETS.items():
        r = gematcht_stetig(s, H, ziel)
        r.update(lauf=lauf, zielgroesse=ziel, matching="+".join(achsen), satz=name)
        zeilen.append(r)
        if r["status"] == "ok":
            print(f"  {name:14s} {r['beobachtet']:7.4f} vs {r['null_mittel']:7.4f}"
                  f"   z {r['z']:+6.2f}   p {r['p']:.4f}   MDE80-d {r['mde80_delta']:.4f}")

E = pd.DataFrame(zeilen)
E.to_csv(AUS / "ws1_p6_publikationsmatching.csv", index=False)

print("\n== decision: what happens to the negative control? ==")
for ziel, a, b in [("expr_rank_med", "P1_expression_ohne_pub", "P1_expression_mit_pub"),
                   ("loeuf", "P3_loeuf_ohne_pub", "P3_loeuf_mit_pub")]:
    t = E[(E.lauf.isin([a, b])) & (E.status == "ok")].pivot(
        index="satz", columns="lauf", values="z")
    print(f"\n  target quantity {ziel}  (z without -> with publication matching)")
    for satz in t.index:
        print(f"    {satz:14s} {t.loc[satz, a]:+6.2f}  ->  {t.loc[satz, b]:+6.2f}")

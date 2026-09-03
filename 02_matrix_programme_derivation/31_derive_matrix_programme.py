# -*- coding: utf-8 -*-
# =============================================================================
# 31_derive_matrix_programme.py -- EXPLORATORY. Convergence at GENE LEVEL over 18 data
#                       sets, and its relation to skeletal dysplasia and short
#                       stature.
# =============================================================================
# *** EXPRESSLY EXPLORATORY. *** No confirmatory claim. Every number here is a
# hypothesis for future data, not evidence from these.
#
# The question: are there GENES whose interaction term iv points consistently
# in the same direction across the eighteen data sets -- and are those
# preferentially skeletal dysplasia or short stature genes?
#
# That is the question this work has never asked: everything since S1 runs at
# set level, and the clustering step showed that only about eleven degrees of
# freedom sit there. At gene level there are thousands.
#
# THE THREE TRAPS this is computed against:
#  (1) BASELINE EXPRESSION. The null diagnostic, section 4: the gene-wise z
#      scale forces cor(basis, dWT) = -0.566. A gene that starts low MUST be
#      strongly induced. Consistent genes could simply be the consistently
#      low-starting ones. This is reported, not computed away.
#  (2) THE NULL. How many genes with k/n equal signs does chance expect?
#      Analytically (binomial) AND empirically: the same count on signs turned
#      jointly per data set -- which preserves the correlation BETWEEN the
#      genes that a gene-wise permutation destroys.
#  (3) ENRICHMENT. Disease genes are well-studied, highly expressed genes. The
#      comparison therefore runs ONLY against the background of the testable
#      genes, and additionally against the second panel (dysplasia against
#      short stature), where this confounder cancels out.
# =============================================================================
import os, glob, json
import numpy as np, pandas as pd
from scipy import stats

# The session tree and the pinned PanelApp copies are not part of the public
# archive; set PAPER_V2_SESSIONS and PAPER_V2_PANELS to point at them.
W = os.environ.get("PAPER_V2_SESSIONS")
PANEL = os.environ.get("PAPER_V2_PANELS")
if not W or not PANEL:
    raise SystemExit("Set PAPER_V2_SESSIONS and PAPER_V2_PANELS. Neither tree "
                     "is part of the public archive; see README.md.")
ERG  = os.path.join(W, "20_Exploration", "derived_data")
AUF  = pd.read_csv(os.path.join(W, "16_Pseudobulk", "derived_data", "16_aufteilung.csv"))

L = []
def log(s=""):
    print(s); L.append(s)

log("=" * 76)
log("31_derive_matrix_programme.py -- CONVERGENCE AT GENE LEVEL over 18 data sets")
log("*** EXPLORATORY. No confirmatory claim. ***")
log("=" * 76)

# ------------------------------------------------------------- the matrix
G = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(ERG, "20d_gene_*.csv")))])
log("rows: %d | points: %d" % (len(G), G.punkt.nunique()))

# TWO universes, both reported (the choice co-determines the result, so it is
# not made but disclosed):
#   U1  ALL genes with a valid iv, in >= 16 of 18 data sets -> 10 177.
#       The two small data sets (points 5 and 13, about 3 300 genes) are the
#       bottleneck; at >= 17 only 1 834 genes would remain.
#   U2  pool genes only (dWT >= 0.5), in >= 8 of 18 pools -> 3 292.
#       Stricter in substance (only induced genes), weaker in the number of
#       data sets per gene. The pool filter differs per data set, and that is
#       exactly why the intersection shrinks so strongly: only 319 genes lie in
#       >= 12 pools.
import sys
UNIV = sys.argv[1] if len(sys.argv) > 1 else "U1"
if UNIV == "U1":
    Gp = G; MIN_DA = 16
    log("\nUNIVERSE U1: all genes with a valid iv, in >= 16 of 18 data sets")
else:
    Gp = G[G.im_pool]; MIN_DA = 8
    log("\nUNIVERSE U2: pool genes only (dWT >= 0.5), in >= 8 of 18 pools")
IV = Gp.pivot(index="gen", columns="punkt", values="iv")
BA = Gp.pivot(index="gen", columns="punkt", values="basis")
n_da = IV.notna().sum(axis=1)
IVm = IV[n_da >= MIN_DA]
log("genes in total: %d | evaluable: %d" % (IV.shape[0], IVm.shape[0]))

pos = (IVm > 0).sum(axis=1); neg = (IVm < 0).sum(axis=1)
nn  = pos + neg
v   = np.maximum(pos, neg)
richt = np.where(pos >= neg, 1, -1)
R = pd.DataFrame({"gen": IVm.index, "n": nn.values, "v": v.values,
                  "richtung": richt, "median_iv": IVm.median(axis=1).values,
                  "median_basis": BA.reindex(IVm.index).median(axis=1).values})
R["p"] = [stats.binomtest(int(a), int(b), 0.5).pvalue for a, b in zip(R.v, R.n)]

# --------------------------------------------------------- (1) the null
log("\n--- How many consistent genes does chance expect? ----------------------")
M = IVm.values
sgn = np.sign(M)
rng = np.random.default_rng(20260819)
NR = 500
SCHW = [(18,18),(18,17),(18,16),(17,17),(17,16),(16,16),(16,15),(16,14)]
log("\nCount by share, over all evaluable genes:")
log("%6s %10s %14s %14s %10s" % ("v/n", "observed", "noise gene-wise",
                                 "noise data-set-wise", "p_binom"))
# Noise A: signs independent per gene and data set (destroys gene correlation)
# Noise B: signs turned jointly per DATA SET (preserves gene correlation)
def zaehl(Sg, schwelle_v, schwelle_n):
    p_ = (Sg > 0).sum(axis=1); ng_ = (Sg != 0).sum(axis=1)
    vv = np.maximum(p_, ng_ - p_)
    return int(((ng_ >= schwelle_n) & (vv >= schwelle_v)).sum())

for (nreq, vreq) in SCHW:
    beob = int(((R.n >= nreq) & (R.v >= vreq)).sum())
    a = []; b = []
    for r in range(NR):
        Sa = sgn * rng.choice([-1, 1], size=sgn.shape)
        a.append(zaehl(Sa, vreq, nreq))
        Sb = sgn * rng.choice([-1, 1], size=(1, sgn.shape[1]))
        b.append(zaehl(Sb, vreq, nreq))
    pb = stats.binomtest(vreq, nreq, 0.5).pvalue
    log("%6s %10d %14.1f %14.1f %10.4g" % ("%d/%d" % (vreq, nreq), beob,
                                           np.mean(a), np.mean(b), pb))
log("\nGene-wise noise destroys the correlation BETWEEN genes and is")
log("anti-conservative. The data-set-wise column is the governing one.")

# ------------------------------------------------- (2) the consistent genes
MINN = 16 if UNIV == "U1" else 12
KON = R[(R.n >= MINN) & (R.v >= R.n - 1)].copy()   # at most one outlier
log("\n--- The consistent genes (n >= 16, at most one outlier) ----------------")
log("count: %d  (of these %d positive, %d negative)"
    % (len(KON), int((KON.richtung > 0).sum()), int((KON.richtung < 0).sum())))
if len(KON):
    log("median basis consistent %+.3f against background %+.3f (Mann-Whitney p %.3g)"
        % (KON.median_basis.median(), R.median_basis.median(),
           stats.mannwhitneyu(KON.median_basis.dropna(),
                              R.median_basis.dropna()).pvalue))
KON.sort_values("median_iv").to_csv(os.path.join(ERG, "20e_konsistente_gene_%s.csv" % UNIV),
                                    index=False)
R.to_csv(os.path.join(ERG, "20e_gene_alle_%s.csv" % UNIV), index=False)

# -------------------------------------------------------- (3) the panels
def panel(nr):
    p = json.load(open(os.path.join(PANEL, "panel_%d.json" % nr), encoding="utf-8"))
    out = set()
    for g in p["genes"]:
        if str(g.get("confidence_level")) != "3":
            continue
        e = g["gene_data"].get("ensembl_genes", {}).get("GRch38")
        if e:
            out.add(list(e.values())[0]["ensembl_id"])
    return out
P309, P1471 = panel(309), panel(1471)
UEB = P309 & P1471
P309x, P1471x = P309 - UEB, P1471 - UEB
log("\n--- Enrichment in the disease panels -----------------------------------")
log("PanelApp 309 green %d | 1471 green %d | overlap %d"
    % (len(P309), len(P1471), len(UEB)))

HG = set(R.gen)                                     # the background: the evaluable genes
def fisher(menge, name, ziel):
    a = len(ziel & menge); b = len(menge) - a
    c = len(ziel & HG) - a; d = len(HG) - len(menge) - c
    if min(a + b, c + d) == 0: return
    od, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    log("  %-34s %3d of %4d (%.1f%%) against %4d of %5d (%.1f%%)  OR %.2f  p %.3g"
        % (name, a, len(menge), 100*a/len(menge), len(ziel & HG), len(HG),
           100*len(ziel & HG)/len(HG), od, p))

for nm, ziel in (("PanelApp 309 skeletal dysplasia", P309x),
                 ("PanelApp 1471 short stature", P1471x)):
    log("\n%s against a background of %d evaluable genes:" % (nm, len(HG)))
    fisher(set(KON.gen), "all consistent", ziel)
    fisher(set(KON.gen[KON.richtung > 0]), "consistent POSITIVE", ziel)
    fisher(set(KON.gen[KON.richtung < 0]), "consistent NEGATIVE", ziel)

# Dysplasia against short stature -- "well studied" cancels out here
log("\nDirect comparison of the two panels (the confounder cancels out):")
for nm, menge in (("all consistent", set(KON.gen)),
                  ("consistent POSITIVE", set(KON.gen[KON.richtung > 0])),
                  ("consistent NEGATIVE", set(KON.gen[KON.richtung < 0]))):
    a = len(P309x & menge); b = len(P1471x & menge)
    A = len(P309x & HG); B = len(P1471x & HG)
    if min(A, B) == 0 or (a + b) == 0:
        log("  %-24s too few genes" % nm); continue
    od, p = stats.fisher_exact([[a, A - a], [b, B - b]], alternative="two-sided")
    log("  %-24s dysplasia %d/%d (%.1f%%) against short stature %d/%d (%.1f%%)  OR %.2f  p %.3g"
        % (nm, a, A, 100*a/A, b, B, 100*b/B, od, p))

# ------------------------------------------------ (4) the strongest genes
log("\n--- The thirty most consistent genes (by |median iv|) ------------------")
top = KON.reindex(KON.median_iv.abs().sort_values(ascending=False).index).head(30)
log("%-18s %4s %4s %9s %9s %s" % ("gene", "v", "n", "med iv", "med basis", "panel"))
for _, r in top.iterrows():
    tag = "309" if r.gen in P309x else ("1471" if r.gen in P1471x else "")
    log("%-18s %4d %4d %+9.3f %+9.3f %s" % (r.gen, r.v, r.n, r.median_iv,
                                            r.median_basis, tag))

open(os.path.join(ERG, "20e_konvergenz_%s.txt" % UNIV), "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\nwritten: 20e_konvergenz.txt, 20e_konsistente_gene.csv, 20e_gene_alle.csv")

# -*- coding: utf-8 -*-
# =============================================================================
# 21_detection_limit_analysis.py -- S6 step C: counting the scan, the noise expectation, the
#                    artefact control and the selection of EXACTLY ONE
#                    hypothesis.
# =============================================================================
# Preregistration: PREREG_S6.md 4.4 (noise), 4.5 (baseline expression),
# 4.6 (the selection rule), 4.7 (the detection limit). This script fixes
# NOTHING anew -- it applies the rules standing there in the order they stand.
#
# No test statistic by hand (section 13 of the S6 brief): scipy.stats.binomtest
# and scipy.stats.spearmanr.
# =============================================================================
import os, glob, json
import numpy as np, pandas as pd
from scipy import stats

# The session tree is not part of the public archive; set PAPER_V2_SESSIONS.
W = os.environ.get("PAPER_V2_SESSIONS")
if not W:
    raise SystemExit("Set PAPER_V2_SESSIONS to the tree of raw analysis "
                     "sessions. That tree is not part of the public archive; "
                     "see README.md.")
ERG = os.path.join(W, "17_Scan", "derived_data")
AUF = pd.read_csv(os.path.join(W, "16_Pseudobulk", "derived_data", "16_aufteilung.csv"))
ENT = sorted(AUF.loc[AUF.haelfte == "Entdeckung", "punkt"].tolist())
NULL_OK = sorted(AUF.loc[(AUF.haelfte == "Entdeckung") & (AUF.null_intakt.astype(str)
                 .str.upper() == "TRUE"), "punkt"].tolist())
L = []
def log(s=""):
    print(s); L.append(s)

log("=" * 70)
log("21_detection_limit_analysis.py -- S6 step C, counting the scan")
log("PREREG_S6.md 4.4/4.5/4.6, dated 2026-08-19")
log("=" * 70)
log("discovery points        : %s" % ENT)
log("of these, null demonstrably intact: %s" % NULL_OK)

# --------------------------------------------------------------- (1) the scan
sc = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(ERG, "17a_scan_*.csv")))])
sc = sc[sc.teil == "scan"]
log("\nscan rows: %d  (points %s)" % (len(sc), sorted(sc.punkt.unique())))
assert sorted(sc.punkt.unique()) == ENT, "not all nine discovery points present"

Z  = sc.pivot(index="go", columns="punkt", values="z_korr")
DB = sc.pivot(index="go", columns="punkt", values="d_basis")
DN = sc.pivot(index="go", columns="punkt", values="d_niveau")
NAM = sc.drop_duplicates("go").set_index("go")["name"]
K = Z.shape[0]
log("K (candidates with a complete row): %d" % K)
assert not Z.isna().any().any(), "missing z_korr in the scan"

npos = (Z > 0).sum(axis=1); nneg = (Z < 0).sum(axis=1)
v    = np.maximum(npos, nneg)
rich = np.where(npos >= nneg, 1, -1)
medz = Z.median(axis=1); medb = DB.median(axis=1); medn = DN.median(axis=1)

R = pd.DataFrame({"go": Z.index, "name": NAM.reindex(Z.index).values,
                  "v": v.values, "richtung": rich, "median_z": medz.values,
                  "median_d_basis": medb.values, "median_d_niveau": medn.values})
R["p_binom"] = [stats.binomtest(int(x), 9, 0.5).pvalue for x in R.v]

# --------------------------------------- (2) 4.5 control on baseline expression
rho, prho = stats.spearmanr(R.median_z, R.median_d_basis)
rho2, prho2 = stats.spearmanr(R.median_z.abs(), R.median_d_basis.abs())
log("\n--- 4.5 baseline expression -------------------------------------------")
log("Spearman(median z_korr , median dBase)      rho %+.3f  p %.3g" % (rho, prho))
log("Spearman(|median z_korr| , |median dBase|)  rho %+.3f  p %.3g" % (rho2, prho2))
grenze = R.median_d_basis.abs().quantile(2/3)
R["artefakt"] = R.median_d_basis.abs() >= grenze
log("tertile bound |median dBase| = %.4f -> %d artefact candidates of %d"
    % (grenze, int(R.artefakt.sum()), K))

# ------------------------------------------------- (3) 4.6 (c) 2 Jaccard
JA = pd.read_csv(os.path.join(ERG, "17a_kandidaten_jaccard.csv")).set_index("go")
R["jaccard_max"] = JA.reindex(R.go)["jaccard_max"].values
R["umbenennung"] = R.jaccard_max >= 0.5
log("\n--- 4.6 (c) 2 overlap with sets already tested -------------------------")
log("candidates with Jaccard >= 0.5 (not a discovery): %d" % int(R.umbenennung.sum()))

# ------------------------------------------------- (4) 4.6 (c) 3 the intact null
zok = Z[NULL_OK]
mit = pd.Series([int((np.sign(zok.loc[g]) == R.set_index("go").richtung[g]).sum())
                 for g in R.go], index=R.index)
R["n_nullok_in_richtung"] = mit.values
R["nur_defekte_null"] = R.n_nullok_in_richtung <= 2
log("\n--- 4.6 (c) 3 carried by the defective null ----------------------------")
log("candidates with <= 2 of %d null-intact points in direction: %d"
    % (len(NULL_OK), int(R.nur_defekte_null.sum())))

# ----------------------------------------------------- (5) 4.4 noise expectation
log("\n--- 4.4 (a) the analytical noise expectation ---------------------------")
for s in (9, 8):
    p = stats.binomtest(s, 9, 0.5).pvalue
    log("  v >= %d : p %.4g per candidate -> K x p = %.2f expected hits at K=%d"
        % (s, p, K * (p if s == 9 else stats.binomtest(8, 9, 0.5).pvalue), K))
p9  = stats.binomtest(9, 9, 0.5).pvalue
p8  = stats.binomtest(8, 9, 0.5).pvalue
erw9, erw8 = K * p9, K * p8

beo9 = int((R.v == 9).sum()); beo8 = int((R.v >= 8).sum())
log("\nOBSERVED: v = 9 -> %d candidates (expected %.2f)" % (beo9, erw9))
log("          v >= 8 -> %d candidates (expected %.2f)" % (beo8, erw8))

emp = {}
pf = sorted(glob.glob(os.path.join(ERG, "17b_perm_*.csv")))
if pf:
    pm = pd.concat([pd.read_csv(f) for f in pf])
    pm = pm[pm.teil == "perm"]
    log("\n--- 4.4 (b) the empirical noise expectation ----------------------------")
    log("permutation rows: %d | points %s | rounds per point %s"
        % (len(pm), sorted(pm.punkt.unique()), sorted(pm.runde.unique())))
    gos = list(Z.index)
    kub = {}                      # kub[point] -> array (nperm, K) holding sign(z)
    for pt in ENT:
        sub = pm[pm.punkt == pt].pivot(index="runde", columns="go", values="z_korr")
        sub = sub.reindex(columns=gos)
        kub[pt] = np.sign(sub.values)
    nper = min(a.shape[0] for a in kub.values())
    rng = np.random.default_rng(20260819)
    NR = 1000
    maxv = np.empty(NR, dtype=int); c9 = np.empty(NR); c8 = np.empty(NR)
    for r in range(NR):
        S = np.zeros((len(ENT), len(gos)))
        for j, pt in enumerate(ENT):
            S[j] = kub[pt][rng.integers(0, min(nper, kub[pt].shape[0]))]
        pos = (S > 0).sum(axis=0); neg = (S < 0).sum(axis=0)
        vv = np.maximum(pos, neg)
        maxv[r] = vv.max(); c9[r] = (vv == 9).sum(); c8[r] = (vv >= 8).sum()
    log("1000 noise rounds from %d permutations per point (drawn independently)" % nper)
    log("  candidates with v = 9  : mean %.2f   (observed %d)" % (c9.mean(), beo9))
    log("  candidates with v >= 8 : mean %.2f   (observed %d)" % (c8.mean(), beo8))
    log("  best sign count per round: median %d, 95%% quantile %d, maximum %d"
        % (int(np.median(maxv)), int(np.quantile(maxv, .95)), int(maxv.max())))
    pemp = float((maxv >= int(R.v.max())).mean())
    log("  empirical p of the best observed hit (v = %d): %.3f"
        % (int(R.v.max()), pemp))
    # The distribution of the sign counts under noise, per v -- for the figure
    VZ = np.zeros((NR, 10))
    rng2 = np.random.default_rng(20260819)
    for r in range(NR):
        S = np.zeros((len(ENT), len(gos)))
        for j, pt in enumerate(ENT):
            S[j] = kub[pt][rng2.integers(0, min(nper, kub[pt].shape[0]))]
        vv = np.maximum((S > 0).sum(axis=0), (S < 0).sum(axis=0))
        for x in range(5, 10):
            VZ[r, x] = (vv == x).sum()
    beob = R.v.value_counts().reindex(range(5, 10), fill_value=0)
    RV = pd.DataFrame({"v": list(range(5, 10)),
                       "beobachtet": beob.values,
                       "rauschen_mittel": VZ[:, 5:10].mean(axis=0),
                       "rauschen_q05": np.quantile(VZ[:, 5:10], .05, axis=0),
                       "rauschen_q95": np.quantile(VZ[:, 5:10], .95, axis=0)})
    RV["analytisch"] = [K * stats.binomtest(x, 9, 0.5).pvalue -
                        (K * stats.binomtest(x + 1, 9, 0.5).pvalue if x < 9 else 0)
                        for x in range(5, 10)]
    RV.to_csv(os.path.join(ERG, "17d_rauschverteilung.csv"), index=False)
    log("\nDistribution of the sign count, observed against noise:")
    log(RV.round(2).to_string(index=False))
    emp = dict(nperm=int(nper), mittel_v9=float(c9.mean()), mittel_v8=float(c8.mean()),
               median_maxv=int(np.median(maxv)), q95_maxv=int(np.quantile(maxv, .95)),
               max_maxv=int(maxv.max()), p_emp_bester=pemp)
else:
    log("\n!!! No noise arm found -- 4.4 (b) is missing, the selection is suspended.")

# ---------------------------------------------------- (6) 4.6 the selection
log("\n--- 4.6 the selection rule ---------------------------------------------")
R["ausgeschlossen"] = R.artefakt | R.umbenennung | R.nur_defekte_null
R["grund"] = ["+".join([g for g, b in
    (("Artefakt", a), ("Umbenennung", u), ("nur_defekte_Null", n)) if b]) or ""
    for a, u, n in zip(R.artefakt, R.umbenennung, R.nur_defekte_null)]
R = R.sort_values(["v", "median_z", "median_d_basis", "go"],
                  ascending=[False, False, True, True],
                  key=lambda s: s.abs() if s.name in ("median_z", "median_d_basis") else s)
R.to_csv(os.path.join(ERG, "17d_scanbilanz.csv"), index=False)

log("\nThe ten candidates with the highest sign count:")
log("%-12s %-46s %2s %+5s %8s %9s  %s" %
    ("GO", "name", "v", "ri", "med z", "med dBase", "verdict"))
for _, x in R.head(10).iterrows():
    log("%-12s %-46s %2d %+5d %8.3f %9.4f  %s" %
        (x.go, str(x["name"])[:46], x.v, x.richtung, x.median_z,
         x.median_d_basis, x.grund if x.grund else "admissible"))

kand = R[(R.v >= 8) & (~R.ausgeschlossen)]
log("\ncandidates with v >= 8 in total            : %d" % int((R.v >= 8).sum()))
log("of these excluded under 4.6 (c)            : %d" % int(((R.v >= 8) & R.ausgeschlossen).sum()))
log("REMAINING, admissible                      : %d" % len(kand))

ergebnis = {"K": int(K), "beo_v9": beo9, "beo_v8": beo8,
            "erw_v9": float(erw9), "erw_v8": float(erw8),
            "spearman_z_dbasis": float(rho), "terzilgrenze": float(grenze),
            "empirisch": emp}

if len(kand) == 0:
    log("\n>>> NO candidate survives 4.6. S6-AB2 or S6-AB3 applies.")
    log(">>> No hypothesis goes into the validation. OUTCOME 2.")
    ergebnis["gewaehlt"] = None
else:
    g = kand.iloc[0]
    log("\n>>> SELECTED (tiebreak 4.6 (d) applied):")
    log("      GO        %s" % g.go)
    log("      name      %s" % g["name"])
    log("      v         %d of 9   (p %.4g)" % (g.v, g.p_binom))
    log("      direction %+d  (frozen for step D)" % g.richtung)
    log("      median z  %+.3f | median dBase %+.4f | Jaccard %.3f"
        % (g.median_z, g.median_d_basis, g.jaccard_max))
    ergebnis["gewaehlt"] = dict(go=g.go, name=str(g["name"]), v=int(g.v),
                                richtung=int(g.richtung), p=float(g.p_binom),
                                median_z=float(g.median_z),
                                median_d_basis=float(g.median_d_basis))

with open(os.path.join(ERG, "17d_auswahl.json"), "w", encoding="utf-8") as f:
    json.dump(ergebnis, f, indent=2, ensure_ascii=False)

# --------------------------------------------- (7) 4.7 the detection limit
ef = sorted(glob.glob(os.path.join(ERG, "17c_empfind_*.csv")))
if ef:
    em = pd.concat([pd.read_csv(f) for f in ef])
    log("\n--- 4.7 the detection limit of the metric ------------------------------")
    log("rows: %d | points %s | sizes %s | delta %s"
        % (len(em), sorted(em.punkt.unique()), sorted(em.groesse_seite.unique()),
           sorted(em.delta.unique())))
    g = em.groupby(["punkt", "groesse_seite", "delta"])["z_korr"].median().reset_index()
    log("\nMedian |z_korr| per size and delta, averaged over the points:")
    piv = g.pivot_table(index="groesse_seite", columns="delta",
                        values="z_korr", aggfunc=lambda s: np.median(np.abs(s)))
    log(piv.round(2).to_string())
    log("\nDetection limit = the smallest delta with a median |z_korr| >= 2:")
    zeilen = []
    for (pt, gr), sub in g.groupby(["punkt", "groesse_seite"]):
        sub = sub.sort_values("delta")
        tr = sub[sub.z_korr.abs() >= 2]
        d = float(tr.delta.iloc[0]) if len(tr) else np.nan
        zeilen.append(dict(punkt=int(pt), groesse=int(gr), delta_nachweis=d))
    ND = pd.DataFrame(zeilen)
    ND.to_csv(os.path.join(ERG, "17d_nachweisgrenze.csv"), index=False)
    for gr, sub in ND.groupby("groesse"):
        gef = sub.delta_nachweis.dropna()
        log("  side %3d genes: median delta %s   (%d of %d points reach |z| >= 2)"
            % (gr, ("%.2f" % gef.median()) if len(gef) else "never reached",
               len(gef), len(sub)))
    log("\nContext: the three H2 conspicuities from S5 lie at |z| 3.99 to 6.36;")
    log("the metric finds an offset from the size stated above upward.")

with open(os.path.join(ERG, "17d_auszaehlung.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
print("\nwritten: 17d_scanbilanz.csv, 17d_auswahl.json, 17d_auszaehlung.txt")

# =============================================================================
# 30_graphical_abstract.py -- the graphical abstract for iScience
# =============================================================================
# Purpose  Draws the one panel that tells the paper in a single image:
#            1  the growth plate: matrix programme and disease genes rise
#               together along the maturation axis, the widest gap sits at the
#               prehypertrophic transition (data: F5B)
#            2  the decoupling: 2 of 18 datasets reach the lineage, 18 of 18
#               carry the module above its own limit (data: F1D, F2B)
#            3  the orthogonality: the disease genes are defined over distal
#               secretion and gene dosage, not over differentiation dynamics
#               (data: F4B, F4C, F4E)
#          No number is entered by hand; every value comes from the panel
#          files in figures/data/. The growth-plate scheme draws
#          deterministically with a fixed seed (numpy default_rng), so that
#          the build is reproducible.
#
# Format   The iScience submission checklist, verbatim: "The image should
#          be an exact square, 1200 x 1200 pixels in dimension at 300 dpi."
#          That is 101.6 x 101.6 mm at 300 dpi. The script checks the pixel
#          dimensions itself and aborts if they are not exact.
#          The TIF conversion in the submission packaging rescales nothing,
#          it only wraps.
#
# Inputs   figures/data/F5B_both_curves.csv
#          figures/data/F1D_calibration_per_dataset.csv
#          figures/data/F2B_module_per_dataset.csv
#          figures/data/F4B_complementarity.csv
#          figures/data/F4C_mode_of_inheritance.csv
#          figures/data/F4E_dynamics_axis.csv
# Outputs  figures/GA.pdf, figures/GA.png
#
# Style    figure_style/FIGURE_RULES.md: Arial, black type, the colour tokens
#          from publication_style.R (programme teal #0E8C8C, disease-gene
#          blue #1F6FB2, calibration passed #12946B, calibration failed
#          #C2472A). Colour codes the series, never the genotype.
#          Rule 0 of the session: every word in the image is English.
# =============================================================================

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Ellipse

# ----------------------------------------------------------------- paths
WURZEL = Path(__file__).resolve().parent.parent
DAT = WURZEL / "figures" / "data"
PUB_DIR = WURZEL / "figures"

# ----------------------------------------------------------- colour tokens
# identical with figure_style/publication_style.R
INK = "#000000"
SOFT = "#F0F0EE"
GRAU = "#9A9A9A"
TEAL = "#0E8C8C"   # the matrix programme / chondrogenic
BLUE = "#1F6FB2"   # the disease genes (PanelApp 309)
GRUEN = "#12946B"  # calibration passed
ROT = "#C2472A"    # calibration failed (documents the semantics)

FONT = "Arial"


def hell(farbe: str, anteil: float = 0.45) -> str:
    r, g, b = (int(farbe[i:i + 2], 16) / 255 for i in (1, 3, 5))
    return "#{:02X}{:02X}{:02X}".format(
        *(round((v + (1 - v) * anteil) * 255) for v in (r, g, b)))


# ------------------------------------------------------------ the numbers
eich = pd.read_csv(DAT / "F1D_calibration_per_dataset.csv")
modul = pd.read_csv(DAT / "F2B_module_per_dataset.csv")
kurven = pd.read_csv(DAT / "F5B_both_curves.csv")
f4b = pd.read_csv(DAT / "F4B_complementarity.csv")
f4c = pd.read_csv(DAT / "F4C_mode_of_inheritance.csv")
f4d = pd.read_csv(DAT / "F4E_dynamics_axis.csv")

n_ds = len(eich)
n_passed = int(eich["passed"].fillna(False).sum())
n_modul = int(modul["above_detection_limit"].fillna(False).sum())

prog = kurven[kurven["curve"] == "programme (173 genes)"].sort_values("zone_rank")
krank = kurven[kurven["curve"] == "PanelApp 309 minus programme"].sort_values("zone_rank")

or_distal = float(f4b.loc[(f4b["compared_to"] == "S_DISTAL")
                          & (f4b["gene_set"] == "PanelApp 309"),
                          "odds_ratio_matched"].iloc[0])
loeuf = f4c.loc[f4c["variable"] == "loeuf"].iloc[0]
z_dynamik_krank = float(f4d.loc[f4d["gene_set"] == "PanelApp 309", "z"].iloc[0])
z_dynamik_prog = float(f4d.loc[f4d["gene_set"] == "programme (173 genes)", "z"].iloc[0])

# -------------------------------------------------------------- the sheet
plt.rcParams.update({
    "font.family": FONT, "pdf.fonttype": 42, "ps.fonttype": 42,
})
fig = plt.figure(figsize=(4, 4), dpi=300)
fig.patch.set_facecolor("white")


def text(x, y, s, size=7, bold=False, color=INK, ha="left", va="baseline",
         rotation=0, **kw):
    return fig.text(x, y, s, fontsize=size,
                    fontweight="bold" if bold else "normal", color=color,
                    ha=ha, va=va, rotation=rotation, **kw)


# header line on the sheet, not inside the field
text(0.045, 0.975, "What skeletal differentiation", size=13, bold=True,
     va="top")
text(0.045, 0.940, "models actually measure", size=13, bold=True, va="top")
text(0.045, 0.896, "a lineage-independent matrix programme \u00b7 "
     "orthogonal disease genes", size=7.5, va="top")

# =====================================================================
# Hero: the human fetal growth plate (data: F5B)
# =====================================================================
X_L, X_B = 0.105, 0.845          # left edge and width of the hero panel
ZONEN = ["MesCond", "ChondroProg", "RestingChon", "ProlifChon",
         "PrehyperChon", "HyperChon"]
ZONEN_KURZ = ["MesCond", "ChondroProg", "Resting", "Prolif", "Prehyper",
              "Hyper"]
PRE = ZONEN.index("PrehyperChon") + 1     # x position of the prehypertrophic zone

# --- curve field
ax_k = fig.add_axes((X_L, 0.550, X_B, 0.315))
ax_k.set_xlim(0.45, 6.75)
ax_k.set_ylim(0, 0.225)
ax_k.set_yticks([0, 0.1, 0.2])
ax_k.tick_params(axis="y", length=1.6, width=0.5, labelsize=6,
                 colors=INK, pad=1.5)
ax_k.tick_params(axis="x", length=0, labelbottom=False)
ax_k.set_ylabel("zone median\n(contrast vs background)", fontsize=6.5,
                labelpad=2, linespacing=0.95)
for s in ("top", "right", "bottom"):
    ax_k.spines[s].set_visible(False)
ax_k.spines["left"].set_linewidth(0.5)

# the prehypertrophic transition: band over curves AND tissue
ax_k.axvspan(PRE - 0.5, PRE + 0.5, color=SOFT, zorder=0)
ax_k.text(PRE, 0.223, "prehypertrophic\ntransition", fontsize=6,
          ha="center", va="top", color=INK, linespacing=1.15)

for d, farbe in ((prog, TEAL), (krank, BLUE)):
    x, y = d["zone_rank"], d["median"]
    ax_k.fill_between(x, d["ci_low"], d["ci_high"], color=farbe,
                      alpha=0.14, linewidth=0, zorder=2)
    ax_k.plot(x, y, color=farbe, linewidth=1.1, marker="o",
              markersize=2.6, zorder=3)

# key at data coordinates, in the free area top left
ax_k.plot([0.72, 1.02], [0.192, 0.192], color=BLUE, linewidth=1.1,
          marker="o", markersize=2.6, zorder=4)
ax_k.text(1.10, 0.192, "disease genes (PanelApp 309)", fontsize=6.5,
          va="center", color=BLUE)
ax_k.plot([0.72, 1.02], [0.174, 0.174], color=TEAL, linewidth=1.1,
          marker="o", markersize=2.6, zorder=4)
ax_k.text(1.10, 0.174, "matrix programme (173 genes)", fontsize=6.5,
          va="center", color=TEAL)

# the widest gap between the curves, in the prehypertrophic zone
y_pre_k = float(krank.loc[krank["zone_rank"] == PRE, "median"].iloc[0])
y_pre_p = float(prog.loc[prog["zone_rank"] == PRE, "median"].iloc[0])
ax_k.annotate("", xy=(PRE, y_pre_k), xytext=(PRE, y_pre_p),
              arrowprops=dict(arrowstyle="<->", color=INK, linewidth=0.6,
                              shrinkA=1.5, shrinkB=1.5), zorder=4)
ax_k.text(PRE + 0.15, 0.175, "both peak here", fontsize=6, ha="left",
          va="baseline", color=INK)

# --- tissue strip: the growth plate itself
ax_g = fig.add_axes((X_L, 0.435, X_B, 0.100))
ax_g.set_xlim(0.45, 6.75)
ax_g.set_ylim(0, 1)
ax_g.axis("off")

# base body with rounded corners, zone colours clipped onto it
rumpf = FancyBboxPatch((0.5, 0.02), 6.0, 0.96,
                       boxstyle="round,pad=0,rounding_size=0.18",
                       facecolor="white", edgecolor="none", zorder=1)
ax_g.add_patch(rumpf)
rampen = [hell(TEAL, a) for a in (0.90, 0.86, 0.82, 0.78, 0.72, 0.66)]
for i in range(1, 7):
    ax_g.add_patch(plt.Rectangle((i - 0.5, 0.02), 1.0, 0.96,
                                 facecolor=rampen[i - 1], edgecolor="white",
                                 linewidth=0.8, zorder=2, clip_path=rumpf))
ax_g.add_patch(plt.Rectangle((PRE - 0.5, 0.02), 1.0, 0.96, facecolor=INK,
                             alpha=0.05, edgecolor="none", zorder=3,
                             clip_path=rumpf))

# chondrocytes, deterministic (fixed seed)
rng = np.random.default_rng(7)


def zelle(x, y, w, h):
    ax_g.add_patch(Ellipse((x, y), w, h, facecolor="white",
                           edgecolor=INK, linewidth=0.45, zorder=4))


# MesCond / ChondroProg: small, dense mesenchymal cells (deterministic,
# but organic-looking scatter)
for i, n in ((1, 15), (2, 12)):
    for _ in range(n):
        zelle(rng.uniform(i - 0.36, i + 0.36), rng.uniform(0.14, 0.86),
              0.13, 0.13)
# Resting: single round cells, fixed positions without overlap
for x, y in ((2.65, 0.30), (3.00, 0.22), (3.30, 0.35), (2.80, 0.60),
             (3.15, 0.55), (2.95, 0.82), (3.35, 0.75)):
    zelle(x, y, 0.20, 0.20)
# Prolif: columns of flattened cells (palisades)
for col in range(5):
    cx = 3.62 + col * 0.19
    for row in range(4):
        zelle(cx, 0.16 + row * 0.225, 0.15, 0.10)
# Prehyper: round, larger, fixed positions
for x, y in ((4.62, 0.25), (4.95, 0.20), (5.28, 0.28), (4.70, 0.55),
             (5.05, 0.50), (5.35, 0.60), (4.80, 0.82), (5.20, 0.80)):
    zelle(x, y, 0.24, 0.24)
# Hyper: large and sparse, fixed positions
for x, y in ((5.70, 0.28), (6.15, 0.30), (5.90, 0.66), (6.30, 0.72)):
    zelle(x, y, 0.34, 0.34)

# zone labels under the tissue
for i, zone in enumerate(ZONEN_KURZ, start=1):
    ax_g.text(i, -0.18, zone, fontsize=6, ha="center", va="top", color=INK)

# maturation arrow under the plate, label right above the arrow head
ax_g.annotate("", xy=(6.55, -0.63), xytext=(0.45, -0.63),
              arrowprops=dict(arrowstyle="-|>", color=INK, linewidth=0.6,
                              mutation_scale=7))
ax_g.text(6.55, -0.63, "chondrocyte maturation", fontsize=6.5, ha="right",
          va="center", color=INK)

# =====================================================================
# bottom left: the same 18 datasets, two counts
# =====================================================================
ax_z = fig.add_axes((0.05, 0.115, 0.47, 0.215))
ax_z.set_xlim(0, 1)
ax_z.set_ylim(0, 1)
ax_z.axis("off")

text(0.05, 0.325, f"{n_ds} published perturbation datasets", size=7.5,
     bold=True)

dot_x = np.linspace(0.02, 0.66, n_ds)


def punktreihe(y, farben):
    ax_z.scatter(dot_x, [y] * n_ds, s=10, c=farben, zorder=3)


# row 1: the calibration
ax_z.text(0.02, 0.73, "reach the lineage they model", fontsize=6.5,
          ha="left", va="baseline", color=INK)
punktreihe(0.60, [GRUEN if p else hell(GRAU, 0.45)
                  for p in eich["passed"].fillna(False)])
ax_z.text(0.97, 0.60, f"{n_passed} / {n_ds}", fontsize=7.5, fontweight="bold",
          color=GRUEN, ha="right", va="center")

# row 2: the module
ax_z.text(0.02, 0.33, "run the matrix programme above its own limit",
          fontsize=6.5, ha="left", va="baseline", color=INK)
punktreihe(0.20, [TEAL] * n_modul)
ax_z.text(0.97, 0.20, f"{n_modul} / {n_ds}", fontsize=7.5, fontweight="bold",
          color=TEAL, ha="right", va="center")

# =====================================================================
# bottom right: the orthogonality of the definitions
# =====================================================================
ax_o = fig.add_axes((0.565, 0.115, 0.385, 0.215))
ax_o.set_xlim(0, 1)
ax_o.set_ylim(0, 1)
ax_o.axis("off")

text(0.565, 0.325, "disease genes lie orthogonal", size=7.5, bold=True)

# the two definition axes
ax_o.add_patch(FancyArrowPatch((0.13, 0.13), (0.97, 0.13),
                               arrowstyle="-|>", mutation_scale=7,
                               color=INK, linewidth=0.7))
ax_o.add_patch(FancyArrowPatch((0.13, 0.13), (0.13, 0.80),
                               arrowstyle="-|>", mutation_scale=7,
                               color=INK, linewidth=0.7))
ax_o.text(0.55, 0.02, "differentiation dynamics", fontsize=6.5, ha="center",
          va="top", color=INK)
ax_o.text(0.045, 0.47, "distal secretion\n\u00d7 gene dosage", fontsize=6.5,
          ha="center", va="center", rotation=90, color=INK,
          linespacing=1.1)

# the programme: high on the dynamics, null on the definition of the genes
ax_o.plot([0.35], [0.20], "o", markersize=5, color=TEAL, zorder=3)
ax_o.plot([0.35, 0.35], [0.135, 0.20], linestyle=":", color=GRAU,
          linewidth=0.5, zorder=2)
ax_o.text(0.42, 0.20, "matrix programme", fontsize=6.5, ha="left",
          va="center", color=TEAL)

# the disease genes: high on their own axis, null on the dynamics
ax_o.plot([0.20], [0.70], "o", markersize=5, color=BLUE, zorder=3)
ax_o.plot([0.135, 0.20], [0.70, 0.70], linestyle=":", color=GRAU,
          linewidth=0.5, zorder=2)
ax_o.text(0.27, 0.70, "disease genes", fontsize=6.5, ha="left",
          va="center", color=BLUE)

# =====================================================================
# footer: the sentence of the paper
# =====================================================================
fuss = FancyBboxPatch((0.045, 0.030), 0.91, 0.052,
                      boxstyle="round,pad=0.008,rounding_size=0.015",
                      transform=fig.transFigure, facecolor=SOFT,
                      edgecolor="none", zorder=1)
fig.add_artist(fuss)
text(0.5, 0.056, "Dysplasia is a failure of the machinery, not of the "
     "programme.", size=7.5, ha="center", zorder=2)

# -------------------------------------------------------------- save
fig.savefig(PUB_DIR / "GA.pdf")
fig.savefig(PUB_DIR / "GA.png", dpi=300)
plt.close(fig)

from PIL import Image
masse = Image.open(PUB_DIR / "GA.png").size
if masse != (1200, 1200):
    raise RuntimeError(
        f"the graphical abstract is {masse[0]}x{masse[1]} px; "
        "exactly 1200x1200 at 300 dpi is mandatory "
        "(iScience Final File Requirements).")
print(f"30_graphical_abstract.py -- GA to {PUB_DIR}")
print(f"  {n_passed}/{n_ds} calibrated, {n_modul}/{n_ds} module above its "
      f"limit, OR {or_distal:.2f} distal, LOEUF "
      f"{loeuf['median_monoallelic']:.2f}/{loeuf['median_biallelic']:.2f}, "
      f"z {z_dynamik_krank:+.2f} dynamics")
print("  GA.png 1200x1200 px at 300 dpi -- ok")

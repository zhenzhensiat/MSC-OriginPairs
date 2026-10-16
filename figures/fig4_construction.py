"""fig4_construction.py —— Figure 4 (classifier construction).

Panels: a schematic of the gene-pair principle (drawing only); b leave-one-
group-out method comparison + 3x3 confusion matrix; c distribution of 100
random gene-set accuracies (random-set size and HOX limitation noted in the
legend); d the 30 gene-pair rules with weights, reported markers and SLC14A1
coverage flagged.

Data: data/loso_v2r2/loso_summary.tsv,
data/fig4b_confusion_v2.tsv, data/loso_v2r2/random_anchor_null.tsv,
data/fig4c_random_anchor_meta_v2.tsv, data/model_v2_rules.tsv +
data/fig4d_rule_weights_v2.tsv, data/fig4d_rule_genes_v2.tsv,
data/qpcr_spec_genes_v2.tsv. Python performs counting only.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (SRC_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

FIG_NAME = "Fig4"
nm = Numbers(FIG_NAME)

T_LS = "data/loso_v2r2/loso_summary.tsv"
T_CF = "data/fig4b_confusion_v2.tsv"
T_NULL = "data/loso_v2r2/random_anchor_null.tsv"
T_META = "data/fig4c_random_anchor_meta_v2.tsv"
T_RULES = "data/model_v2_rules.tsv"
T_W = "data/fig4d_rule_weights_v2.tsv"
T_G = "data/fig4d_rule_genes_v2.tsv"
T_QP = "data/qpcr_spec_genes_v2.tsv"

ls = read_table(T_LS)
cf = read_table(T_CF)
null = read_table(T_NULL)
meta = read_table(T_META)
rules = read_table(T_RULES)
wt = read_table(T_W)
genes = read_table(T_G)
qp = read_table(T_QP)

C_AD, C_BM, C_UC = SRC_COLORS["AD"], SRC_COLORS["BM"], SRC_COLORS["UC"]
CLS = ["AD", "BM", "UC"]
CLASSES_ORDER = CLS_ORD = {"AD": 0, "BM": 1, "UC": 2}

# =====================================================================
# 0. Hard checks (each aligned to the result table / known values)
# =====================================================================
assert len(ls) == 11 and abs(float(ls.loc[ls.model == "genome_main", "ba"].iloc[0]) - 0.970605759682224) < 1e-9
assert abs(float(ls.loc[ls.model == "published", "ba"].iloc[0]) - 0.979675604104601) < 1e-9
assert abs(float(ls.loc[ls.model == "HOX", "ba"].iloc[0]) - 0.964316451506124) < 1e-9
assert abs(float(ls.loc[ls.model == "glmnet", "ba"].iloc[0]) - 0.951737835153923) < 1e-9
assert abs(float(ls.loc[ls.model == "centroid", "ba"].iloc[0]) - 0.703172663152802) < 1e-9
assert int(ls["n_samples"].iloc[0]) == 174

cm = cf.pivot(index="true", columns="pred", values="n").reindex(index=CLS, columns=CLS).to_numpy()
assert cm.sum() == 174 and cm[0, 0] == 50 and cm[1, 1] == 92 and cm[2, 2] == 26
assert cm[0, 1] == 3 and cm[1, 0] == 3
assert cm[0, 2] == 0 and cm[2, 0] == 0 and cm[1, 2] == 0 and cm[2, 1] == 0

assert len(null) == 100 and null["set"].nunique() == 100
null_ba = null["ba"].to_numpy() * 100
null_mean = float(null["ba"].mean() * 100)
null_q95 = float(np.percentile(null["ba"], 95) * 100)
assert abs(null_mean - 94.7695897945153) < 1e-9 or abs(null_mean - 94.77) < 0.01
n_sets = int(meta.loc[meta.metric == "n_random_sets", "value"].iloc[0])
assert n_sets == 100

assert len(wt) == 30 and len(rules) == 30
wj = rules.merge(wt[["class", "rule_rank", "weight"]], on=["class", "rule_rank"], validate="1:1")
assert len(wj) == 30
assert abs(wj["weight"].max() - 1.0002536867679) < 1e-9

assert len(genes) == 56 and int(genes["is_hox"].astype(str).str.upper().eq("TRUE").sum()) == 0
MARKERS = set(genes.loc[genes["is_published_marker"].astype(str).str.upper().eq("TRUE"), "gene"])
assert MARKERS == {"EDIL3", "ENPP2", "EPB41L3", "FOXF1", "TGFB2", "VCAM1"}, MARKERS
assert float(qp.loc[qp["gene"] == "SLC14A1", "frac_studies_present"].iloc[0]) < 1.0
assert bool(qp.loc[qp["gene"] == "SLC14A1", "coverage_risk"].iloc[0])
rom = meta.loc[(meta.metric == "n_universe") & (meta.key == "G01"), "value"].iloc[0]
assert int(rom) == 4563

# =====================================================================
# Layout
# =====================================================================
fig = plt.figure(figsize=(180 * MM, 208 * MM))
outer = fig.add_gridspec(5, 1, height_ratios=[64, 4, 46, 12, 72], hspace=0.0,
                         left=0.132, right=0.992, top=0.975, bottom=0.025)
sec2 = outer[2].subgridspec(1, 3, width_ratios=[1.32, 0.62, 1.06], wspace=0.30)

# =====================================================================
# 4a Schematic (pure vector; illustrative values are not data)
# =====================================================================
axa = fig.add_subplot(outer[0])
axa.set_xlim(0, 180)
axa.set_ylim(0, 64)
axa.axis("off")
BAR_GREY = "#4D4D4D"


def _box(x0, x1, y0=4, y1=58):
    axa.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="white", ec="#BBBBBB",
                            lw=0.7, zorder=1))


def _arrow(x0, x1, y=30, label=""):
    axa.annotate("", xy=(x1, y), xytext=(x0, y),
                 arrowprops=dict(arrowstyle="-|>", lw=0.9, color=DGREY, shrinkA=0, shrinkB=0))
    if label:
        axa.text((x0 + x1) / 2, y + 1.6, label, ha="center", va="bottom", fontsize=6,
                 color=DGREY, linespacing=1.1)


_box(1, 52)
_box(58, 106)
_box(112, 179)
_arrow(52, 57.4, label="30 gene-pair\nrules")
_arrow(106, 111.4, label="weighted\nvote")
axa.text(26.5, 54.4, "1 · Compare gene pairs within one sample", ha="center", va="center",
         fontsize=6)
axa.text(82, 54.4, "2 · Rules vote per source", ha="center", va="center", fontsize=6)
axa.text(145.5, 54.4, "3 · Weighted score → call + confidence", ha="center", va="center",
         fontsize=6)
axa.text(90, 61.8, "MSC-OriginPairs classifier", ha="center", va="center", fontsize=7,
         color=DGREY)
# --- mini panels: two studies, one gene pair ---
for x0, hA, hB, sd in [(6, 26, 12, "Study 1"), (28, 17, 8, "Study 2")]:
    axa.text(x0 + 10, 50.5, sd, ha="center", va="center", fontsize=6, color=DGREY)
    axa.add_patch(Rectangle((x0 + 1, 19), 5, hA, fc=BAR_GREY, ec="none"))
    axa.add_patch(Rectangle((x0 + 11, 19), 5, hB, fc="#9E9E9E", ec="none"))
    axa.plot([x0, x0 + 18], [19, 19], color="#666666", lw=0.7)
    axa.text(x0 + 3.5, 17.2, "Gene A", ha="center", va="top", fontsize=6, color=DGREY)
    axa.text(x0 + 13.5, 17.2, "Gene B", ha="center", va="top", fontsize=6, color=DGREY)
    axa.text(x0 + 9, 13.2, "A > B", ha="center", va="center", fontsize=6)
axa.text(26.5, 10.9, "Levels differ between studies;\nwithin-sample order is unchanged.",
         ha="center", va="top", fontsize=6, color=DGREY, linespacing=1.3)
# --- vote tracks ---
for i, (cl, frac) in enumerate([("AD", 0.35), ("BM", 0.5), ("UC", 0.85)]):
    y = 40 - i * 9
    axa.add_patch(Rectangle((64, y - 2.2), 36, 4.4, fc="#F0F0F0", ec="none"))
    axa.add_patch(Rectangle((64, y - 2.2), 36 * frac, 4.4, fc=SRC_COLORS[cl], ec="none"))
    axa.text(62.5, y, cl, ha="right", va="center", fontsize=6)
axa.text(82, 13, "fraction of each source's\nrules satisfied by the sample",
         ha="center", va="center", fontsize=6, color=DGREY, linespacing=1.2)
# --- scores + margin ---
for x0, cl, h in [(116, "AD", 18), (127, "BM", 26), (138, "UC", 38)]:
    axa.add_patch(Rectangle((x0, 12), 7, h, fc=SRC_COLORS[cl], ec="none"))
    axa.text(x0 + 3.5, 9.6, cl, ha="center", va="top", fontsize=6, color=DGREY)
axa.annotate("", xy=(136.5, 12 + 38), xytext=(136.5, 12 + 26),
             arrowprops=dict(arrowstyle="<|-|>", lw=0.7, color="black", shrinkA=0, shrinkB=0))
axa.text(137.8, 12 + 32, "margin", fontsize=6, rotation=90, ha="left", va="center")
axa.text(152, 42, "highest score\n→ predicted source", ha="left", va="center",
         fontsize=6, color=DGREY, linespacing=1.3)
axa.text(152, 22, "margin (top − second)\n→ confidence", ha="left", va="center",
         fontsize=6, color=DGREY, linespacing=1.3)
panel_label(axa, "A", x=0.002, y=1.002)
nm.add("a", "schematic", "MSC-OriginPairs workflow (illustrative; not data)", "—",
       "pure vector schematic; color = source semantics (AD/BM/UC); no data values")

# =====================================================================
# 4b left: method comparison (BA+CI / recall)
# =====================================================================
axb = fig.add_subplot(sec2[0])
METH = [("genome_main", "Genome-wide pairs"),
        ("published", "Reported-marker pairs"),
        ("HOX", "HOX pairs"),
        ("glmnet", "Elastic net"),
        ("centroid", "Nearest centroid")]
ys = {}
for k, (model, label) in enumerate(METH):
    y = 4 - k
    r = ls.loc[ls.model == model].iloc[0]
    ys[model] = y
    axb.errorbar(float(r.ba) * 100, y + 0.24,
                 xerr=[[(float(r.ba) - float(r.ba_lo)) * 100], [(float(r.ba_hi) - float(r.ba)) * 100]],
                 fmt="o", ms=3.4, color="black", ecolor="black", elinewidth=0.8,
                 capsize=1.6, zorder=4)
    rec = [float(r.recall_AD) * 100, float(r.recall_BM) * 100, float(r.recall_UC) * 100]
    axb.plot([min(rec), max(rec)], [y - 0.24, y - 0.24], color="#CCCCCC", lw=0.8, zorder=2)
    for cl, v in zip(CLS, rec):
        axb.scatter(v, y - 0.24, s=9, c=SRC_COLORS[cl], linewidths=0, zorder=3)
    axb.text(float(r.ba_lo) * 100 - 0.8, y + 0.24, f"{float(r.ba) * 100:.1f}", fontsize=6,
             ha="right", va="center", color=DGREY)
axb.set_yticks([4 - k for k in range(5)])
axb.set_yticklabels([m[1] for m in METH], fontsize=6)
axb.tick_params(axis="y", length=0, pad=2)
axb.set_ylim(-0.72, 4.9)
axb.set_xlim(50, 103)
axb.set_xticks([50, 60, 70, 80, 90, 100])
axb.set_xlabel("Balanced accuracy / recall (%)", fontsize=6.5)
handles = [Line2D([], [], marker="o", ls="", ms=3.4, color="black",
                  label="Balanced accuracy (95% CI)"),
           Line2D([], [], marker="o", ls="", ms=3.2, color=C_AD, label="Recall: AD"),
           Line2D([], [], marker="o", ls="", ms=3.2, color=C_BM, label="Recall: BM"),
           Line2D([], [], marker="o", ls="", ms=3.2, color=C_UC, label="Recall: UC")]
axb.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, 1.02), fontsize=6,
           handletextpad=0.4, labelspacing=0.3, borderpad=0.2)
panel_label(axb, "B", x=-0.40, y=1.05)
for model, label in METH:
    r = ls.loc[ls.model == model].iloc[0]
    nm.add("b", f"ba_{model}", f"{float(r.ba) * 100:.1f}% (95%CI {float(r.ba_lo) * 100:.1f}-{float(r.ba_hi) * 100:.1f})",
           T_LS, label)
    nm.add("b", f"recall_{model}", f"AD {float(r.recall_AD) * 100:.1f} / BM {float(r.recall_BM) * 100:.1f} / UC {float(r.recall_UC) * 100:.1f}",
           T_LS, label)
    nm.add("b", f"label_{model}", label, T_LS, "panel b y-axis label (method name)")

# =====================================================================
# 4b right: 3×3 confusion matrix (LOSO readout, n=174)
# =====================================================================
axc = fig.add_subplot(sec2[1])
axc.imshow(cm, cmap="Greys", vmin=0, vmax=cm.max(), aspect="equal")
for i in range(3):
    for j in range(3):
        v = cm[i, j]
        axc.text(j, i, str(v), ha="center", va="center", fontsize=6,
                 color="white" if v > 60 else "black")
axc.set_xticks(range(3))
axc.set_xticklabels(CLS)
axc.set_yticks(range(3))
axc.set_yticklabels(CLS)
axc.set_xlabel("Predicted", fontsize=6.5)
axc.set_ylabel("Reported", fontsize=6.5)
axc.set_title("Leave-one-group-out (n = 174)", fontsize=6, color=DGREY, pad=2)
axc.tick_params(length=0)
nm.add("b", "confusion_n", 174, T_CF, "3×3 confusion (all-gene gene-pair classifier)")
for i, t in enumerate(CLS):
    for j, p in enumerate(CLS):
        nm.add("b", f"confusion_{t}2{p}", int(cm[i, j]), T_CF, "")

# =====================================================================
# 4c Random-gene null distribution
# =====================================================================
axc2 = fig.add_subplot(sec2[2])
bins = np.arange(85, 101, 1.0)
axc2.hist(null_ba, bins=bins, color="#C9C9C9", edgecolor="white", lw=0.3)
for v, lab, ls_ in [(float(ls.loc[ls.model == "published", "ba"].iloc[0]) * 100, "Reported markers", "-"),
                    (float(ls.loc[ls.model == "HOX", "ba"].iloc[0]) * 100, "HOX", "--")]:
    axc2.axvline(v, color="black", lw=1.0, ls=ls_, zorder=4)
axc2.text(85.4, axc2.get_ylim()[1] * 0.995,
         "— Reported markers 98.0%\n– – HOX 96.4%\n"
         f"Random gene sets (n = {n_sets})\nmean {null_mean:.1f}%",
         ha="left", va="top", fontsize=6, color=DGREY, linespacing=1.35)
axc2.set_xlim(85, 101)
axc2.set_xticks([85, 90, 95, 100])
axc2.yaxis.set_major_locator(plt.matplotlib.ticker.MaxNLocator(integer=True))
axc2.set_xlabel("Balanced accuracy (%)", fontsize=6.5)
axc2.set_ylabel("Random gene sets", fontsize=6.5)
panel_label(axc2, "C", x=-0.24, y=1.05)
nm.add("c", "n_random_sets", n_sets, T_NULL, "number of random gene sets (size = HOX survivors per fold)")
nm.add("c", "null_mean", f"{null_mean:.1f}%", T_NULL, "")
nm.add("c", "null_range", f"{null_ba.min():.1f}–{null_ba.max():.1f}%", T_NULL, "")
nm.add("c", "null_q95", f"{null_q95:.1f}%", T_NULL, "95th percentile of distribution")
nm.add("c", "line_published", "98.0%", T_LS, "model=published (vertical line)")
nm.add("c", "line_hox", "96.4%", T_LS, "model=HOX (dashed line)")
nm.add("c", "caption_roma", "G01-fold shared genes 4,563; HOX 9", T_META, "caption item")

# =====================================================================
# 4d Overview of 30 gene pairs
# =====================================================================
axd = fig.add_subplot(outer[4])
axd.set_xlim(0, 180)
axd.set_ylim(0, 68)
axd.axis("off")
PITCH = 3.95
TOP = 61.8
col_x = {"AD": 1.0, "BM": 61.0, "UC": 121.0}
for cl in CLS:
    sub = wj[wj["class"] == cl].sort_values("rule_rank")
    axd.add_patch(Rectangle((col_x[cl], TOP + 2.2), 2.8, 2.8, fc=SRC_COLORS[cl], ec="none"))
    axd.text(col_x[cl] + 4.2, TOP + 3.6, f"{cl} rules (n = {len(sub)})", fontsize=6.5,
             va="center")
    for i, (_, r) in enumerate(sub.iterrows()):
        y = TOP - i * PITCH
        g1, g2 = r["gene1"], r["gene2"]
        w1 = "bold" if g1 in MARKERS else "normal"
        w2 = "bold" if g2 in MARKERS else "normal"
        t1 = g1 + ("†" if g1 == "SLC14A1" else "")
        axd.text(col_x[cl] + 13.5, y, t1, fontsize=6, ha="right", va="center",
                 fontweight=w1)
        axd.text(col_x[cl] + 14.8, y, ">", fontsize=6, ha="center", va="center", color=DGREY)
        axd.text(col_x[cl] + 16.2, y, g2, fontsize=6, ha="left", va="center",
                 fontweight=w2)
        bw = 15.0 * float(r["weight"]) / float(wj["weight"].max())
        axd.add_patch(Rectangle((col_x[cl] + 33, y - 0.85), bw, 1.7,
                                fc=SRC_COLORS[cl], ec="none", alpha=0.85))
        axd.add_patch(Rectangle((col_x[cl] + 33, y - 0.85), 15.0, 1.7, fc="none",
                                ec="#DDDDDD", lw=0.4))
        axd.text(col_x[cl] + 33 + bw + 0.8, y, f"{float(r['weight']):.3f}", fontsize=6,
                 va="center", color=DGREY)
        nm.add("d", f"{cl}#{int(r['rule_rank'])}", f"{g1} > {g2}; weight={float(r['weight']):.4f}",
               T_W, ("includes reported marker " + ",".join([g for g in (g1, g2) if g in MARKERS])) if
               (g1 in MARKERS or g2 in MARKERS) else ("SLC14A1 (not detected on some microarray platforms)" if "SLC14A1" in (g1, g2) else ""))
axd.text(1.0, 2.2,
         "Bold: previously reported marker gene · † SLC14A1: not detected on some microarray platforms"
         " · bar = rule weight (normalised to 1.0)",
         fontsize=6, va="center", color=DGREY)
panel_label(axd, "D", x=-0.006, y=1.03)
nm.add("d", "summary", "56 genes; 0 HOX; 6 reported markers; 8 rules contain a marker",
       T_G, "caption item: no HOX among 56 genes")

# =====================================================================
save_fig(fig, "Fig4_construction")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "Fig4_construction.png"),
           os.path.join(FIG, "qc", "Fig4_construction_deutan.png"))
print("done")

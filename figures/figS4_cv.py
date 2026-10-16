"""figS4_cv.py —— Supplementary Figure S4: per-group cross-validation accuracy
and the accuracy-vs-k curve (k = 3-15 rules per source).

Data: data/loso_v2r2/loso_predictions.tsv, data/loso_v2r2/loso_summary.tsv.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

nm = Numbers("FigS4")
T_LP = "data/loso_v2r2/loso_predictions.tsv"
T_LS = "data/loso_v2r2/loso_summary.tsv"
lp = read_table(T_LP)
ls = read_table(T_LS)

gm = lp[lp["method"] == "genome_main"]
assert gm["study"].nunique() == 18 and len(gm) == 174 and gm["fold"].nunique() == 13
per = gm.assign(ok=(gm["pred"] == gm["true"])).groupby("study", sort=False).agg(
    n=("gsm", "size"), acc=("ok", "mean")).reset_index()
assert len(per) == 18 and abs(per["acc"].min() - 1 / 3) < 1e-9
assert per.loc[per["acc"].idxmin(), "study"] == "GSE37065"
assert abs(float(ls.loc[ls.model == "genome_main", "study_acc_min"].iloc[0]) - 1 / 3) < 1e-9
assert abs(float(ls.loc[ls.model == "genome_main", "study_acc_median"].iloc[0]) - 1.0) < 1e-9

fig = plt.figure(figsize=(180 * MM, 88 * MM))
g = fig.add_gridspec(1, 2, width_ratios=[1.12, 1.0], wspace=0.30,
                     left=0.075, right=0.99, top=0.90, bottom=0.16)

# ---- a: per-study accuracy ----
axa = fig.add_subplot(g[0])
per = per.sort_values(["acc", "n"], ascending=[True, False]).reset_index(drop=True)
xs = np.arange(len(per))
axa.scatter(xs, per["acc"] * 100, s=4 + per["n"] * 2.2, c="#333333", linewidths=0)
for i, r in per.iterrows():
    axa.text(i, r["acc"] * 100 + 3.5, f"{int(r['n'])}", fontsize=6, ha="center",
             va="bottom", color=DGREY)
axa.set_xticks(xs)
axa.set_xticklabels(per["study"], rotation=90, fontsize=6)
axa.set_ylim(20, 112)
axa.set_yticks([25, 50, 75, 100])
axa.set_xlim(-0.7, len(per) - 0.3)
axa.set_ylabel("Accuracy per held-out study (%)", fontsize=6.5)
panel_label(axa, "A", x=-0.09, y=1.03)
for _, r in per.iterrows():
    nm.add("a", f"acc_{r['study']}", f"{r['acc'] * 100:.1f}% (n={int(r['n'])})", T_LP,
           "genome_main within-fold accuracy")
nm.add("a", "note", "point size ∝ n; numbers are n; lowest = GSE37065 (1/3, n=3)", T_LP,
       "consistent with study_acc_median=1.0 / min=1/3; 13 held-out clusters map to 18 studies")

# ---- b: k flat line ----
axb = fig.add_subplot(g[1])
KS = [("genome_k3", 3), ("genome_k5", 5), ("genome_k8", 8), ("genome_k10", 10),
      ("genome_k15", 15)]
main = ls.loc[ls.model == "genome_main"].iloc[0]
for model, k in KS:
    r = ls.loc[ls.model == model].iloc[0]
    axb.errorbar(k, float(r["ba"]) * 100,
                 yerr=[[(float(r["ba"]) - float(r["ba_lo"])) * 100],
                       [(float(r["ba_hi"]) - float(r["ba"])) * 100]],
                 fmt="o", ms=3.2, color="black", ecolor="black", elinewidth=0.8,
                 capsize=1.8, zorder=4)
axb.axhline(float(main["ba"]) * 100, color=DGREY, lw=0.8, ls="--", zorder=2)
axb.text(15.4, float(main["ba"]) * 100, "auto-k\n(main)", fontsize=6, va="center",
         ha="left", color=DGREY, linespacing=1.2)
axb.set_xticks([3, 5, 8, 10, 15])
axb.set_xlim(1.2, 19.5)
axb.set_ylim(93, 100.5)
axb.set_yticks([94, 96, 98, 100])
axb.set_xlabel("Rules per class (k)", fontsize=6.5)
axb.set_ylabel("Balanced accuracy (%)", fontsize=6.5)
panel_label(axb, "B", x=-0.14, y=1.03)
for model, k in KS:
    r = ls.loc[ls.model == model].iloc[0]
    nm.add("b", f"ba_k{k}", f"{float(r['ba']) * 100:.2f}%"
           f" (95%CI {float(r['ba_lo']) * 100:.1f}-{float(r['ba_hi']) * 100:.1f})", T_LS,
           f"k={k}")
nm.add("b", "main", f"{float(main['ba']) * 100:.2f}% (auto-k)", T_LS, "dashed line")

save_fig(fig, "FigS4_cv")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "FigS4_cv.png"), os.path.join(FIG, "qc", "FigS4_cv_deutan.png"))
print("done")

"""figS5_glmnet.py —— Supplementary Figure S5: the elastic-net comparator in
the external cohort under two conventions (unclassifiable samples counted as
errors vs restricted to classifiable samples).

Data: data/confirm_secondary_summary_v2.tsv, data/confirm_results_v2.tsv.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (SRC_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

nm = Numbers("FigS5")
T_G = "data/confirm_secondary_summary_v2.tsv"
T_R = "data/confirm_results_v2.tsv"
g2 = read_table(T_G)
res = read_table(T_R)

gm_cons = g2.loc[g2["model"] == "glmnet_all_conservative"].iloc[0]
gm_eval = g2.loc[g2["model"] == "glmnet_evaluable_only"].iloc[0]
assert abs(float(gm_cons["ba"]) - 0.775292091421124) < 1e-9
assert abs(float(gm_eval["ba"]) - 1.0) < 1e-9
assert int(gm_cons["n_studies_not_evaluable"]) == 5
ba_main = float(res.loc[res["metric"] == "balanced_accuracy", "value"].iloc[0])
assert abs(ba_main - 0.957575757575758) < 1e-12

C_AD, C_BM, C_UC = SRC_COLORS["AD"], SRC_COLORS["BM"], SRC_COLORS["UC"]

rows = [("Gene-pair classifier (n = 176)", ba_main, None,
         (1.0, 0.872727272727273, 1.0)),
        ("Elastic net, conservative (n = 176)", float(gm_cons["ba"]), None,
         (float(gm_cons["recall_AD"]), float(gm_cons["recall_BM"]), float(gm_cons["recall_UC"]))),
        ("Elastic net, evaluable only (n = 146)", float(gm_eval["ba"]), None,
         (float(gm_eval["recall_AD"]), float(gm_eval["recall_BM"]), float(gm_eval["recall_UC"])))]

fig = plt.figure(figsize=(180 * MM, 74 * MM))
ax = fig.add_subplot(fig.add_gridspec(1, 1, left=0.235, right=0.985, top=0.88,
                                      bottom=0.18)[0])
for k, (lab, ba, _, rec) in enumerate(rows):
    y = 2 - k
    ax.plot([ba * 100], [y + 0.22], marker="o", ms=4, color="black", zorder=4)
    ax.text(ba * 100 - 1.2, y + 0.22, f"{ba * 100:.1f}", fontsize=6, ha="right",
            va="center", color=DGREY)
    ax.plot([min(rec) * 100, max(rec) * 100], [y - 0.24, y - 0.24], color="#CCCCCC",
            lw=0.8, zorder=2)
    for cl, v, dv in zip(["AD", "BM", "UC"], rec, [0.075, 0.0, -0.075]):
        ax.scatter([v * 100], [y - 0.24 + dv], s=9, c=SRC_COLORS[cl], linewidths=0,
                   zorder=3)                                   # slight vertical offset where values coincide
ax.text(100.8, -0.24 - 0.62, "", fontsize=6)
ax.set_yticks([2, 1, 0])
ax.set_yticklabels([r[0] for r in rows], fontsize=6)
ax.tick_params(axis="y", length=0, pad=2)
ax.set_ylim(-1.15, 2.75)
ax.set_xlim(50, 106)
ax.set_xticks([50, 60, 70, 80, 90, 100])
ax.set_xlabel("Balanced accuracy / recall (%)", fontsize=6.5)
ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, color="black",
                          label="Balanced accuracy"),
                   Line2D([], [], marker="o", ls="", ms=3.2, color=C_AD, label="Recall: AD"),
                   Line2D([], [], marker="o", ls="", ms=3.2, color=C_BM, label="Recall: BM"),
                   Line2D([], [], marker="o", ls="", ms=3.2, color=C_UC, label="Recall: UC")],
          loc="lower left", bbox_to_anchor=(0.0, 1.005), fontsize=6, ncol=4,
          handletextpad=0.4, columnspacing=1.2, borderpad=0.2)
ax.text(50, -1.02,
        "Conservative: 5 studies not evaluable (platform-missing genes) counted as errors. "
        "Recall markers offset vertically where values coincide.",
        fontsize=6, va="center", ha="left", color=DGREY)
panel_label(ax, "A", x=-0.55, y=1.05)
for lab, ba, _, rec in rows:
    nm.add("a", lab.split(" (")[0].replace(" ", "_").replace(",", ""),
           f"BA {ba * 100:.1f}% (n={lab.split('n = ')[1].rstrip(')')}); recall AD {rec[0] * 100:.1f} / "
           f"BM {rec[1] * 100:.1f} / UC {rec[2] * 100:.1f}", T_G if "Elastic" in lab else T_R,
           "glmnet descriptive (post-hoc)" if "Elastic" in lab else "gene-pair primary analysis (reference)")
nm.add("a", "glmnet_note", "5 studies not evaluable", T_G, "counted as errors under the conservative convention")

save_fig(fig, "FigS5_glmnet")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "FigS5_glmnet.png"), os.path.join(FIG, "qc", "FigS5_glmnet_deutan.png"))
print("done")

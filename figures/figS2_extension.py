"""figS2_extension.py —— Supplementary Figure S2: three-way outcomes for the
extension marker set (57 testable) and the combined set (156).

Data: data/marker_three_class_v2r.tsv, data/marker_repro_classification_ext_v2.tsv,
data/marker_repro_combined_v2.tsv.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (CL_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

nm = Numbers("FigS2")
T_OR = "data/marker_three_class_v2r.tsv"
T_EX = "data/marker_repro_classification_ext_v2.tsv"
T_CB = "data/marker_repro_combined_v2.tsv"
orig = read_table(T_OR)
ext = read_table(T_EX)
comb = read_table(T_CB)

KEYS = ["replicated", "reversed", "not_replicated"]
COLMAP = {"replicated": "reproduced_claimed", "reversed": "reproduced_opposite",
          "not_replicated": "not_reproduced"}
CL_TEXT = {"replicated": "white", "reversed": "black", "not_replicated": "black"}
CONTRASTS = ["BM-UC", "AD-UC", "BM-AD"]

cnt_o = orig["category"].value_counts().to_dict()
assert cnt_o == {"claimed": 36, "reversed": 16, "not_replicated": 47}, cnt_o
ext_t = ext[ext["class"] != "not_testable"]
assert len(ext_t) == 57
cnt_e = ext_t["class"].value_counts().to_dict()
assert cnt_e == {"not_reproduced": 46, "reproduced_opposite": 9, "reproduced_claimed": 2}
row = comb.set_index("scope")
assert list(row.columns[:4]) == ["n_testable", "reproduced_claimed", "reproduced_opposite",
                                 "not_reproduced"]
for scope, o in [("extension", cnt_e), ("combined", None)]:
    pass
assert (int(row.loc["combined", "n_testable"]) == 156 and
        int(row.loc["combined", "reproduced_claimed"]) == 38 and
        int(row.loc["combined", "reproduced_opposite"]) == 25 and
        int(row.loc["combined", "not_reproduced"]) == 93)

# totals per contrast (original + extension)
by_c = {}
for c in CONTRASTS:
    o = orig[orig["contrast"] == c]["category"].value_counts()
    e = ext_t[ext_t["contrast"] == c]["class"].value_counts()
    by_c[c] = {"replicated": int(o.get("claimed", 0)) + int(e.get("reproduced_claimed", 0)),
               "reversed": int(o.get("reversed", 0)) + int(e.get("reproduced_opposite", 0)),
               "not_replicated": int(o.get("not_replicated", 0)) + int(e.get("not_reproduced", 0))}
assert sum(v["replicated"] for v in by_c.values()) == 38
assert sum(sum(v.values()) for v in by_c.values()) == 156

fig = plt.figure(figsize=(180 * MM, 88 * MM))
row2 = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.24], wspace=0.30,
                        left=0.075, right=0.99, top=0.86, bottom=0.14)


def stacked(ax, labels, counts_list, xgap=None):
    xs = []
    for i, (lab, cnts) in enumerate(zip(labels, counts_list)):
        n = sum(cnts.values())
        bottom = 0.0
        for k in KEYS:
            v = cnts[k] / n * 100
            if v <= 0:
                continue
            ax.bar([i], [v], bottom=bottom, width=0.62, color=CL_COLORS[k],
                   edgecolor="white", lw=0.4)
            if v >= 6:
                ax.text(i, bottom + v / 2, f"{cnts[k]}", ha="center", va="center",
                        fontsize=6, color=CL_TEXT[k])
            bottom += v
        ax.text(i, 101.5, f"{n}", ha="center", va="bottom", fontsize=6, color=DGREY)
        xs.append(i)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=6)
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, 110)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-0.6, len(labels) - 0.4)


axa = fig.add_subplot(row2[0])
stacked(axa, ["Original\n(99)", "Extension\n(57)", "Combined\n(156)"],
        [{"replicated": 36, "reversed": 16, "not_replicated": 47},
         {"replicated": 2, "reversed": 9, "not_replicated": 46},
         {"replicated": 38, "reversed": 25, "not_replicated": 93}])
axa.set_ylabel("Marker tests (%)", fontsize=6.5)
panel_label(axa, "A", x=-0.20, y=1.06)
for scope, cc in [("Original (99)", {"replicated": 36, "reversed": 16, "not_replicated": 47}),
                  ("Extension (57)", {"replicated": 2, "reversed": 9, "not_replicated": 46}),
                  ("Combined (156)", {"replicated": 38, "reversed": 25, "not_replicated": 93})]:
    for k in KEYS:
        nm.add("a", f"{k}_{scope.split()[0]}", cc[k], T_OR if "Original" in scope else
               (T_EX if "Extension" in scope else T_CB), scope)

axb = fig.add_subplot(row2[1])
stacked(axb, [c.replace("-", "–") for c in CONTRASTS], [by_c[c] for c in CONTRASTS])
axb.set_ylabel("Marker tests (%)", fontsize=6.5)
axb.text(-0.6, 108.5, "n =", fontsize=6, color=DGREY, ha="right", va="bottom")
panel_label(axb, "B", x=-0.14, y=1.06)
for c in CONTRASTS:
    for k in KEYS:
        nm.add("b", f"{k}_{c}", by_c[c][k], T_CB, "combined set by contrast (original + extension)")
nm.add("a", "combined_rate", "38/156 = 24.4% (reported direction replicated)", T_CB, "figure legend / main text item")

fig.legend(handles=[Patch(fc=CL_COLORS[k], ec="none",
                          label={"replicated": "Replicated, reported direction",
                                 "reversed": "Replicated, opposite direction",
                                 "not_replicated": "Not replicated"}[k]) for k in KEYS],
           loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False,
           fontsize=6, handlelength=1.0, columnspacing=1.6)

save_fig(fig, "FigS2_extension")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "FigS2_extension.png"),
           os.path.join(FIG, "qc", "FigS2_extension_deutan.png"))
print("done")

"""figS7_confidence.py —— Supplementary Figure S7: confidence distribution of
the external-validation samples with the low-confidence threshold (10.7% below
threshold; all correct).

Data: data/confirm_predictions_v2.tsv.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

nm = Numbers("FigS7")
T_CP = "data/confirm_predictions_v2.tsv"
cp = read_table(T_CP)
m = pd.to_numeric(cp["margin"], errors="coerce")
assert int(m.notna().sum()) == 169 and int((m < 0.2494).sum()) == 18
low = cp.loc[(m < 0.2494), :]
assert (low["pred"] == low["true"]).all()

vals = m.dropna().to_numpy()
THR = 0.2494
fig = plt.figure(figsize=(180 * MM, 72 * MM))
ax = fig.add_subplot(fig.add_gridspec(1, 1, left=0.085, right=0.985, top=0.88,
                                      bottom=0.16)[0])
bins = np.arange(0, 1.001, 0.05)
n_lo, _ = np.histogram(vals[vals < THR], bins=bins)
n_hi, _ = np.histogram(vals[vals >= THR], bins=bins)
ctr = bins[:-1]
ax.bar(ctr, n_hi, width=0.045, align="edge", color="#4D4D4D", edgecolor="white", lw=0.3)
ax.bar(ctr, n_lo, width=0.045, align="edge", color="#C9C9C9", edgecolor="white", lw=0.3)
ax.axvline(THR, color="black", lw=0.9, ls="--", zorder=4)
ax.text(THR + 0.008, ax.get_ylim()[1] * 0.97, f"{THR}", fontsize=6, va="top")
ax.text(0.99, ax.get_ylim()[1] * 0.97,
        f"{int((m < THR).sum())}/{int(m.notna().sum())} = 10.7% low-confidence calls; all correct",
        fontsize=6, color=DGREY, ha="right", va="top")
ax.set_xlim(0, 1)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xlabel("Confidence (margin)", fontsize=6.5)
ax.set_ylabel("External validation samples", fontsize=6.5)
panel_label(ax, "A", x=-0.075, y=1.03)
nm.add("a", "n", 169, T_CP, "margin non-NA")
nm.add("a", "low_conf", "18/169 = 10.7%", T_CP, "margin<0.2494; all correct")
nm.add("a", "threshold", THR, "R/utils_models.R choose_tau", "dashed line")

save_fig(fig, "FigS7_confidence")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "FigS7_confidence.png"),
           os.path.join(FIG, "qc", "FigS7_confidence_deutan.png"))
print("done")

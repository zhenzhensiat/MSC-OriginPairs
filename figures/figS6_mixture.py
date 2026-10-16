"""figS6_mixture.py —— Supplementary Figure S6: admixture detection based on
changed or low-confidence calls (80% detection at 50% contamination; pure-sample
false-positive rate 3.9%).

Data: data/mixture_curve_v2.tsv.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

nm = Numbers("FigS6")
T_M = "data/mixture_curve_v2.tsv"
mc = read_table(T_M)
assert len(mc) == 6
r0 = mc.iloc[0]
assert abs(float(r0["flagged"]) - 0.0388888888888889) < 1e-9
r5 = mc.iloc[-1]
assert abs(float(r5["flagged"]) - 0.842956349206349) < 1e-9
assert abs(float(mc["detection_limit_80"].iloc[0]) - 0.5) < 1e-9

x = mc["minor_frac"].to_numpy() * 100
fig = plt.figure(figsize=(180 * MM, 76 * MM))
ax = fig.add_subplot(fig.add_gridspec(1, 1, left=0.075, right=0.985, top=0.90,
                                      bottom=0.15)[0])
ax.plot(x, mc["flagged"] * 100, color="black", lw=1.2, marker="o", ms=3.2, zorder=4)
ax.plot(x, mc["rejected"] * 100, color="#666666", lw=1.0, ls="--", marker="s", ms=2.8,
        zorder=3)
ax.axhline(80, color="#BBBBBB", lw=0.8, ls=":", zorder=1)
ax.text(15, 82.5, "80% detection", fontsize=6, va="bottom", ha="center", color=DGREY)
ax.text(50.8, float(mc["flagged"].iloc[-1]) * 100, "Flagged\n(low-confidence\nor changed call)",
        fontsize=6, va="center", ha="left", linespacing=1.25)
ax.text(50.8, float(mc["rejected"].iloc[-1]) * 100 - 6, "Rejected", fontsize=6, va="center",
        ha="left", color="#666666")
ax.annotate(f"pure-sample false positive {float(r0['flagged']) * 100:.1f}%", xy=(0.4, 4.2),
            xytext=(4, 32), fontsize=6, color=DGREY,
            arrowprops=dict(arrowstyle="-", lw=0.6, color=DGREY))
ax.set_xticks([0, 5, 10, 20, 30, 50])
ax.set_xticklabels(["0", "5", "10", "20", "30", "50"])
ax.set_xlim(-1, 64)
ax.set_ylim(0, 105)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel("Minor component (%)", fontsize=6.5)
ax.set_ylabel("Simulated samples detected (%)", fontsize=6.5)
panel_label(ax, "A", x=-0.075, y=1.03)
for _, r in mc.iterrows():
    nm.add("a", f"detect_{int(r['minor_frac'] * 100)}", f"{float(r['flagged']) * 100:.1f}%",
           T_M, "flagged (changed or low-confidence call)")
    nm.add("a", f"rejected_{int(r['minor_frac'] * 100)}", f"{float(r['rejected']) * 100:.1f}%",
           T_M, "rejected")
nm.add("a", "fp0", f"{float(r0['flagged']) * 100:.1f}%", T_M, "pure-sample false positive (x=0)")
nm.add("a", "limit80", "80% detection at 50% contamination", T_M, "detection_limit_80=0.5")

save_fig(fig, "FigS6_mixture")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "FigS6_mixture.png"), os.path.join(FIG, "qc", "FigS6_mixture_deutan.png"))
print("done")

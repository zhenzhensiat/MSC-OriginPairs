"""figS8_passage_lmm.py —— Supplementary Figure S8: passage-associated change
of individual rules (linear mixed models; descriptive only; 66 of 279 fits
singular).

Data: data/pair_stability_lmm_v2.tsv (per rule x study slope +/- 95% CI).
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

nm = Numbers("FigS8")
T_LM = "data/pair_stability_lmm_v2.tsv"
lm = read_table(T_LM)
assert len(lm) == 279 and lm["study"].nunique() == 10
drift = lm["drift"].astype(str).str.upper().eq("TRUE")
n_drift = int(drift.sum())
assert 0 < n_drift < 279

studies = sorted(lm["study"].unique())
sidx = {s: i for i, s in enumerate(studies)}
fig = plt.figure(figsize=(180 * MM, 86 * MM))
ax = fig.add_subplot(fig.add_gridspec(1, 1, left=0.07, right=0.985, top=0.89,
                                      bottom=0.17)[0])
rng = np.random.default_rng(19)
xs = np.array([sidx[s] for s in lm["study"]], dtype=float) + rng.uniform(-0.16, 0.16, len(lm))
ax.vlines(xs, lm["lo"], lm["hi"], color="#CCCCCC", lw=0.5, zorder=1)
ax.scatter(xs, lm["slope"], s=np.where(drift, 5, 3),
           c=np.where(drift, "black", "#999999"), linewidths=0, zorder=3)
ax.axhline(0, color="black", lw=0.8, ls="--", zorder=2)
ax.set_xticks(range(10))
ax.set_xticklabels(studies, rotation=90, fontsize=6)
ax.set_xlim(-0.6, 9.6)
ax.set_ylim(-1.05, 1.05)
ax.set_yticks([-1, -0.5, 0, 0.5, 1])
ax.set_xlabel("Passage-series study", fontsize=6.5)
ax.set_ylabel("Per-passage slope ± 95% CI", fontsize=6.5)
ax.text(-0.55, 0.93,
        "Descriptive only: 66 of 279 fits (23.7%) are singular (random effect variance → 0); no multiplicity correction.",
        fontsize=6, va="top", ha="left", color=DGREY)
ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=4, color="black",
                          label="Drift (95% CI excludes 0)"),
                   Line2D([], [], marker="o", ls="", ms=4, color="#999999",
                          label="No drift")],
          loc="upper right", bbox_to_anchor=(1.0, 1.03), fontsize=6,
          handletextpad=0.4, labelspacing=0.3, borderpad=0.2)
panel_label(ax, "A", x=-0.062, y=1.03)
nm.add("a", "n_rows", 279, T_LM, "rule × study (10 studies)")
nm.add("a", "n_drift", int(drift.sum()), T_LM, "CI excludes 0")
nm.add("a", "note", "66/279 (23.7%) fits singular (re-verified 2026-10-14); descriptive only", T_LM,
       "previous 'mostly singular' wording corrected on re-check")
for s, d in lm.groupby("study"):
    nm.add("a", f"study_{s}", f"n={len(d)}; drift={int(d['drift'].astype(str).str.upper().eq('TRUE').sum())}",
           T_LM, "")

save_fig(fig, "FigS8_passage_lmm")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "FigS8_passage_lmm.png"),
           os.path.join(FIG, "qc", "FigS8_passage_lmm_deutan.png"))
print("done")

"""figS3_marker_direction.py —— Supplementary Figure S3: direction-agnostic
replication, effect-size distributions and an exploratory predictor model
(flagged "exploratory" in the figure).

Data: data/sensitivity_published_two_sided_v2.tsv,
data/marker_transfer_v2.tsv, data/marker_predictors_v2.tsv.
"""
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (CL_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

nm = Numbers("FigS3")
T_TS = "data/sensitivity_published_two_sided_v2.tsv"
T_MT = "data/marker_transfer_v2.tsv"
T_PR = "data/marker_predictors_v2.tsv"
ts = read_table(T_TS)
mt = read_table(T_MT)
pr = read_table(T_PR)

# verification
r = ts.iloc[0]
assert int(r["transfer"]) == 52 and abs(float(r["rate"]) - 0.525252525252525) < 1e-9
assert abs(float(r["null_mean"]) - 0.0874141414141414) < 1e-9
tr = mt.loc[mt["transfer"].astype(str).str.upper().eq("TRUE")]
assert len(tr) == 36
eff = tr["median_logFC"].abs().to_numpy()
assert abs(np.median(eff) - 2.23384453273294) < 1e-9
summary = pr["model_summary"].iloc[0]
terms = re.findall(r"([^:;]+?): (-?\d+\.\d+) \(SE (\d+\.\d+)\)", summary)
assert len(terms) == 4, terms
coef = [(t[0].strip(), float(t[1]), float(t[2])) for t in terms]
coef = [t for t in coef if t[0] != "(Intercept)"]          # intercept row dropped
assert len(coef) == 3, coef

fig = plt.figure(figsize=(180 * MM, 84 * MM))
g = fig.add_gridspec(1, 3, width_ratios=[0.78, 1.06, 1.16], wspace=0.42,
                     left=0.065, right=0.99, top=0.90, bottom=0.16)

# ---- a: directed vs direction-agnostic ----
axa = fig.add_subplot(g[0])
C = CL_COLORS["replicated"]
axa.bar([0], [36.4], width=0.58, color=C, edgecolor="none")
axa.bar([1], [52.5], width=0.58, color=C, edgecolor="none")
axa.text(0, 37.6, "36/99", ha="center", va="bottom", fontsize=6)
axa.text(1, 53.7, "52/99", ha="center", va="bottom", fontsize=6)
axa.set_xticks([0, 1])
axa.set_xticklabels(["Reported\ndirection", "Either\ndirection"], fontsize=6)
axa.tick_params(axis="x", length=0)
axa.set_ylim(0, 64)
axa.set_yticks([0, 20, 40, 60])
axa.set_ylabel("Replicated marker tests (%)", fontsize=6.5)
axa.text(1, -7.5, "Matched null: 4.6% / 8.7%; p = 0.001", ha="center", va="top",
         fontsize=6, color=DGREY, clip_on=False)
panel_label(axa, "A", x=-0.30, y=1.05)
nm.add("a", "directed", "36/99 = 36.4%", "data/marker_transfer_v2.tsv", "claimed direction (one-sided)")
nm.add("a", "two_sided", "52/99 = 52.5%", T_TS, "direction-agnostic (two-sided)")
nm.add("a", "nulls", "4.6±2.6% / 8.7±3.4%; p=0.001 (1000 sets each)", T_TS, "")

# ---- b: effect-size distribution of replicated units ----
axb = fig.add_subplot(g[1])
bins = np.arange(0, 6.5, 0.5)
axb.hist(eff, bins=bins, color="#C9C9C9", edgecolor="white", lw=0.3)
med = float(np.median(eff))
axb.axvline(med, color="black", lw=1.1)
axb.text(med + 0.12, axb.get_ylim()[1] * 0.96, f"median {med:.2f}", fontsize=6, va="top")
axb.set_xlabel("|log2 fold change| (replicated)", fontsize=6.5)
axb.set_ylabel(f"Marker tests (n = {len(eff)})", fontsize=6.5)
axb.set_xlim(0, 6.2)
panel_label(axb, "B", x=-0.22, y=1.05)
nm.add("b", "n", len(eff), T_MT, "transfer=TRUE")
nm.add("b", "median", f"{med:.2f}", T_MT, "median |log2FC|")

# ---- c: exploratory model (coefficients ± SE) ----
axc = fig.add_subplot(g[2])
LBL = {"log2_n": "log2(min group n)", "expr_pct": "Expression percentile",
       "K": "K (studies)"}
for i, (name, val, se) in enumerate(coef):
    y = 2 - i
    axc.errorbar(val, y, xerr=se, fmt="o", ms=3.2, color="black", ecolor="black",
                 elinewidth=0.8, capsize=1.8)
    axc.text(val, y + 0.32, f"{val:.2f}", fontsize=6, ha="center", va="bottom",
             color=DGREY)
axc.axvline(0, color="#999999", lw=0.7, ls="--")
axc.set_yticks([2, 1, 0])
axc.set_yticklabels([LBL[t[0]] for t in coef], fontsize=6)
axc.tick_params(axis="y", length=0, pad=2)
axc.set_ylim(-0.7, 2.9)
axc.set_xlim(-2.6, 4.4)
axc.set_xticks([-2, 0, 2, 4])
axc.set_xlabel("Coefficient (log-odds) ± SE", fontsize=6.5)
axc.text(4.35, 2.75, "exploratory", fontsize=6, ha="right", va="top", color=DGREY)
panel_label(axc, "C", x=-0.30, y=1.05)
nm.add("c", "intercept_removed", "intercept row removed", T_PR, "")
nm.add("c", "log2_n_definition", "n_orig_min = minimum samples/donors per group in the original study (log2)", T_PR,
       "scripts/22_marker_reproducibility.R line 81; definition verified")
for name, val, se in coef:
    nm.add("c", f"coef_{name}", f"{val:+.3f} (SE {se:.3f})", T_PR, "glm exploratory; most coefficient CIs include 0")

save_fig(fig, "FigS3_direction")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
# (mostly monochrome/neutral; deutan output still produced)
deutan_sim(os.path.join(FIG, "FigS3_direction.png"),
           os.path.join(FIG, "qc", "FigS3_direction_deutan.png"))
print("done")

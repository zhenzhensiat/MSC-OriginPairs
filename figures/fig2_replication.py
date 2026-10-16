"""fig2_replication.py —— Figure 2 (marker replication; 2 rows x 5 panels,
180 x 118 mm).

Panels:
  a  per-gene variance explained by tissue of origin vs study
  b  replication rate of published markers vs expression-matched random gene
     sets (marker-specific null; discrete k/99 bars with z)
  c  three-way outcome: all tests + by source comparison
  d  three-way outcome by assay type of the original study (small categories
     grouped as "Other")
  e  per test, fraction of studies agreeing with the reported direction
     (dot strip + median line; median study count per group below the axis)

Style: three-way colours (replicated/reversed/not replicated) with contrast-safe
label backgrounds; minimum font size 6 pt. Python performs counting only.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (CL_COLORS, DGREY, LGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

FIG_NAME = "Fig2"
nm = Numbers(FIG_NAME)
rng = np.random.default_rng(7)

C_SOURCE, C_STUDY = "#1A1A1A", "#A0A0A0"      # neutral colors, not reusing source colors / three-class colors
CL_TEXT = {"replicated": "white", "reversed": "black", "not_replicated": "black"}
CL_LEGEND = {"replicated": "Replicated, reported direction",
             "reversed": "Replicated, opposite direction",
             "not_replicated": "Not replicated"}
KEYS = ["replicated", "reversed", "not_replicated"]
COLMAP = {"replicated": "claimed", "reversed": "reversed", "not_replicated": "not_replicated"}
DASH = "–"

# ---------------- data ----------------
T_VP = "data/variance_partition_v2.tsv"
T_TC = "data/marker_three_class_v2r.tsv"
T_BS = "data/marker_by_stratum_v2.tsv"
T_NR = "data/fig2c_null_rates_published_v2.tsv"
T_NS = "data/null_summary_v2.tsv"
T_F2 = "data/fig2f_direction_consistency_v2.tsv"
vp = read_table(T_VP)
tc = read_table(T_TC)
bs = read_table(T_BS)
nrate = read_table(T_NR)
nsum = read_table(T_NS)
f2 = read_table(T_F2)

cat_map = {"claimed": "replicated", "replicated": "replicated",
           "reversed": "reversed", "not_replicated": "not_replicated"}
cats = tc["category"].map(cat_map)
assert cats.notna().all(), f"unknown category: {sorted(tc['category'].unique())}"
cnt = {k: int((cats == k).sum()) for k in KEYS}
N_ALL = sum(cnt.values())
assert N_ALL == 99, cnt


def pct(x):
    return f"{x * 100:.1f}%"


# ---------------- layout ----------------
fig = plt.figure(figsize=(180 * MM, 118 * MM))
outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.62,
                         left=0.075, right=0.985, top=0.955, bottom=0.115)
row1 = outer[0].subgridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.30)
row2 = outer[1].subgridspec(1, 3, width_ratios=[1.0, 1.15, 1.0], wspace=0.42)

# ======== a: variance partition ========
ax = fig.add_subplot(row1[0])
grp_all = vp.copy()
groups = [("All genes", grp_all),
          ("HOX genes", vp[vp["group"] == "HOX"]),
          ("Reported markers", vp[vp["group"] == "published"])]
xt, xl = [], []
for gi, (lab, d) in enumerate(groups):
    n = len(d)
    for col, dx, colr, tag in [("R2_source", -0.19, C_SOURCE, "source"),
                               ("R2_study", 0.19, C_STUDY, "study")]:
        x0 = gi + dx
        v = d[col].to_numpy() * 100
        if n > 200:
            # Large samples: quantile intervals (thin line 5–95%, thick line IQR), more readable than violins for skewed data
            q05, q25, q75, q95 = np.percentile(v, [5, 25, 75, 95])
            ax.plot([x0, x0], [q05, q95], color=colr, lw=0.8, solid_capstyle="butt", zorder=2)
            ax.plot([x0, x0], [q25, q75], color=colr, lw=5.0, solid_capstyle="butt", zorder=2)
            nm.add("a", f"IQR_{col}_{lab}", f"{q25 / 100:.3f}–{q75 / 100:.3f}", T_VP, "interquartile range")
        else:
            ax.scatter(x0 + rng.uniform(-0.07, 0.07, n), v, s=4, c=colr, alpha=0.8,
                       linewidths=0, zorder=3)
        med = float(np.median(v))
        if n > 200:
            ax.scatter([x0], [med], s=9, c="white", edgecolors=colr, linewidths=0.8, zorder=6)
        else:
            ax.plot([x0 - 0.13, x0 + 0.13], [med, med], color=colr, lw=1.2, zorder=5)
        ax.text(x0, 106, f"{med:.1f}%", ha="center", va="bottom", fontsize=6)
        nm.add("a", f"median_{col}_{lab}", f"{med / 100:.3f}", T_VP,
               f"{lab} (n={n}) per-gene median; All genes = all rows in table")
    nm.add("a", f"n_genes_{lab}", n, T_VP, "")
    xt.append(gi)
    xl.append(f"{lab}\n(n = {n:,})")
ax.set_xticks(xt)
ax.set_xticklabels(xl)
ax.tick_params(axis="x", length=0)
ax.set_ylim(0, 115)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel("Expression variance\nexplained per gene (%)")
ax.text(-0.55, 106, "Median:", fontsize=6, color=DGREY, ha="right", va="bottom")
ax.legend(handles=[Patch(fc=C_SOURCE, label="Tissue source"),
                   Patch(fc=C_STUDY, label="Study")],
          loc="upper center", bbox_to_anchor=(0.36, 0.93), handlelength=0.9,
          handleheight=0.8, fontsize=6)
ax.set_xlim(-0.55, 2.55)
panel_label(ax, "A", x=-0.13, y=1.0)

# ======== b: matched null ========
ax = fig.add_subplot(row1[1])
r = nsum.loc[nsum["group"] == "published"].iloc[0]
obs, nmean, nsd, nset = float(r["observed_rate"]), float(r["null_mean"]), float(r["null_sd"]), int(r["n_sets"])
z = (obs - nmean) / nsd
vals, counts = np.unique(np.round(nrate["rate"].to_numpy(), 6), return_counts=True)
ax.bar(vals * 100, counts, width=0.75, color=LGREY, edgecolor="none")
ax.axvline(obs * 100, color=CL_COLORS["replicated"], lw=1.2)
ymax = counts.max() * 1.18
ax.set_ylim(0, ymax)
ax.text(obs * 100 - 1.0, ymax * 0.97, f"Reported markers\n{pct(obs)} (36/99)", ha="right",
        va="top", fontsize=6, color=CL_COLORS["replicated"])
ax.text(nmean * 100 + 4.5, ymax * 0.62,
        f"Matched random gene sets\n(n = {nset:,}): {pct(nmean)} ± {pct(nsd)}",
        ha="left", va="top", fontsize=6, color=DGREY)
ax.text(obs * 100 - 1.0, ymax * 0.70, f"z = {z:.1f}", ha="right", va="top", fontsize=6)
ax.set_xlim(0, 45)
ax.set_xlabel("Marker tests replicated in the reported direction (%)")
ax.set_ylabel("Random gene sets")
panel_label(ax, "B", x=-0.13, y=1.0)
nm.add("b", "observed_rate_published", f"{obs:.4f}", T_NS, "group=published")
nm.add("b", "null_mean_sd_published", f"{nmean:.4f} ± {nsd:.4f}", T_NS, f"n_sets={nset}")
nm.add("b", "z_vs_null", f"{z:.2f}", T_NS, "(obs-mean)/sd as shown")
nm.add("b", "null_rates_n", len(nrate), T_NR, "bar-plot data (counts of discrete values)")


# ======== stacked-bar helper ========
def stacked(ax, rows, labels, gap_after=None):
    xs = []
    x = 0.0
    for i, (_, rr) in enumerate(rows.iterrows()):
        xs.append(x)
        bottom = 0.0
        n = int(rr["n_rows"])
        for k in KEYS:
            c = int(rr[COLMAP[k]])
            v = c / n * 100
            ax.bar([x], [v], bottom=bottom, width=0.7, color=CL_COLORS[k],
                   edgecolor="white", lw=0.4)
            if v >= 6:
                ax.text(x, bottom + v / 2, f"{c}", ha="center", va="center", fontsize=6,
                        color=CL_TEXT[k])
            bottom += v
        ax.text(x, 101.5, f"{n}", ha="center", va="bottom", fontsize=6, color=DGREY)
        x += 1.0 + (0.45 if gap_after is not None and i == gap_after else 0.0)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, 110)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(xs[0] - 0.6, xs[-1] + 0.6)
    return xs


# ======== c: All + by comparison ========
ax = fig.add_subplot(row2[0])
crow = bs[bs["stratum"] == "contrast"].set_index("level").loc[["BM-UC", "AD-UC", "BM-AD"]].reset_index()
allrow = pd.DataFrame([{"level": "All", "n_rows": N_ALL, "claimed": cnt["replicated"],
                        "reversed": cnt["reversed"], "not_replicated": cnt["not_replicated"]}])
crows = pd.concat([allrow, crow[allrow.columns]], ignore_index=True)
assert int(crow["n_rows"].sum()) == N_ALL
stacked(ax, crows, ["All", f"BM{DASH}UC", f"AD{DASH}UC", f"BM{DASH}AD"], gap_after=0)
ax.set_ylabel("Marker tests (%)")
ax.text(-0.6, 108.5, "n =", fontsize=6, color=DGREY, ha="right", va="bottom")
panel_label(ax, "C", x=-0.30, y=1.0)
for _, rr in crows.iterrows():
    for k in KEYS:
        nm.add("c", f"{k}_{rr['level']}", int(rr[COLMAP[k]]),
               T_TC if rr["level"] == "All" else T_BS, f"level={rr['level']}, n={int(rr['n_rows'])}")

# ======== d: by original assay ========
ax = fig.add_subplot(row2[1])
ar = bs[bs["stratum"] == "assay_orig"].copy()
MAIN = ["RNAseq", "array", "qPCR", "protein"]
other = ar[~ar["level"].isin(MAIN)]
orow = pd.DataFrame([{"level": "Other", "n_rows": int(other["n_rows"].sum()),
                      "claimed": int(other["claimed"].sum()),
                      "reversed": int(other["reversed"].sum()),
                      "not_replicated": int(other["not_replicated"].sum())}])
arows = pd.concat([ar.set_index("level").loc[MAIN].reset_index()[orow.columns], orow],
                  ignore_index=True)
assert int(arows["n_rows"].sum()) == N_ALL
stacked(ax, arows, ["RNA-seq", "Array", "qPCR", "Protein", "Other"])
ax.set_ylabel("Marker tests (%)")
panel_label(ax, "D", x=-0.25, y=1.0)
for _, rr in arows.iterrows():
    for k in KEYS:
        nm.add("d", f"{k}_{rr['level']}", int(rr[COLMAP[k]]), T_BS,
               f"assay={rr['level']}, n={int(rr['n_rows'])}"
               + ("; caption Other = transcription-factor enrichment (2 claims; 4 tests)"
                  " and meta-analysis-derived claim (1 testable)"
                  if rr["level"] == "Other" else ""))

# ======== e: direction consistency (dot strip) ========
ax = fig.add_subplot(row2[2])
order = ["BM-UC", "AD-UC", "BM-AD"]
for i, c in enumerate(order):
    d = f2.loc[f2["contrast"] == c]
    v = d["frac_same_direction"].to_numpy() * 100
    key = np.round(v, 0)
    for val in np.unique(key):
        m = int((key == val).sum())
        step = min(0.045, 0.72 / max(m, 1))
        xs = i + (np.arange(m) - (m - 1) / 2) * step
        ax.scatter(xs, np.full(m, val), s=5.5, c="#4D4D4D", linewidths=0, zorder=3)
    med = float(np.median(v))
    ax.plot([i - 0.38, i + 0.38], [med, med], color="black", lw=1.1, zorder=4)
    kmed = float(np.median(d["K"]))
    ax.text(i, -24, f"{kmed:g}", ha="center", va="top", fontsize=6, color=DGREY, clip_on=False)
    nm.add("e", f"n_units_{c}", len(v), T_F2, "")
    nm.add("e", f"median_pct_same_{c}", f"{med:.1f}", T_F2, "")
    nm.add("e", f"median_K_{c}", f"{kmed:g}", T_F2, "median number of studies available per test")
ax.text(-0.62, -24, "Studies per test\n(median)", ha="right", va="top", fontsize=6,
        color=DGREY, clip_on=False)
ax.set_xticks(range(3))
ax.set_xticklabels([f"BM{DASH}UC", f"AD{DASH}UC", f"BM{DASH}AD"])
ax.tick_params(axis="x", length=0)
ax.set_ylim(-5, 105)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xlim(-0.6, 2.6)
ax.set_ylabel("Studies agreeing with\nreported direction (%)")
panel_label(ax, "E", x=-0.30, y=1.0)

# ---- shared legend for c–d (above row 2) ----
handles = [Patch(fc=CL_COLORS[k], ec="none", label=CL_LEGEND[k]) for k in KEYS]
fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.43, 0.505), ncol=3,
           frameon=False, handlelength=1.0, handleheight=0.8, columnspacing=1.4, fontsize=6)

save_fig(fig, "Fig2_replication")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "Fig2_replication.png"),
           os.path.join(FIG, "qc", "Fig2_replication_deutan.png"))
print("done")

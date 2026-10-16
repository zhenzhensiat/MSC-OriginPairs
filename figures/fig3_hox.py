"""fig3_hox.py —— Figure 3 (HOX family).

Panels: a HOX replication vs its own matched null distribution; b the 269 gene
families ranked by z (HOX highlighted; rank range across five random seeds);
c direction heat map by HOXA-D cluster and comparison; d effect-size comparison
against replicated published markers.

Data:
  data/fig3a_null_rates_hox_v2.tsv, data/null_summary_v2.tsv,
  data/hox_transfer_v2.tsv, data/marker_transfer_v2.tsv,
  data/family_ranking_v2.tsv + family_seed_summary_v2.tsv,
  data/marker_measurability_v2.tsv
Colours: source colours AD/BM/UC; HOX highlight purple #4B3F72; observation
line purple. Python performs counting only; statistics are asserted against the
result tables.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (CL_COLORS, SRC_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

FIG_NAME = "Fig3"
nm = Numbers(FIG_NAME)

T_NR = "data/fig3a_null_rates_hox_v2.tsv"
T_NS = "data/null_summary_v2.tsv"
T_HT = "data/hox_transfer_v2.tsv"
T_MT = "data/marker_transfer_v2.tsv"
T_FR = "data/family_ranking_v2.tsv"
T_FS = "data/family_seed_summary_v2.tsv"
T_MM = "data/marker_measurability_v2.tsv"

nrate = read_table(T_NR)
nsum = read_table(T_NS)
ht = read_table(T_HT)
mt = read_table(T_MT)
fr = read_table(T_FR)
fs = read_table(T_FS)
mm = read_table(T_MM)

C_AD, C_BM, C_UC = SRC_COLORS["AD"], SRC_COLORS["BM"], SRC_COLORS["UC"]
C_HOX = "#4B3F72"
C_OBS = CL_COLORS["replicated"]
CONTRASTS = ["BM-UC", "AD-UC", "BM-AD"]

# =====================================================================
# 0. Hard checks
# =====================================================================
r = nsum.loc[nsum["group"] == "hox"].iloc[0]
obs, nmean, nsd, nset = (float(r["observed_rate"]), float(r["null_mean"]),
                         float(r["null_sd"]), int(r["n_sets"]))
pval = float(r["empirical_p"])
z_obs = (obs - nmean) / nsd
assert len(nrate) == 1000 and abs(obs - 0.196581196581197) < 1e-9
assert abs(nmean - 0.0673760683760684) < 1e-9 and abs(nsd - 0.0283398818092925) < 1e-9
assert abs(pval - 0.000999000999000999) < 1e-12
n_tr = int(ht["transfer"].astype(str).str.upper().eq("TRUE").sum())
n_units = len(ht)
assert n_units == 117 and n_tr == 23
assert int(mt["transfer"].astype(str).str.upper().eq("TRUE").sum()) == 36

hox_row = fr.loc[fr["family"] == "HOX_39"].iloc[0]
assert int(hox_row["rank_z"]) == 13 and abs(float(hox_row["z"]) - 4.70633531562324) < 1e-9
assert int((fr["included"].astype(str).str.upper().eq("TRUE")).sum()) == 269
seed_ranks = sorted(int(v) for v in fs["rank_z"])
assert seed_ranks == [9, 13, 14, 14, 15], seed_ranks

med_hox = float(mm.loc[mm["group"] == "HOX", "median_abs_logFC_transferred"].iloc[0])
med_pub = float(mm.loc[mm["group"] == "published", "median_abs_logFC_transferred"].iloc[0])
assert abs(med_hox - 1.0735175535444) < 1e-9 and abs(med_pub - 2.23384453273294) < 1e-9
lfc_hox = ht.loc[ht["transfer"].astype(str).str.upper().eq("TRUE"), "median_logFC"].abs()
lfc_pub = mt.loc[mt["transfer"].astype(str).str.upper().eq("TRUE"), "median_logFC"].abs()
assert abs(float(lfc_hox.median()) - med_hox) < 1e-9
assert abs(float(lfc_pub.median()) - med_pub) < 1e-9

# =====================================================================
# Layout
# =====================================================================
fig = plt.figure(figsize=(180 * MM, 176 * MM))
outer = fig.add_gridspec(3, 1, height_ratios=[46, 15, 104], hspace=0.0,
                         left=0.088, right=0.992, top=0.965, bottom=0.03)
row1 = outer[0].subgridspec(1, 2, width_ratios=[0.96, 1.16], wspace=0.34)
row2 = outer[2].subgridspec(1, 2, width_ratios=[1.0, 1.28], wspace=0.30)

# =====================================================================
# 3a HOX vs self-matched null distribution (discrete k/117 bars; annotation moved to left blank area)
# =====================================================================
axa = fig.add_subplot(row1[0])
vals_k = nrate["rate"].to_numpy() * 117
assert np.max(np.abs(vals_k - np.round(vals_k))) < 1e-6, "replication rates should take discrete k/117 values"
ks = np.round(vals_k).astype(int)
kk = np.arange(ks.min(), ks.max() + 1)
cnts = np.bincount(ks, minlength=ks.max() + 1)[kk]
axa.bar(kk / 117 * 100, cnts, width=(100 / 117) * 0.88, color="#C9C9C9",
        edgecolor="white", lw=0.3)
axa.axvline(obs * 100, color=C_OBS, lw=1.2)
ymax = axa.get_ylim()[1]
# Three-part annotation: in blank area above-left of the observed line (x≈15.5–19), left-aligned, narrow width; gap from the observed line, does not cover bars
xa = 15.45 / 24.0
axa.text(xa, 0.97, f"HOX\n{obs * 100:.1f}%\n(23/117)", transform=axa.transAxes,
         ha="left", va="top", fontsize=6, color=C_OBS, linespacing=1.35, zorder=5)
axa.text(xa, 0.755, f"z = {z_obs:.1f},\np = {pval:.3f}", transform=axa.transAxes,
         ha="left", va="top", fontsize=6, linespacing=1.35, zorder=5)
axa.text(xa, 0.60, "Matched\nrandom\nsets\n"
         f"{nmean * 100:.1f}%±{nsd * 100:.1f}%\n(n={nset:,})",
         transform=axa.transAxes, ha="left", va="top", fontsize=6, color=DGREY,
         linespacing=1.35, zorder=5)
axa.set_xlim(0, 24)
axa.set_xticks([0, 5, 10, 15, 20])
axa.set_xlabel("HOX units replicated (%)", fontsize=6.5, labelpad=1)
axa.set_ylabel("Random gene sets", fontsize=6.5)
panel_label(axa, "A", x=-0.20, y=1.03)
nm.add("a", "observed", f"{obs * 100:.1f}% (23/117)", T_NS, "group=hox")
nm.add("a", "null_mean_sd", f"{nmean * 100:.1f}% ± {nsd * 100:.1f}%", T_NS, f"n_sets={nset}")
nm.add("a", "z", f"{z_obs:.2f}", T_NS, "(obs-mean)/sd")
nm.add("a", "empirical_p", f"{pval:.3f}", T_NS, "1000 null distributions")
nm.add("a", "null_n", len(nrate), T_NR, "bar-plot data")
nm.add("a", "discrete_bars", "bars drawn at discrete k/117 values", T_NR, "consistent with Fig2b convention")
nm.add("a", "text_sensitivity_no_gse38947", "24/117 = 20.5%",
       "data/sensitivity_no_gse38947_v2.tsv", "cited in main text §3 (GSE38947-excluded sensitivity; verified 24/117)")

# =====================================================================
# 3b Ranking of 269 families (HOX purple point + visible seed-range segment)
# =====================================================================
axb = fig.add_subplot(row1[1])
inc = fr[fr["included"].astype(str).str.upper().eq("TRUE")].copy()
axb.scatter(inc["rank_z"], inc["z"], s=4, c="#8C8C8C", linewidths=0, alpha=0.85, zorder=2)
hx = int(hox_row["rank_z"])
z_hox = float(hox_row["z"])
axb.plot([min(seed_ranks), max(seed_ranks)], [z_hox, z_hox], color=C_HOX, lw=1.6,
         solid_capstyle="butt", zorder=5)
for xe in (min(seed_ranks), max(seed_ranks)):
    axb.plot([xe, xe], [z_hox - 0.14, z_hox + 0.14], color=C_HOX, lw=0.9, zorder=5)
axb.plot([hx], [z_hox], marker="o", ms=4.2, color=C_HOX, zorder=6)
axb.annotate(f"HOX ({hx}th; seeds {min(seed_ranks)}–{max(seed_ranks)})",
             xy=(hx + 2, z_hox - 0.18), xytext=(58, 6.35),
             fontsize=6, color=C_HOX,
             arrowprops=dict(arrowstyle="-", lw=0.6, color=C_HOX))
axb.set_xlim(0, 276)
z_lo, z_hi = float(inc["z"].min()), float(inc["z"].max())
axb.set_ylim(z_lo - 0.3, z_hi + 0.3)
axb.set_xlabel("Family rank (z vs own matched null)", fontsize=6.5, labelpad=1)
axb.set_ylabel("z score", fontsize=6.5)
axb.text(276, z_lo - 0.15, f"n = {int(inc.shape[0])} families", ha="right", va="bottom",
         fontsize=6, color=DGREY)
panel_label(axb, "B", x=-0.11, y=1.03)
nm.add("b", "n_families", 269, T_FR, "included=TRUE; ylim covers full set")
nm.add("b", "hox_rank", 13, T_FR, "HOX_39; primary seed")
nm.add("b", "hox_seed_range", "9–15", T_FS, "5 seed rank_z: 9/13/14/14/15 (shown as horizontal segment)")
nm.add("b", "hox_z", f"{z_hox:.3f}", T_FR, "")
nm.add("b", "hox_z_seeds", "4.23–5.17", T_FS, "z range across seeds")
nm.add("b", "ylim", f"{z_lo - 0.3:.2f} – {z_hi + 0.3:.2f}", T_FR, "covers all 269 families")
nm.add("b", "text_top6pct", "9–15 / 269 = top 6%", T_FS, "cited in main text §3 (15/269 = 5.6%)")

# =====================================================================
# 3c Direction heatmap (39 genes × 3 contrasts; clustered by A–D)
# =====================================================================
axc = fig.add_subplot(row2[0])


def cluster_of(g):
    return {"A": 0, "B": 1, "C": 2, "D": 3}[g[3]]


def num_of(g):
    import re
    return int(re.search(r"\d+", g).group())


genes = sorted(ht["gene"].unique(), key=lambda g: (cluster_of(g), num_of(g)))
CLUSTERS = ["HOXA", "HOXB", "HOXC", "HOXD"]
rows_g = []
gap_idx = []
for cl in CLUSTERS:
    gs = [g for g in genes if g.startswith(cl)]
    rows_g.extend(gs)
    gap_idx.append(len(rows_g))
rows_g2 = list(rows_g)
gap_idx = gap_idx[:-1]
# Build matrix by inserting blank rows
rowpos = {}
pos = 0
for i, g in enumerate(rows_g2):
    if i in gap_idx:
        pos += 1
    rowpos[g] = pos
    pos += 1
NR = pos

cell_code = np.full((NR, 3), -1, dtype=int)
rep = np.zeros((NR, 3), dtype=bool)
for _, rr in ht.iterrows():
    i, j = rowpos[rr["gene"]], CONTRASTS.index(rr["contrast"])
    lfc = float(rr["median_logFC"])
    first, second = rr["contrast"].split("-")
    higher = first if lfc > 0 else second
    cell_code[i, j] = {"AD": 0, "BM": 1, "UC": 2}[higher]
    rep[i, j] = str(rr["transfer"]).upper() == "TRUE"
axc.set_ylim(NR, 0)
axc.set_xlim(0, 4.1)
for i in range(NR):
    for j in range(3):
        if cell_code[i, j] >= 0:
            colr = [C_AD, C_BM, C_UC][cell_code[i, j]]
            axc.add_patch(Rectangle((j, i), 1, 1, fc=colr, ec="none",
                                    alpha=(1.0 if rep[i, j] else 0.30), zorder=2))
        if rep[i, j]:
            axc.plot(j + 0.5, i + 0.5, marker="o", ms=1.6, mec="black", mfc="black", zorder=4)
axc.set_xticks([0.5, 1.5, 2.5])
axc.set_xticklabels(["BM" + "–" + "UC", "AD–UC", "BM–AD"], fontsize=6)
axc.xaxis.set_ticks_position("top")
axc.set_yticks([rowpos[g] + 0.5 for g in rows_g2])
axc.set_yticklabels(rows_g2, fontsize=6)
axc.tick_params(length=0, pad=1.5)
for s in axc.spines.values():
    s.set_visible(False)
# Cluster labels (right side)
for cl in CLUSTERS:
    gs = [g for g in rows_g2 if g.startswith(cl)]
    cy = np.mean([rowpos[g] + 0.5 for g in gs])
    axc.text(3.25, cy, cl, fontsize=6, color=DGREY, ha="left", va="center")
axc.set_xlim(0, 4.1)
panel_label(axc, "C", x=-0.30, y=1.075)
nm.add("c", "n_units", "39 genes × 3 contrasts = 117", T_HT, "all testable")
nm.add("c", "n_replicated", 23, T_HT, "dot = replicated in majority of studies")
for c in CONTRASTS:
    n_c = int(ht[(ht["contrast"] == c) & ht["transfer"]].shape[0])
    nm.add("c", f"replicated_{c}", n_c, T_HT, "")
nm.add("c", "clusters", "HOXA 11 / HOXB 10 / HOXC 9 / HOXD 9", T_HT, "rows ordered by cluster, whitespace between groups")
nm.add("c", "cell_color", "cell = source with higher median expression", T_HT,
       "direction=sign(median_logFC), logFC=former minus latter")

# =====================================================================
# 3d Effect-size comparison (replication units |log2FC|; panel shrunk, explanation moved to caption)
# =====================================================================
axd = fig.add_subplot(row2[1].subgridspec(2, 1, height_ratios=[0.72, 0.28], hspace=0)[0])
groups = [("HOX", lfc_hox.to_numpy(), C_HOX), ("Reported markers", lfc_pub.to_numpy(), "#333333")]
rng = np.random.default_rng(3)
for yi, (lab, vals_g, colr) in enumerate(groups):
    y0 = 1 - yi
    axd.add_patch(Rectangle((0, y0 - 0.30), 6.0, 0.60, fc="#F5F5F5", ec="none", zorder=1))
    jit = rng.uniform(-0.19, 0.19, len(vals_g))
    axd.scatter(vals_g, y0 + jit, s=6, c=colr, linewidths=0, alpha=0.85, zorder=3)
    med = float(np.median(vals_g))
    axd.plot([med, med], [y0 - 0.30, y0 + 0.30], color="black", lw=1.1, zorder=4)
    axd.text(med + 0.08, y0 + 0.36, f"median {med:.2f}", fontsize=6, va="bottom",
             ha="left")
    axd.text(5.95, y0 + 0.40, f"n = {len(vals_g)}", fontsize=6, va="bottom", ha="right",
             color=DGREY)
axd.set_yticks([1, 0])
axd.set_yticklabels(["HOX", "Reported\nmarkers"], fontsize=6)
axd.tick_params(length=0, pad=2)
axd.set_ylim(-0.55, 1.62)
axd.set_xlim(0, 6.0)
axd.set_xticks([0, 1, 2, 3, 4, 5, 6])
axd.set_xlabel("|log2 fold change| (replicated units)", fontsize=6.5, labelpad=1)
panel_label(axd, "D", x=-0.16, y=1.05)
nm.add("d", "median_hox", f"{med_hox:.2f}", T_MM, "median |log2FC| of replication units")
nm.add("d", "median_published", f"{med_pub:.2f}", T_MM, "median |log2FC| of replication units")
nm.add("d", "n_hox", len(lfc_hox), T_HT, "transfer=TRUE")
nm.add("d", "n_published", len(lfc_pub), T_MT, "transfer=TRUE")
nm.add("d", "note_scope", "replication units only; explanatory text moved to caption", T_MM, "3d panel shrunk")

# =====================================================================
save_fig(fig, "Fig3_hox")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "Fig3_hox.png"), os.path.join(FIG, "qc", "Fig3_hox_deutan.png"))
print("done")

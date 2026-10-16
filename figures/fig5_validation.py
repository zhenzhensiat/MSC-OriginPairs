"""fig5_validation.py —— Figure 5 (external validation; core figure).

Panels: a sample-by-rule vote heat map (30 rules x 176 samples ordered by
source and study; holds = source colour, fails = #EEEEEE, not computable =
#9E9E9E; three top annotation bars; 180 mm full width); b per-study accuracy;
c label-permutation null distribution; d forest plot of sensitivity and
subgroup analyses.

Data: fig5a_vote_matrix_v2.tsv,
confirm_predictions_v2.tsv, confirm_null_v2.tsv, confirm_results_v2.tsv,
confirm_sensitivity_v2.tsv, confirm_strata_curated_v2.tsv,
fig5d_platform_accuracy_v2.tsv.
Python performs counting only; statistics come from the result tables and are
asserted in-script.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (SRC_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

FIG_NAME = "Fig5"
nm = Numbers(FIG_NAME)

T_VM = "data/fig5a_vote_matrix_v2.tsv"
T_CP = "data/confirm_predictions_v2.tsv"
T_NULL = "data/confirm_null_v2.tsv"
T_RES = "data/confirm_results_v2.tsv"
T_SENS = "data/confirm_sensitivity_v2.tsv"
T_STRATA = "data/confirm_strata_curated_v2.tsv"
T_PLAT = "data/fig5d_platform_accuracy_v2.tsv"

vm = read_table(T_VM)
cp = read_table(T_CP)
nn = read_table(T_NULL)
res = read_table(T_RES)
sens = read_table(T_SENS)
strata = read_table(T_STRATA)
plat = read_table(T_PLAT)

C_AD, C_BM, C_UC = SRC_COLORS["AD"], SRC_COLORS["BM"], SRC_COLORS["UC"]
C_NOT, C_MISS, C_WHITE = "#EEEEEE", "#9E9E9E", "#FFFFFF"
CELL_COLORS = [C_AD, C_BM, C_UC, C_NOT, C_MISS, C_WHITE]
CLS_ORD = {"AD": 0, "BM": 1, "UC": 2}
THR = 0.2494  # low-confidence threshold (τ; fixed value from utils_models.R choose_tau)

# =====================================================================
# 0. Hard checks on the result tables
# =====================================================================
assert len(vm) == 176 * 30 and vm["gsm"].nunique() == 176, "vote matrix dimensions"
assert dict(vm["state"].value_counts()) == {"not": 3082, "holds": 1866, "missing": 332}, "state counts"
assert vm.groupby("class")["rule_rank"].max().to_dict() == {"AD": 15, "BM": 10, "UC": 5}, "rule counts"
assert (vm["class_direction"] == "gene1>gene2").all(), "direction assumption"

rules = (vm[["class", "rule_rank", "gene1", "gene2"]].drop_duplicates()
         .assign(ord=lambda d: d["class"].map(CLS_ORD))
         .sort_values(["ord", "rule_rank"]).reset_index(drop=True))
assert len(rules) == 30

meta = cp[["gsm", "study", "true", "pred", "margin", "donor"]].copy()
assert len(meta) == 176 and set(meta["gsm"]) == set(vm["gsm"]), "sample sets identical"
assert dict(meta["true"].value_counts()) == {"UC": 93, "BM": 55, "AD": 28}, "source counts"
assert meta["study"].nunique() == 22, "study count"
n_nocall = int(meta["pred"].isna().sum())
assert n_nocall == 7 and set(meta.loc[meta["pred"].isna(), "study"]) == {"GSE24871"}, "no call"
mrg = pd.to_numeric(meta["margin"], errors="coerce")
assert int((mrg < THR).sum()) == 18, "low-confidence count"
assert (meta.loc[mrg < THR, "pred"] == meta.loc[mrg < THR, "true"]).all(), "all low-confidence calls correct"
assert int(cp["rejected"].astype(str).str.upper().eq("TRUE").sum()) == 18, "rejected convention consistent"

meta["cls_ord"] = meta["true"].map(CLS_ORD)
meta["study_num"] = meta["study"].str.replace("GSE", "", regex=False).astype(int)
meta = meta.sort_values(["cls_ord", "study_num", "gsm"]).reset_index(drop=True)

# Per-study accuracy (no calls counted as errors per the fixed convention)
per = (meta.assign(ok=(meta["pred"] == meta["true"]))
       .groupby("study", sort=False)
       .agg(n=("gsm", "size"), acc=("ok", "mean"), cls=("true", "first"), ordn=("cls_ord", "first"),
            sn=("study_num", "first")).reset_index())
assert per["n"].sum() == 176 and len(per) == 22
assert per.loc[per["study"] == "GSE24871", "acc"].iloc[0] == 0.0
assert (per.loc[per["study"] != "GSE24871", "acc"] == 1.0).all()

# Values for each forest-plot row
def pick(df, name, metric):
    v = df.loc[(df["analysis"] == name) & (df["metric"] == metric), "value"]
    assert len(v) == 1, (name, metric)
    return float(v.iloc[0])

ba_main = float(res.loc[res["metric"] == "balanced_accuracy", "value"].iloc[0])
ci_lo = float(res.loc[res["metric"] == "ba_ci_lo", "value"].iloc[0])
ci_hi = float(res.loc[res["metric"] == "ba_ci_hi", "value"].iloc[0])
ba_excl = pick(sens, "excl_GSE24871", "balanced_accuracy")
ba_loo_min = pick(sens, "leave_one_cluster_out", "ba_min")
ba_donor = pick(sens, "donor_level", "balanced_accuracy")
assert abs(ba_main - 0.957575757575758) < 1e-12 and abs(ci_lo - 0.857142857142857) < 1e-12
assert ba_excl == 1.0 and abs(ba_loo_min - 0.946969696969697) < 1e-12
assert abs(ba_donor - 0.941666666666667) < 1e-12

# Donor-level unit count (convention: majority vote after grouping by study × donor; consistent with R utils_metrics.donor_level)
donor_rows = []
for _, d in meta.groupby(["study", "donor"], sort=False):
    vc = d["pred"].value_counts()
    top = vc[vc == vc.max()]
    donor_rows.append((d["true"].iloc[0], top.index[0] if len(top) == 1 else np.nan))
dn = len(donor_rows)
assert dn == 128, dn

# Fold with the minimum under leave-one-team-out: n = 176 - cluster sample count (GSE111495, n=11 → 165)
loo = []
for cl in sorted(meta["study"].unique()):
    sub = meta[meta["study"] != cl]
    rec = [((sub.loc[sub["true"] == c, "pred"] == c).mean()) for c in ("AD", "BM", "UC")]
    loo.append((float(np.mean(rec)), cl, len(sub)))
loo.sort(key=lambda t: (t[0], t[1]))
assert abs(loo[0][0] - ba_loo_min) < 1e-12 and loo[0][2] == 165, loo[0]

med = strata[strata["stratum"] == "medium_cur"].set_index("level")
FBS, HPL, XENO, MEDUNK = (float(med.loc["FBS", "ba"]), float(med.loc["hPL", "ba"]),
                          float(med.loc["xeno_free", "ba"]), float(med.loc["unknown", "ba"]))
assert (FBS, HPL, XENO) == (1.0, 1.0, 1.0) and abs(MEDUNK - 0.933333333333333) < 1e-12
PL = plat.set_index("platform")
BA_ARR, N_ARR = float(PL.loc["microarray", "ba"]), int(PL.loc["microarray", "n"])
BA_RNA, N_RNA = float(PL.loc["RNA-seq", "ba"]), int(PL.loc["RNA-seq", "n"])
assert (BA_ARR, N_ARR, BA_RNA, N_RNA) == (0.933333333333333, 78, 1.0, 98)

# Permutation null distribution
perm = nn["perm_ba"].to_numpy() * 100
assert len(perm) == 1000
perm_mean = float(nn["perm_ba"].mean() * 100)
n_ge = int((nn["perm_ba"] >= ba_main).sum())
assert n_ge == 0

# =====================================================================
# 1. 5a matrix and column order
# =====================================================================
col_of = {g: j for j, g in enumerate(meta["gsm"])}
row_of = {(r["class"], r["rule_rank"]): i for i, r in rules.iterrows()}
Cm = np.full((30, 176), -1, dtype=int)
for rec in vm[["gsm", "class", "rule_rank", "state"]].to_dict("records"):
    i = row_of[(rec["class"], rec["rule_rank"])]
    j = col_of[rec["gsm"]]
    s = rec["state"]
    Cm[i, j] = CLS_ORD[rec["class"]] if s == "holds" else (3 if s == "not" else 4)
assert (Cm >= 0).all()

# Add inter-group whitespace (rows: 1 blank row after each of AD|BM|UC; columns: 3 blank columns between source blocks)
def with_gaps(arr, row_bounds, col_bounds, gap_rows=1, gap_cols=3):
    rows_a, cols_a = [], []
    row_of_new = {}
    ri = 0
    for i in range(arr.shape[0]):
        if i in row_bounds:
            for _ in range(gap_rows):
                rows_a.append(np.full(arr.shape[1], 5))
                ri += 1
        rows_a.append(arr[i])
        row_of_new[i] = ri
        ri += 1
    top = np.vstack(rows_a)
    cols_b = []
    ci = 0
    col_of_new = {}
    for j in range(top.shape[1]):
        if j in col_bounds:
            for _ in range(gap_cols):
                cols_b.append(np.full(top.shape[0], 5))
                ci += 1
        cols_b.append(top[:, j])
        col_of_new[j] = ci
        ci += 1
    return np.column_stack(cols_b), row_of_new, col_of_new

row_bounds = {15, 25}                     # before BM, UC blocks
col_bounds = {28, 83}                     # before BM, UC blocks (AD 28, BM 55)
Cg, row_new, col_new = with_gaps(Cm, row_bounds, col_bounds)
NR, NC = Cg.shape                         # 32 × 182
aug = np.full((3, NC), 5, dtype=int)      # annotation bar (A: labeled source / B: predicted source)
cmargs = np.full(NC, np.nan)              # annotation bar C: margin
study_aug = np.full(NC, "", dtype=object) # study of each column (including blank columns), used for separators and ranges
for g, j in col_of.items():
    jj = col_new[j]
    m = meta.loc[meta["gsm"] == g].iloc[0]
    aug[0, jj] = CLS_ORD[m["true"]]
    if pd.notna(m["pred"]):
        aug[1, jj] = CLS_ORD[m["pred"]]
    if pd.notna(m["margin"]):
        cmargs[jj] = float(m["margin"])
    study_aug[jj] = m["study"]

# Study boundaries (column coordinates, for white separator lines)
sep_x = []
for jj in range(1, NC):
    if study_aug[jj] != "" and study_aug[jj] != study_aug[jj - 1] and study_aug[jj - 1] != "":
        sep_x.append(jj)

# no-call runs (× drawn on annotation bar B)
runs, j = [], 0
while j < NC:
    if study_aug[j] == "GSE24871" and j > 0 and study_aug[j - 1] != "GSE24871":
        k = j
        while k < NC and study_aug[k] == "GSE24871":
            k += 1
        runs.append((j, k))
        j = k
    else:
        j += 1

# =====================================================================
# 2. Layout
# =====================================================================
fig = plt.figure(figsize=(180 * MM, 234 * MM))
outer = fig.add_gridspec(5, 1, height_ratios=[100, 6, 42, 26, 58], hspace=0.0,
                         left=0.115, right=0.992, top=0.975, bottom=0.035)
sec1 = outer[0].subgridspec(5, 1, height_ratios=[2.6, 2.6, 2.6, 68, 7.2], hspace=0.12)
sec2 = outer[2].subgridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.24)

cmap6 = ListedColormap(CELL_COLORS)

# ---- 5a annotation bar ----
strip_axes = []
strip_labels = ["Reported source", "Predicted source", "Confidence"]
for si in range(3):
    axs = fig.add_subplot(sec1[si])
    strip_axes.append(axs)
    if si < 2:
        axs.imshow(aug[si:si + 1], cmap=cmap6, vmin=-0.5, vmax=5.5, aspect="auto",
                   interpolation="nearest", extent=(0, NC, 1, 0))
    else:
        cmc = plt.get_cmap("Greys").copy()
        cmc.set_bad("white")
        arr = np.ma.masked_invalid(cmargs.reshape(1, -1))
        axs.imshow(arr, cmap=cmc, vmin=0, vmax=1, aspect="auto",
                   interpolation="nearest", extent=(0, NC, 1, 0))
    axs.set_xlim(0, NC)
    for s in axs.spines.values():
        s.set_visible(False)
    axs.set_xticks([])
    axs.set_yticks([])
    axs.text(-0.008, 0.5, strip_labels[si], transform=axs.transAxes, ha="right",
             va="center", fontsize=6, color=DGREY)
    for x in sep_x:
        axs.axvline(x, color="white", lw=0.4, zorder=4)
for (a0, a1) in runs:
    strip_axes[1].text((a0 + a1) / 2, 0.5, "×", ha="center", va="center",
                       fontsize=7, color=DGREY, zorder=6)
nm.add("a", "strips", "Reported/Predicted source + Confidence (Greys, 0-1)", T_CP,
       "annotation bar A colored by true; B colored by pred, NA=white ×; C is margin grayscale")

# ---- 5a heatmap ----
axh = fig.add_subplot(sec1[3])
axh.imshow(Cg, cmap=cmap6, vmin=-0.5, vmax=5.5, aspect="auto",
           interpolation="nearest", extent=(0, NC, NR, 0))
for x in sep_x:
    axh.plot([x, x], [0, NR], color="white", lw=0.5, zorder=5)   # heatmap area only
ypos = [row_new[i] for i in range(30)]
ylab = [f"{r['gene1']} > {r['gene2']}" for _, r in rules.iterrows()]
axh.set_yticks([p + 0.5 for p in ypos])
axh.set_yticklabels(ylab, fontsize=6)
axh.tick_params(axis="y", length=0, pad=1.5)
axh.set_xticks([])
axh.set_ylim(NR + 3.3, 0)
for s in axh.spines.values():
    s.set_visible(False)
# Block labels
for cls, txt in [("AD", "AD (n = 28)"), ("BM", "BM (n = 55)"), ("UC", "UC (n = 93)")]:
    nm.add("a", f"block_{cls}", txt, T_CP, "column-block sample count/source")
    cx = np.mean([col_new[j] for j in range(176) if meta.iloc[j]["cls_ord"] == CLS_ORD[cls]])
    axh.text(cx, NR + 1.6, txt, ha="center", va="top", fontsize=6, color=DGREY)
panel_label(strip_axes[0], "A", x=-0.028, y=1.5)
nm.add("a", "n_samples", 176, T_CP, "external validation cohort")
nm.add("a", "n_rules", "30 (AD 15 / BM 10 / UC 5)", T_VM, "grouped by source of the rule")
nm.add("a", "cell_states", "holds 1866 / not 3082 / missing 332", T_VM, "all 176×30 cells")
nm.add("a", "no_call_samples", "7 (GSE24871)", T_CP, "pred=NA; counted as errors per the fixed convention")
nm.add("a", "low_conf_calls", "18/169 (10.7%)", T_CP, f"margin<{THR}", )
nm.add("a", "threshold", THR, "R/utils_models.R choose_tau",
       "5th percentile of margins for correctly decided training cases (fixed)")

# ---- legend + confidence colorbar (sec1 row 5) ----
axl = fig.add_subplot(sec1[4])
axl.set_xlim(0, 180)
axl.set_ylim(0, 1)
axl.axis("off")


def _prect(ax, x, w, color):
    ax.add_patch(Rectangle((x, 0.30), w, 0.40, fc=color, ec="none", clip_on=False))


def _ptext(ax, x, s, **kw):
    ax.text(x, 0.5, s, fontsize=6, va="center", ha="left", **kw)


x = 0.0
_ptext(axl, x, "Rule holds:")
x += 18
for c, lab in [(C_AD, "AD"), (C_BM, "BM"), (C_UC, "UC")]:
    _prect(axl, x, 2.6, c)
    _ptext(axl, x + 3.4, lab)
    x += 14
_prect(axl, x, 2.6, C_NOT)
_ptext(axl, x + 3.4, "Rule does not hold")
x += 28
_prect(axl, x, 2.6, C_MISS)
_ptext(axl, x + 3.4, "Not computable")
x += 24
_ptext(axl, x, "×  No call")
x += 14
_ptext(axl, x, "Confidence:", color=DGREY)
x += 15
grad = np.linspace(0, 1, 256).reshape(1, -1)
axl.imshow(grad, cmap="Greys", aspect="auto", extent=(x, x + 34, 0.30, 0.70))
for t in (0, 0.5, 1):
    axl.plot([x + t * 34, x + t * 34], [0.24, 0.30], color=DGREY, lw=0.5)
    axl.text(x + t * 34, 0.12, f"{t:g}", fontsize=6, color=DGREY,
             ha="center", va="top")
xt = x + THR * 34
axl.plot([xt, xt], [0.22, 0.78], color="black", lw=0.8)
axl.text(xt, 0.86, f"{THR}", fontsize=6, ha="center", va="bottom")

# =====================================================================
# 3. 5b per-study accuracy
# =====================================================================
axb = fig.add_subplot(sec2[0])
order = per.sort_values(["ordn", "sn"]).reset_index(drop=True)
xs = np.arange(len(order))
for i, r in order.iterrows():
    if r["study"] == "GSE24871":
        axb.scatter(i, 0, marker="x", s=26, c=SRC_COLORS[r["cls"]], linewidths=1.0, zorder=4)
    else:
        axb.scatter(i, r["acc"] * 100, s=4 + r["n"] * 1.35, c=SRC_COLORS[r["cls"]],
                    linewidths=0, zorder=4)
bounds = [i for i in range(1, len(order)) if order.loc[i, "ordn"] != order.loc[i - 1, "ordn"]]
for b in bounds:
    axb.axvline(b - 0.5, color="#DDDDDD", lw=0.6, zorder=1)
axb.set_xticks(xs)
axb.set_xticklabels(order["study"], rotation=90, fontsize=6)
axb.set_xlim(-0.8, len(order) - 0.2)
axb.set_ylim(-2, 105)
axb.set_yticks([0, 25, 50, 75, 100])
axb.set_ylabel("Accuracy (%)")
panel_label(axb, "B", x=-0.075, y=1.02)
# Marker-size legend (blank area, lower left)
for (yy, nn_) in [(28, 3), (18, 10), (4, 30)]:
    axb.scatter(0.4, yy, s=4 + nn_ * 1.35, c="#808080", linewidths=0, zorder=3)
    axb.text(0.75, yy, f"n = {nn_}", fontsize=6, va="center", color=DGREY)
nm.add("b", "marker", "× = no call (GSE24871; counted as 7 errors)", T_CP, "")
for _, r in order.iterrows():
    nm.add("b", f"acc_{r['study']}", f"{r['acc'] * 100:.1f}% (n={r['n']})", T_CP,
           "per-study accuracy=mean(pred==true), no calls counted as errors")

# =====================================================================
# 4. 5c permutation null distribution
# =====================================================================
axc = fig.add_subplot(sec2[1])
bins = np.arange(0, 72.5, 2.5)
axc.hist(perm, bins=bins, color="#C9C9C9", edgecolor="white", lw=0.3)
axc.axvline(ba_main * 100, color="black", lw=1.1)
_ymax = axc.get_ylim()[1]
axc.set_ylim(0, _ymax * 1.38)                                   # raise to leave headroom
axc.text(ba_main * 100 - 1.5, _ymax * 1.34, f"Observed\n{ba_main * 100:.1f}%",
         ha="right", va="top", fontsize=6)
axc.text(2.5, _ymax * 1.34,
         f"Label permutations (n = {len(perm):,})\nmean {perm_mean:.1f}%\np < 0.001 (0/{len(perm):,} ≥ observed)",
         ha="left", va="top", fontsize=6, color=DGREY)
axc.set_xlim(0, 100)
axc.set_xlabel("Balanced accuracy (%)")
axc.set_ylabel("Permutations")
panel_label(axc, "C", x=-0.16, y=1.02)
nm.add("c", "n_perm", 1000, T_NULL, "number of label permutations")
nm.add("c", "perm_mean", f"{perm_mean:.1f}%", T_NULL, "permutation null mean")
nm.add("c", "observed", f"{ba_main * 100:.1f}%", T_RES, "observed BA (main analysis)")
nm.add("c", "perm_ge_obs", n_ge, T_NULL, "permuted BA ≥ observed (0 → p<0.001)")

# =====================================================================
# 5. 5d forest plot
# =====================================================================
sec3 = outer[4].subgridspec(1, 2, width_ratios=[41, 117], wspace=0.05)
axdl = fig.add_subplot(sec3[0])   # row-label column
axd = fig.add_subplot(sec3[1])
rows = [
    ("Main analysis", ba_main, 176, (ci_lo, ci_hi)),
    ("Excluding GSE24871", ba_excl, 169, None),
    ("Leave-one-group-out, min.", ba_loo_min, 165, None),
    ("Donor level", ba_donor, dn, None),
    ("Medium: FBS", FBS, int(med.loc["FBS", "n"]), None),
    ("Medium: hPL", HPL, int(med.loc["hPL", "n"]), None),
    ("Medium: xeno-free", XENO, int(med.loc["xeno_free", "n"]), None),
    ("Medium: not reported", MEDUNK, int(med.loc["unknown", "n"]), None),
    ("Platform: microarray", BA_ARR, N_ARR, None),
    ("Platform: RNA-seq", BA_RNA, N_RNA, None),
]
Y = np.arange(len(rows))[::-1]
axdl.set_xlim(0, 1)
axdl.set_ylim(-0.8, len(rows) - 0.2)
axdl.axis("off")
for y, (lab, v, n, ci) in zip(Y, rows):
    axdl.text(0.995, y, f"{lab} (n = {n})", fontsize=6, ha="right", va="center")
    axd.scatter(v * 100, y, s=13, c="black", zorder=4)
    if ci:
        axd.errorbar(v * 100, y, xerr=[[(v - ci[0]) * 100], [(ci[1] - v) * 100]],
                     fmt="none", ecolor="black", elinewidth=0.8, capsize=1.8, zorder=3)
    lx = (ci[1] * 100 + 0.5) if ci else (v * 100 + 0.5)         # right of point, same row; main analysis placed beyond CI right end
    axd.text(lx, y, f"{v * 100:.1f}", fontsize=6, va="center", ha="left", color=DGREY)
for sy in [6.5, 2.5]:
    axd.axhline(sy, color="#DDDDDD", lw=0.6, zorder=1)
axd.set_ylim(-0.8, len(rows) - 0.2)
axd.set_yticks([])
axd.tick_params(axis="y", length=0)
axd.set_xlim(84, 103.5)
axd.set_xticks([85, 90, 95, 100])
axd.set_xlabel("Balanced accuracy (%)")
panel_label(axdl, "D", x=-0.015, y=1.04)
for (lab, v, n, ci) in rows:
    if lab == "Main analysis":
        note = f"main analysis (with 95% CI {ci_lo * 100:.1f}-{ci_hi * 100:.1f})"
    elif lab == "Leave-one-group-out, min.":
        note = "post-hoc; lowest fold = n=165 after leaving out GSE111495 (n=11) (tied with GSE90149 fold)"
    elif lab == "Donor level":
        note = "post-hoc; 128 donor units (study×donor majority vote, convention=R utils_metrics.donor_level)"
    else:
        note = "post-hoc"
    if lab == "Main analysis":
        src = T_RES
    elif lab.startswith("Medium"):
        src = T_STRATA
    elif lab.startswith("Platform"):
        src = T_PLAT
    else:
        src = T_SENS
    nm.add("d", lab, f"{v * 100:.1f}% (n={n})", src, note)
nm.add("d", "medium_source", "confirm_strata_curated_v2", T_STRATA, "FBS/hPL/xeno/unknown four rows")
nm.add("d", "platform_source", "fig5d_platform_accuracy_v2", T_PLAT, "microarray/RNA-seq two rows")

# =====================================================================
save_fig(fig, "Fig5_validation")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "Fig5_validation.png"),
           os.path.join(FIG, "qc", "Fig5_validation_deutan.png"))
print("done")

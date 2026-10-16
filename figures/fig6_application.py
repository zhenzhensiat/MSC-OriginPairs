"""fig6_application.py —— Figure 6 (application scope).

Panels: a non-baseline categories: accuracy + low-confidence proportion;
b other-tissue call composition with confidence distributions and threshold;
c simulated admixture: AUC and detection rate vs contamination fraction;
d passage: relative passage x confidence, one line per donor.

Data: arm_robustness_summary_v2.tsv, openset_summary_v2.tsv /
openset_pred_v2.tsv / openset_auc_v2.tsv, admixture_roc_by_combo_v2.tsv +
admixture_roc_v2.tsv, passage_accuracy_pred_v2.tsv, confirm_predictions_v2.tsv.
Colours: source colours; low confidence #C9C9C9; misclassified points red
#D62728 (semantics: error). Python performs counting/proportions only.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (SRC_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

FIG_NAME = "Fig6"
nm = Numbers(FIG_NAME)

T_ARM = "data/arm_robustness_summary_v2.tsv"
T_OS = "data/openset_summary_v2.tsv"
T_OSP = "data/openset_pred_v2.tsv"
T_OSA = "data/openset_auc_v2.tsv"
T_ADB = "data/admixture_roc_by_combo_v2.tsv"
T_AD = "data/admixture_roc_v2.tsv"
T_PG = "data/passage_accuracy_pred_v2.tsv"
T_CP = "data/confirm_predictions_v2.tsv"

arm = read_table(T_ARM)
osu = read_table(T_OS)
osp = read_table(T_OSP)
osa = read_table(T_OSA)
adb = read_table(T_ADB)
ad = read_table(T_AD)
pg = read_table(T_PG)
cp = read_table(T_CP)

C_AD, C_BM, C_UC = SRC_COLORS["AD"], SRC_COLORS["BM"], SRC_COLORS["UC"]
C_LC = "#C9C9C9"      # low confidence (neutral light gray)
C_WRONG = "#D62728"
THR = 0.2494

# =====================================================================
# 0. Hard checks
# =====================================================================
assert len(arm) == 5 and set(arm["arm_type"]) == {"disease_donor", "treated", "primed",
                                                  "culture_variant", "other"}
aa = arm.set_index("arm_type")
assert int(aa.loc["disease_donor", "n"]) == 85 and abs(float(aa.loc["disease_donor", "accuracy"]) - 1) < 1e-12
assert int(aa.loc["treated", "n"]) == 27 and int(aa.loc["culture_variant", "n"]) == 22
assert abs(float(aa.loc["culture_variant", "accuracy"]) - 0.954545454545455) < 1e-9
assert abs(float(aa.loc["culture_variant", "rejected"]) - 0.409090909090909) < 1e-9
assert abs(float(aa.loc["disease_donor", "rejected"]) - 0.0588235294117647) < 1e-9

assert len(osu) == 8 and int(osu["n"].sum()) == 65
assert int(osp.shape[0]) == 65
auc_overall = float(osa.loc[osa["comparison"] == "C_inclass_vs_all_ooc", "auc"].iloc[0])
assert abs(auc_overall - 0.748065543923532) < 1e-9

assert len(adb) == 210 and adb.groupby(["major", "minor"]).ngroups == 6
for _, r in ad.iterrows():
    sub = adb[adb["minor_frac"] == r["minor_frac"]]
    assert len(sub) == 42
    assert abs(sub["auc"].median() - float(r["auc"])) < 1e-9
assert abs(float(ad.loc[ad["minor_frac"] == 0.1, "sens_at_spec95"].iloc[0]) - 0.475) < 1e-9

assert len(pg) == 200 and pg["study"].nunique() == 10
wrong = ~pg["correct"].astype(str).str.upper().eq("TRUE")
assert int(wrong.sum()) == 6
assert int(wrong[pg["study"] == "GSE46019"].sum()) == 4
g73 = pg[pg["study"] == "GSE139073"]
assert len(g73) == 33 and int(g73["rejected"].astype(str).str.upper().eq("TRUE").sum()) == 33

m_in = pd.to_numeric(cp["margin"], errors="coerce")
assert int(m_in.notna().sum()) == 169 and int((m_in < THR).sum()) == 18
assert int((osp["margin"] < THR).sum()) == 32

# =====================================================================
# Layout
# =====================================================================
fig = plt.figure(figsize=(180 * MM, 190 * MM))
outer = fig.add_gridspec(3, 1, height_ratios=[72, 26, 82], hspace=0.0,
                         left=0.082, right=0.992, top=0.97, bottom=0.03)
row1 = outer[0].subgridspec(1, 3, width_ratios=[1.02, 0.98, 0.64], wspace=0.40)
row2 = outer[2].subgridspec(1, 2, width_ratios=[0.92, 1.3], wspace=0.30)

# =====================================================================
# 6a Non-baseline five classes: accuracy + low confidence
# =====================================================================
axa = fig.add_subplot(row1[0])
ARMS = [("disease_donor", "Disease donors"), ("treated", "Treated"), ("primed", "Primed"),
        ("culture_variant", "Culture variants"), ("other", "Fetal BM")]
for k, (key, lab) in enumerate(ARMS):
    y0 = 4 - k
    r = aa.loc[key]
    axa.add_patch(Rectangle((0, y0 + 0.02), 100, 0.30, fc="#F0F0F0", ec="none", zorder=1))
    axa.add_patch(Rectangle((0, y0 + 0.02), float(r["accuracy"]) * 100, 0.30,
                            fc="#333333", ec="none", zorder=2))
    axa.text(float(r["accuracy"]) * 100 + 1.2, y0 + 0.17, f"{float(r['accuracy']) * 100:.1f}",
             fontsize=6, va="center", color=DGREY)
    axa.add_patch(Rectangle((0, y0 - 0.36), 100, 0.30, fc="#F0F0F0", ec="none", zorder=1))
    if float(r["rejected"]) > 0:
        axa.add_patch(Rectangle((0, y0 - 0.36), float(r["rejected"]) * 100, 0.30,
                                fc=C_LC, ec="none", zorder=2))
        axa.text(float(r["rejected"]) * 100 + 1.2, y0 - 0.21, f"{float(r['rejected']) * 100:.1f}",
                 fontsize=6, va="center", color=DGREY)
    axa.text(131, y0, f"n = {int(r['n'])}", fontsize=6, va="center", ha="right",
             color=DGREY)
axa.set_yticks([4, 3, 2, 1, 0])
axa.set_yticklabels([l for _, l in ARMS], fontsize=6)
axa.tick_params(axis="y", length=0, pad=2)
axa.set_ylim(-0.7, 4.85)
axa.set_xlim(0, 132)
axa.set_xticks([0, 25, 50, 75, 100])
axa.set_xlabel("Samples (%)", fontsize=6.5, labelpad=1)
axa.legend(handles=[Patch(fc="#333333", label="Correct"),
                    Patch(fc=C_LC, label="Low-confidence")],
           loc="upper left", bbox_to_anchor=(0.0, 1.03), fontsize=6, ncol=2,
           handlelength=1.0, columnspacing=1.0, handletextpad=0.4, borderpad=0.2)
panel_label(axa, "A", x=-0.30, y=1.05)
for key, lab in ARMS:
    r = aa.loc[key]
    nm.add("a", f"acc_{key}", f"{float(r['accuracy']) * 100:.1f}% (n={int(r['n'])})", T_ARM, lab)
    nm.add("a", f"lowconf_{key}", f"{float(r['rejected']) * 100:.1f}%", T_ARM, lab)
nm.add("a", "note_other", "arm_type=other (3 samples, arm detail = fetal bovine marrow, fetal BM)", T_ARM,
       "v0.2(3) caption lists only four classes; fifth needs adding")

# =====================================================================
# 6b left: call composition for 8 other tissues (sorted by low-confidence fraction)
# =====================================================================
axb = fig.add_subplot(row1[1])
ordb = osu.sort_values(["rejection_rate", "n"], ascending=[False, False])
xpos = np.arange(len(ordb))
for i, (_, r) in enumerate(ordb.iterrows()):
    n = int(r["n"])
    segs = [("called_UC", C_UC), ("called_AD", C_AD), ("called_BM", C_BM),
            ("rejection_rate", C_LC)]
    bottom = 0.0
    for col, colr in segs:
        v = float(r[col]) * 100
        if v <= 0:
            continue
        axb.add_patch(Rectangle((i - 0.32, bottom), 0.64, v, fc=colr, ec="white", lw=0.3,
                                zorder=2))
        cnt = v * n / 100.0
        if v >= 8:
            axb.text(i, bottom + v / 2, f"{cnt:.0f}", ha="center", va="center", fontsize=6,
                     color=("white" if colr == C_BM else "black"))
        bottom += v
    assert abs(bottom - 100) < 0.01
axb.set_xticks(xpos)
_labs = []
for _, r in ordb.iterrows():
    t = r["tissue"].replace("_", " ")
    t = "iPSC-derived" if t == "ipsc derived" else t.capitalize()
    _labs.append(f"{t} ({int(r['n'])})")
axb.set_xticklabels(_labs, rotation=90, fontsize=6)
axb.set_xlim(-0.6, len(ordb) - 0.4)
axb.set_ylim(0, 100)
axb.set_yticks([0, 25, 50, 75, 100])
axb.set_ylabel("Samples (%)", fontsize=6.5)
axb.legend(handles=[Patch(fc=C_AD, label="AD"),
                    Patch(fc=C_BM, label="BM"), Patch(fc=C_LC, label="Low-confidence")],
           loc="lower center", bbox_to_anchor=(0.5, 1.005), fontsize=6, ncol=3,
           handlelength=1.0, columnspacing=0.9, handletextpad=0.35, borderpad=0.15)
panel_label(axb, "B", x=-0.22, y=1.10)
nm.add("b", "order", " → ".join(ordb["tissue"]), T_OS, "descending by low-confidence fraction")
for _, r in ordb.iterrows():
    nm.add("b", f"comp_{r['tissue']}",
           f"rej {float(r['rejection_rate']) * 100:.0f}% / AD {float(r['called_AD']) * 100:.0f}% / "
           f"BM {float(r['called_BM']) * 100:.0f}% / UC {float(r['called_UC']) * 100:.0f}% (n={int(r['n'])})",
           T_OS, "call composition")

# =====================================================================
# 6b right: confidence distribution (in-class vs other tissues) + threshold
# =====================================================================
axbm = fig.add_subplot(row1[2])
rng = np.random.default_rng(11)
groups = [("BM / AD / UC\n(n = 169)", m_in.dropna().to_numpy(), "#333333"),
          ("Other tissues\n(n = 65)", osp["margin"].to_numpy(), "#9E9E9E")]
for yi, (lab, vals, colr) in enumerate(groups):
    y0 = 1 - yi
    axbm.add_patch(Rectangle((0, y0 - 0.28), 1.0, 0.56, fc="#F5F5F5", ec="none", zorder=1))
    jit = rng.uniform(-0.17, 0.17, len(vals))
    axbm.scatter(vals, y0 + jit, s=4.5, c=colr, linewidths=0, alpha=0.75, zorder=3)
axbm.plot([THR, THR], [-0.6, 1.42], color="black", lw=0.9, ls="--", zorder=4)
axbm.text(THR - 0.012, 1.30, f"{THR}", fontsize=6, va="center", ha="right")
axbm.text(0.02, 1.55, f"AUC {auc_overall:.2f}", fontsize=6, va="center", ha="left",
          color=DGREY)
axbm.set_yticks([1, 0])
axbm.set_yticklabels([g[0] for g in groups], fontsize=6)
axbm.tick_params(axis="y", length=0, pad=2)
axbm.set_ylim(-0.6, 1.72)
axbm.set_xlim(0, 1)
axbm.set_xticks([0, 0.5, 1])
axbm.set_xlabel("Confidence (margin)", fontsize=6.5, labelpad=1)
nm.add("bm", "auc_overall", f"{auc_overall:.4f}", T_OSA, "C_inclass_vs_all_ooc (confidence discrimination)")
nm.add("bm", "n_inclass", 169, T_CP, "evaluable external validation samples (margin not NA)")
nm.add("bm", "n_ooc", 65, T_OSP, "other-tissue samples")
nm.add("bm", "lowconf_inclass", 18, T_CP, "margin<0.2494")
nm.add("bm", "lowconf_ooc", 32, T_OSP, "margin<0.2494")

# =====================================================================
# 6c Mixed samples: AUC and detection rate vs mixing fraction
# =====================================================================
axc = fig.add_subplot(row2[0])
fr = sorted(ad["minor_frac"].unique())
x = np.arange(len(fr))
for metric, colr, style, lab, ly in [("auc", "black", "-", "AUC", 1.075),
                                     ("sens_at_spec95", "#666666", "--", "Detection\n(5% FPR)", 0.93)]:
    meds, los, his = [], [], []
    for f in fr:
        sub = adb.loc[adb["minor_frac"] == f, metric]
        meds.append(np.median(sub))
        los.append(np.percentile(sub, 25))
        his.append(np.percentile(sub, 75))
    axc.fill_between(x, los, his, color=colr, alpha=0.15, lw=0)
    axc.plot(x, meds, color=colr, ls=style, lw=1.2, marker="o", ms=3)
    axc.text(x[-1] + 0.14, ly, lab, fontsize=6, va="center", ha="left",
             color=colr, linespacing=1.2)
axc.set_xticks(x)
axc.set_xticklabels([f"{int(f * 100)}%" for f in fr])
axc.set_xlim(-0.35, 5.55)
axc.set_ylim(0, 1.14)
axc.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axc.set_xlabel("Minor component (%)", fontsize=6.5, labelpad=1)
axc.set_ylabel("AUC or detection rate", fontsize=6.5)
panel_label(axc, "C", x=-0.20, y=1.04)
for f in fr:
    sub = adb.loc[adb["minor_frac"] == f]
    nm.add("c", f"auc_{int(f * 100)}", f"{np.median(sub['auc']):.3f} "
           f"(IQR {np.percentile(sub['auc'], 25):.3f}-{np.percentile(sub['auc'], 75):.3f})",
           T_ADB, "42 combinations; 0–1 scale")
    nm.add("c", f"det_{int(f * 100)}", f"{np.median(sub['sens_at_spec95']):.3f} "
           f"(IQR {np.percentile(sub['sens_at_spec95'], 25):.3f}-"
           f"{np.percentile(sub['sens_at_spec95'], 75):.3f})", T_ADB,
           "5% false positives on pure samples; 0–1 scale")
nm.add("c", "n_combos", 42, T_ADB, "number of study×major×minor component combinations")

# =====================================================================
# 6d Passaging: relative passage × confidence (one line per donor; point color=right/wrong)
# =====================================================================
axd = fig.add_subplot(row2[1])
STAGE = {"early": 0, "mid": 1, "late": 2}
pg2 = pg.copy()
pg2["xc"] = pg2["band_research"].map(STAGE)
assert pg2["xc"].notna().all()
for (st, dn), d in pg2.groupby(["study", "donor"]):
    med = d.groupby("xc")["margin"].median().sort_index()
    if len(med) >= 2:
        axd.plot(med.index, med.values, color="#999999", lw=0.55, alpha=0.6, zorder=2)
ok = pg2["correct"].astype(str).str.upper().eq("TRUE")
axd.scatter(pg2.loc[ok, "xc"], pg2.loc[ok, "margin"], s=4.5, c="#333333", linewidths=0,
            alpha=0.7, zorder=3)
axd.scatter(pg2.loc[~ok, "xc"], pg2.loc[~ok, "margin"], s=11, c="black", marker="x",
            linewidths=0.8, zorder=4)
axd.axhline(THR, color="black", lw=0.9, ls="--", zorder=1)
axd.text(2.42, THR, f"{THR}", fontsize=6, va="center", ha="left",
         bbox=dict(fc="white", ec="none", pad=0.6))
axd.text(0.42, 0.085, "GSE139073\n(100% low-confidence)", fontsize=6, color=DGREY,
         linespacing=1.25, bbox=dict(fc="white", ec="none", pad=0.6))
axd.text(2.06, 0.33, "GSE46019\n(4 errors)", fontsize=6, color="black", linespacing=1.25,
         bbox=dict(fc="white", ec="none", pad=0.6))
axd.set_xticks([0, 1, 2])
axd.set_xticklabels(["Early", "Mid", "Late"])
axd.set_xlim(-0.45, 2.6)
axd.set_ylim(-0.03, 1.0)
axd.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axd.set_xlabel("Within-study passage stage", fontsize=6.5, labelpad=1)
axd.set_ylabel("Confidence (margin)", fontsize=6.5)
axd.legend(handles=[Patch(fc="#333333", label="Correct call"),
                    plt.matplotlib.lines.Line2D([], [], color="black", marker="x",
                                                lw=0, ms=3.5, label="Wrong call"),
                    plt.matplotlib.lines.Line2D([], [], color="#999999", lw=0.7,
                                                label="Per-donor line")],
           loc="lower left", bbox_to_anchor=(0.0, 1.015), fontsize=6, ncol=3,
           handlelength=1.2, columnspacing=1.3, handletextpad=0.4, borderpad=0.2)
panel_label(axd, "D", x=-0.10, y=1.04)
nm.add("d", "n_samples", 200, T_PG, "10 passaging studies")
nm.add("d", "n_donors", int(pg2.groupby(['study', 'donor']).ngroups), T_PG, "number of donors (study×donor)")
nm.add("d", "wrong_total", 6, T_PG, "total misjudged samples")
nm.add("d", "wrong_gse46019", 4, T_PG, "GSE46019 (2 early / 2 late)")
nm.add("d", "gse139073", "33/33 low-confidence", T_PG, "remaining 2 mid-passage errors see main text (GSE309666/GSE34303)")
nm.add("d", "threshold", THR, "R/utils_models.R choose_tau", "dashed line")

# =====================================================================
save_fig(fig, "Fig6_application")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "Fig6_application.png"),
           os.path.join(FIG, "qc", "Fig6_application_deutan.png"))
print("done")

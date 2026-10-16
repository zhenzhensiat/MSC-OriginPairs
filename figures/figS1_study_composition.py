"""figS1_study_composition.py —— Supplementary Figure S1: per-study sample and
platform composition for each cohort, plus the complete research-group network
(76 datasets; same assertions as Figure 1).

Data: data/platform_composition_v2.tsv,
data/independence_network_v2_r2.tsv.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (SRC_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

nm = Numbers("FigS1")
T_PC = "data/platform_composition_v2.tsv"
T_NET = "data/independence_network_v2_r2.tsv"
pc = read_table(T_PC)
net = read_table(T_NET)
la = read_table("data/layer_allocation_v2.tsv")
cb = read_table("data/c_layer_build_report_v2.tsv")
pg = read_table("data/passage_accuracy_pred_v2.tsv")

C_AD, C_BM, C_UC = SRC_COLORS["AD"], SRC_COLORS["BM"], SRC_COLORS["UC"]
C_ARR, C_RNAS = "#4D4D4D", "#B3B3B3"
C_DIS, C_EXT, C_PAS, C_OTHN = "#404040", "#9E9E9E", "#D9D9D9", "#FFFFFF"

D = set(la.loc[la["layer"] == "D", "gse"])
C = set(cb.loc[cb["status"] == "ok", "gse"])
P = set(pg["study"].unique())
assert len(D) == 19 and len(C) == 22 and len(P) == 10


def plat_of(s):
    s = str(s)
    if "microarray" in s:
        return "array"
    if "RNA-seq" in s:
        return "rna"
    return "na"


psg = pc.loc[pc["passage_series"].astype(str).str.upper().eq("TRUE")]
cols = [("Discovery cohort (n = 19)", pc[pc["layer"] == "D"]),
        ("External validation (n = 22)", pc[pc["layer"] == "C"]),
        ("Passage-series (n = 10)", psg)]
for lab, sub in cols:
    assert int(sub["n_other"].sum()) == 0

fig = plt.figure(figsize=(180 * MM, 250 * MM))
outer = fig.add_gridspec(5, 1, height_ratios=[100, 8, 86, 8, 56], hspace=0.0,
                         left=0.02, right=0.985, top=0.975, bottom=0.025)

# ---------------- a: per-study composition ----------------
axa = fig.add_subplot(outer[0])
axa.set_xlim(0, 180)
axa.set_ylim(0, 104)
axa.axis("off")
X0 = [2.0, 63.0, 124.0]
for (title, sub), x0 in zip(cols, X0):
    sub = sub.copy()
    sub["k"] = sub["gse"].str.replace("GSE", "", regex=False).astype(int)
    sub = sub.sort_values("k")
    axa.text(x0, 101.5, title, fontsize=6.5, ha="left", va="center")
    axa.text(x0, 97.0, "GSE", fontsize=6, ha="left", va="center", color=DGREY)
    axa.text(x0 + 14, 97.0, "Samples by source", fontsize=6, ha="left", va="center",
             color=DGREY)
    axa.text(x0 + 45, 97.0, "Plat.", fontsize=6, ha="left", va="center", color=DGREY)
    y = 92.5
    for _, r in sub.iterrows():
        axa.text(x0, y, r["gse"], fontsize=6, ha="left", va="center")
        tot = int(r["n_samples"])
        xb = x0 + 14.0
        for key, colr, cnt in [("n_AD", C_AD, int(r["n_AD"])), ("n_BM", C_BM, int(r["n_BM"])),
                               ("n_UC", C_UC, int(r["n_UC"]))]:
            if cnt <= 0:
                continue
            w = 24.0 * cnt / tot
            axa.add_patch(Rectangle((xb, y - 1.5), w, 3.0, fc=colr, ec="white", lw=0.3))
            if cnt / tot >= 0.12:
                axa.text(xb + w / 2, y, f"{cnt}", ha="center", va="center", fontsize=6,
                         color=("white" if colr == C_BM else "black"))
            xb += w
        p = plat_of(r["platform_plot"])
        axa.add_patch(Rectangle((x0 + 45.5, y - 1.3), 2.6, 2.6,
                                fc=(C_ARR if p == "array" else C_RNAS if p == "rna" else "#F0F0F0"),
                                ec="#DDDDDD", lw=0.3))
        axa.text(x0 + 49.5, y, f"{tot}", fontsize=6, ha="left", va="center", color=DGREY)
        y -= 4.28
        nm.add("a", f"{r['gse']}", f"AD {int(r['n_AD'])} / BM {int(r['n_BM'])} / "
               f"UC {int(r['n_UC'])}; n={tot}; {p}", T_PC, title)
axa.legend(handles=[Patch(fc=C_AD, label="AD"), Patch(fc=C_BM, label="BM"),
                    Patch(fc=C_UC, label="UC"), Patch(fc=C_ARR, label="Microarray"),
                    Patch(fc=C_RNAS, label="RNA-seq")],
           loc="lower right", bbox_to_anchor=(0.999, -0.01), fontsize=6, ncol=5,
           handlelength=1.0, columnspacing=1.0, handletextpad=0.35, borderpad=0.2)
panel_label(axa, "A", x=0.002, y=1.0)

# ---------------- b: research-group network (full version, 76 nodes) ----------------
axb = fig.add_subplot(outer[2])
OO = set(pc.loc[pc["layer"] == "O_ooc", "gse"])
OT = set(net["gse"]) - D - C - P
assert len(OO) == 9 and len(OT) == 16
FULL = D | C | P | OO | OT
assert len(FULL) == 76
rank = {"D": 0, "C": 1, "P": 2, "O": 3, "OOC": 4}
cat_of = {}
for g in FULL:
    cat_of[g] = "D" if g in D else "C" if g in C else "P" if g in P else "OOC" if g in OO else "O"
order = sorted(FULL, key=lambda g: (rank[cat_of[g]], int(g.replace("GSE", ""))))
num_of = {g: i + 1 for i, g in enumerate(order)}
assert [num_of[g] for g in order[:19]] == list(range(1, 20))
assert [num_of[g] for g in order[19:41]] == list(range(20, 42))

team = read_table("curation/team_numbering.tsv")
tid2 = dict(zip(team["gse"], team["team_id"]))
g02grp = {"GSE12274", "GSE189895", "GSE194298", "GSE37065", "GSE57151", "GSE9520"}
ucpair2 = {"GSE72597", "GSE76295"}
g44grp = {"GSE291956", "GSE309666"}
assert all(g in D for g in g02grp) and all(g in C for g in ucpair2) and all(g in P for g in g44grp)
assert len({tid2[g] for g in g02grp}) == 1 and len({tid2[g] for g in ucpair2}) == 1 \
    and len({tid2[g] for g in g44grp}) == 1

nt_edges2 = set()
for g, lk in zip(net["gse"], net["linked_gse"].fillna("")):
    for l in str(lk).split(";"):
        if l:
            nt_edges2.add(tuple(sorted([g, l])))
nt_edges2.add(("GSE72597", "GSE76295"))          # A1.2 fix: same-cluster edge added in data layer
assert len(nt_edges2) == 17                       # 15 series + 1 uc + 1 same-series


def _hex6(cx, cy, rx, ry):
    ang = np.deg2rad([90, 30, -30, -90, -150, 150])
    return [(cx + rx * np.cos(a), cy + ry * np.sin(a)) for a in ang]


def _grid(gs, xs, ys):
    pos = {}
    for i, g in enumerate(gs):
        r, c = divmod(i, len(xs))
        pos[g] = (xs[c], ys[r])
    return pos


pos2 = {}
knum = lambda g: int(g.replace("GSE", ""))  # noqa: E731
pos2.update(dict(zip(sorted(g02grp), _hex6(0.070, 0.735, 0.042, 0.115))))
pos2.update(_grid(sorted(D - g02grp, key=knum), [0.175, 0.222, 0.269, 0.316],
                  [0.90, 0.775, 0.65, 0.49]))                # last row shifted down: separate nodes 19 and 13
pos2.update({sorted(ucpair2, key=knum)[0]: (0.395, 0.86),
             sorted(ucpair2, key=knum)[1]: (0.442, 0.86)})
pos2.update(_grid(sorted(C - ucpair2, key=knum), [0.50, 0.612, 0.724, 0.836, 0.948],
                  [0.90, 0.775, 0.65, 0.525]))
pos2.update({sorted(g44grp, key=knum)[0]: (0.030, 0.30), sorted(g44grp, key=knum)[1]: (0.078, 0.30)})
pos2.update(_grid(sorted(P - g44grp, key=knum), [0.140, 0.190, 0.240, 0.290], [0.345, 0.185]))
pos2.update(_grid(sorted(OT, key=knum), [0.360, 0.418, 0.476, 0.534, 0.592, 0.650],
                  [0.345, 0.205, 0.065]))
pos2.update(_grid(sorted(OO, key=knum), [0.720, 0.783, 0.846, 0.909, 0.962], [0.32, 0.12]))

for cx0, cy0, cx1, cy1, lab in [
        (0.008, 0.565, 0.132, 0.925, f"Team {tid2['GSE12274']} · 6"),
        (0.365, 0.77, 0.475, 0.925, f"Team {tid2['GSE72597']} · 2"),
        (0.006, 0.185, 0.104, 0.395, f"Team {tid2['GSE291956']} · 2")]:
    axb.add_patch(FancyBboxPatch((cx0, cy0), cx1 - cx0, cy1 - cy0,
                                 boxstyle="round,pad=0.005,rounding_size=0.012",
                                 fc="#F2F2F2", ec="#BBBBBB", lw=0.6, zorder=0))
    axb.text(cx0 + 0.003, cy1 + 0.007, lab, fontsize=6, color=DGREY, ha="left",
             va="bottom")
for a, b in sorted(nt_edges2):
    (x1, y1), (x2, y2) = pos2[a], pos2[b]
    axb.plot([x1, x2], [y1, y2], color="#AAAAAA", lw=0.6, zorder=1)
for cat, colr, ec, ls in [("D", C_DIS, "none", "-"), ("C", C_EXT, "none", "-"),
                          ("P", C_PAS, "none", "-"), ("O", C_OTHN, "#999999", "-"),
                          ("OOC", "#FFFFFF", "#999999", (0, (2, 1.4)))]:
    gs = [g for g in order if cat_of[g] == cat]
    axb.scatter([pos2[g][0] for g in gs], [pos2[g][1] for g in gs], s=26, c=colr,
                edgecolors=ec, linewidths=0.7, zorder=3,
                linestyles=("solid" if ls == "-" else ls))
for g in order:
    x, y = pos2[g]
    axb.text(x, y + 0.030, str(num_of[g]), fontsize=6, ha="center", va="bottom",
             zorder=4)
for txt, x, y in [("Discovery cohort (19)", 0.008, 0.978),
                  ("External validation cohort (22)", 0.365, 0.978),
                  ("Passage-series studies (10)", 0.006, 0.458),
                  ("Reviewed, not included (16)", 0.360, 0.458),
                  ("Other-tissue datasets, application scope (9)", 0.720, 0.458)]:
    axb.text(x, y, txt, fontsize=6, color=DGREY, ha="left", va="center")
axb.text(0.5, 1.045, "Three shared-author groups in total (6 + 2 + 2); no links between cohorts",
         ha="center", va="bottom", fontsize=6, color=DGREY)
axb.set_xlim(0, 1)
axb.set_ylim(0, 1.1)
axb.axis("off")
axb.legend(handles=[Patch(fc=C_DIS, label="Discovery cohort"),
                    Patch(fc=C_EXT, label="External validation cohort"),
                    Patch(fc=C_PAS, label="Passage-series studies"),
                    Patch(fc=C_OTHN, ec="#999999", label="Reviewed, not included"),
                    Patch(fc="#FFFFFF", ec="#999999", ls=(0, (2, 1.4)),
                          label="Other tissues (application scope)")],
           loc="lower center", bbox_to_anchor=(0.5, -0.035), fontsize=6, ncol=5,
           handlelength=1.0, columnspacing=1.1, handletextpad=0.4, borderpad=0.2)
panel_label(axb, "B", x=-0.005, y=1.02)
nm.add("b", "n_nodes", 76, T_NET, "Full version = D19+C22+P10+Reviewed not included 16+Other-tissue datasets, application scope 9")
nm.add("b", "links", "3 groups (6+2+2); 0 cross-cohort", T_NET, "")
nm.add("b", "fix_edge", "GSE72597–GSE76295", T_NET, "same-cluster edge added in data layer")
nm.add("b", "fix_node", "GSE222263 (node number 25)", T_NET, "restored")
for g in order:
    nm.add("b", f"node_{num_of[g]}", g, T_NET,
           {"D": "Discovery", "C": "External validation", "P": "Passage",
            "O": "Reviewed, not included", "OOC": "Application-scope only"}[cat_of[g]])

# ---------------- c: node number ↔ GEO mapping ----------------
axc = fig.add_subplot(outer[4])
axc.set_xlim(0, 180)
axc.set_ylim(0, 56)
axc.axis("off")
axc.text(2, 53, "Numbering key (Fig 1B and Fig S1B): node number → GEO accession",
         fontsize=6.5, ha="left", va="center")
ncol = 4
per = int(np.ceil(76 / ncol))
for i, g in enumerate(order):
    c, r = divmod(i, per)
    axc.text(6 + c * 44, 48 - r * 3.05, f"{num_of[g]:>2}  {g}", fontsize=6, ha="left",
             va="center", family="monospace")
panel_label(axc, "C", x=0.002, y=0.99)

save_fig(fig, "FigS1_composition")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "FigS1_composition.png"),
           os.path.join(FIG, "qc", "FigS1_composition_deutan.png"))
print("done")

"""fig1_design.py —— Figure 1 (study design).

Layout:
- a  screening funnel: "1,089 screened" at top; "60 datasets contributed samples"
     (51 main + 9 application-scope only); three parallel branches (discovery /
     external validation / application scope); exclusion list (tissue outside
     scope 32, cell line or iPSC-derived 209, etc.); the 47 eligible-but-unused
     datasets annotated in place.
- b  research-group network (41 nodes = D19 + C22): studies from the same group
     enclosed in shaded boxes (one 6-study group, one 2-study pair); no edges
     between cohorts; group labels use the published neutral codes (Team G##).
- c  cohort x source x platform composition.

Data: census_adjudication_v2_r3, layer_allocation_v2,
c_layer_build_report_v2, passage_accuracy_pred_v2, platform_composition_v2,
independence_network_v2_r2, dataset_uses_v2, team numbering table.
Python performs counting only.
"""
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_style import (SRC_COLORS, DGREY, MM, Numbers,  # noqa: E402
                       panel_label, read_table, save_fig, deutan_sim, FIG)

FIG_NAME = "Fig1"
nm = Numbers(FIG_NAME)

T_R3 = "data/census_adjudication_v2_r3.tsv"
T_LA = "data/layer_allocation_v2.tsv"
T_CB = "data/c_layer_build_report_v2.tsv"
T_PG = "data/passage_accuracy_pred_v2.tsv"
T_PC = "data/platform_composition_v2.tsv"
T_NET = "data/independence_network_v2_r2.tsv"
T_USE = "data/dataset_uses_v2.tsv"
T_TEAM = "curation/team_numbering.tsv"

r3 = read_table(T_R3)
la = read_table(T_LA)
cb = read_table(T_CB)
pg = read_table(T_PG)
pc = read_table(T_PC)
net = read_table(T_NET)
use = read_table(T_USE)
team = read_table(T_TEAM)

C_AD, C_BM, C_UC = SRC_COLORS["AD"], SRC_COLORS["BM"], SRC_COLORS["UC"]
C_OTH = "#C9C9C9"
C_DIS, C_EXT = "#404040", "#9E9E9E"
C_ARR, C_RNAS = "#4D4D4D", "#B3B3B3"
GDSTYPE_ARR = {"Expression profiling by array"}
GDSTYPE_SEQ = {"Expression profiling by high throughput sequencing"}

# =====================================================================
# 0. Hard checks
# =====================================================================
assert len(r3) == 1089
D = set(la.loc[la["layer"] == "D", "gse"])
C = set(cb.loc[cb["status"] == "ok", "gse"])
P = set(pg["study"].unique())
OO = set(pc.loc[pc["layer"] == "O_ooc", "gse"])
OA = set(pc.loc[pc["layer"] == "O_arm", "gse"])
assert len(D) == 19 and len(C) == 22 and len(P) == 10 and len(OO) == 9 and len(OA) == 15
assert use.shape[0] == 60
assert set(use["gse"]) == D | C | P | OO
assert (use["uses_en"].str.contains("Discovery").sum(),
        use["uses_en"].str.contains("External").sum(),
        use["uses_en"].str.contains("Passage").sum(),
        use["uses_en"].str.contains("application scope only").sum()) == (19, 22, 10, 9)
assert int(use["n_samples_baseline"].sum()) == 356
assert int(use["n_samples_application"].sum()) == 405

G = {"D": D, "C": C, "P": P, "OOC": OO}
OT = set(net["gse"]) - D - C - P
assert len(OT) == 16, len(OT)
FULL = list(D | C | P | OO | OT)
assert len(FULL) == 76
rank = {"D": 0, "C": 1, "P": 2, "O": 3, "OOC": 4}
cat_of = {}
for g in FULL:
    cat_of[g] = "D" if g in D else "C" if g in C else "P" if g in P else "OOC" if g in OO else "O"
order_all = sorted(FULL, key=lambda g: (rank[cat_of[g]], int(g.replace("GSE", ""))))
num_all = {g: i + 1 for i, g in enumerate(order_all)}
assert [num_all[g] for g in order_all[:19]] == list(range(1, 20))
assert [num_all[g] for g in order_all[19:41]] == list(range(20, 42))

# Exclusion list (after splitting out 9 sets): 982 + 47 + 60 = 1089
cnt = Counter(r3["human_verdict"])
excl = [("Treated, no untreated control", 233),
        ("Cell line or iPSC-derived", 209),
        ("Not MSC", 198),
        ("No expression profiling", 178),
        ("Insufficient metadata", 49),
        ("Non-human", 45),
        ("Tissue outside BM/AD/UC", 32),
        ("Single-cell dataset (exploratory)", 18),
        ("Fetal or embryonic tissue", 5),
        ("Duplicate or super-series", 5),
        ("Too few donors", 4),
        ("Not culture-expanded", 4),
        ("Expression matrix unavailable", 2)]
assert sum(v for _, v in excl) == 982
assert cnt["X_SOURCE_OOC"] - 7 == 32 and cnt["X_CELLLINE_IPSC"] - 2 == 209
assert 982 + 47 + 60 == 1089

# Network: 41-node backbone, two intra-group edges (series group K6=15 edges + uc pair 1 edge)
nt_edges = set()
for g, lk in zip(net["gse"], net["linked_gse"].fillna("")):
    for l in str(lk).split(";"):
        if l:
            nt_edges.add(tuple(sorted([g, l])))
g02grp = {"GSE12274", "GSE189895", "GSE194298", "GSE37065", "GSE57151", "GSE9520"}
ucpair = {"GSE72597", "GSE76295"}
g02_edges = {tuple(sorted(p)) for p in
                [(a, b) for a in g02grp for b in g02grp if a < b]}
assert g02_edges <= nt_edges
assert ("GSE72597", "GSE76295") not in nt_edges            # edge to be added (data-layer evidence)
NET_FIX_EDGE = ("GSE72597", "GSE76295")
main_nodes = D | C
assert "GSE222263" in C                                     # restored node
assert "GSE222263" not in set(net["gse"])                   # this node missing from original network table (presentation-layer fix)
assert all(not str(lk).strip() for g, lk in zip(net["gse"], net["linked_gse"].fillna(""))
           if g == "GSE222263") if "GSE222263" in set(net["gse"]) else True

# Team (G##)
tid = dict(zip(team["gse"], team["team_id"]))
assert len({tid[g] for g in g02grp}) == 1
assert len({tid[g] for g in ucpair}) == 1
assert tid["GSE12274"] != tid["GSE72597"]

# Composition (c): O_ooc platforms filled in from gdstype
gdstype = dict(zip(r3["accession"], r3["gdstype"]))


def plat_from_gdstype(g):
    dt = gdstype.get(g)
    if dt in GDSTYPE_SEQ:
        return "RNA-seq"
    if dt in GDSTYPE_ARR:
        return "microarray"
    return None


assert all(plat_from_gdstype(g) for g in OO), "all O_ooc gdstype entries should be classifiable"


def plat_of(s):
    s = str(s)
    if "microarray" in s:
        return "microarray"
    if "RNA-seq" in s:
        return "RNA-seq"
    return None


# =====================================================================
# Layout: a funnel (full width) → b network (full width) → c composition (full width)
# =====================================================================
fig = plt.figure(figsize=(180 * MM, 218 * MM))
outer = fig.add_gridspec(5, 1, height_ratios=[64, 8, 92, 8, 44], hspace=0.0,
                         left=0.05, right=0.993, top=0.972, bottom=0.022)

# =====================================================================
# 1a Screening funnel (branching version)
# =====================================================================
axa = fig.add_subplot(outer[0])
axa.set_xlim(0, 168)
axa.set_ylim(48, 106)
axa.axis("off")


def _box(ax, x, y, w, h, lines, fs=6, pitch=3.3):
    ax.add_patch(Rectangle((x, y), w, h, fc="white", ec="#999999", lw=0.7, zorder=2))
    n = len(lines)
    for i, (txt, bold) in enumerate(lines):
        yy = y + h / 2 + (n - 1) / 2 * pitch - i * pitch
        ax.text(x + w / 2, yy, txt, ha="center", va="center", fontsize=fs,
                fontweight=("bold" if bold else "normal"), zorder=3)


_box(axa, 4, 99, 86, 6.5, [("1,089 datasets screened", True)])
axa.annotate("", xy=(47, 95.5), xytext=(47, 99),
             arrowprops=dict(arrowstyle="-|>", lw=0.8, color=DGREY))
_box(axa, 4, 86.5, 86, 9,
     [("60 datasets contributed samples", True),
      ("19 discovery · 22 external validation · 10 passage-series · 9 application-scope only",
       False)])
# Three branches (parallel, not serial)
axa.plot([47, 47], [86.5, 81], color=DGREY, lw=0.8)
axa.plot([19, 81], [81, 81], color=DGREY, lw=0.8)
for xc in (19, 50, 81):
    axa.annotate("", xy=(xc, 78), xytext=(xc, 81),
                 arrowprops=dict(arrowstyle="-|>", lw=0.8, color=DGREY))
_box(axa, 2, 58, 34, 20,
     [("Discovery cohort", True),
      ("19 studies / 180 samples", False),
      ("evaluation: 18 studies,", False),
      ("174 samples", False),
      ("final model: 17 studies,", False),
      ("166 samples", False)])
_box(axa, 38, 68, 24, 10,
     [("External validation", True),
      ("22 studies / 176 samples", False),
      ("(no shared authors)", False)])
_box(axa, 64, 56, 34, 22,
     [("Application scope", True),
      ("65 other-tissue samples,", False),
      ("9 datasets (used only here)", False),
      ("140 non-baseline samples,", False),
      ("arms of 15 studies", False),
      ("10 passage-series studies,", False),
      ("200 samples", False)])
# Exclusion list (right side)
axa.plot([97, 99, 99], [101.5, 101.5, 49], color="#BBBBBB", lw=0.7, zorder=1)
axa.text(102, 101.5, "Not included (n = 982)", ha="left", va="center", fontsize=6,
         fontweight="bold")
yy = 96.8
for lab, v in excl:
    axa.text(102, yy, f"{v}  {lab}", ha="left", va="center", fontsize=6)
    yy -= 3.2
axa.text(102, 54.8, "Eligible but not used (n = 47)", ha="left", va="center", fontsize=6,
         fontweight="bold")
axa.text(102, 51.8, "(v1-used 26 · not built 14 · passage-series not required 7)",
         ha="left", va="center", fontsize=6, color=DGREY)
axa.text(102, 48.7, "The 9 application-scope-only datasets are counted under “contributed”",
         ha="left", va="center", fontsize=6, color=DGREY)
panel_label(axa, "A", x=0.002, y=0.99)
for lab, v in excl:
    nm.add("a", f"excl_{lab[:26]}", v, T_R3, f"Excluded: {lab}")
nm.add("a", "not_carried", 47, T_R3,
       "Eligible but not used (v1-used 26 + not constructed 14 + not passaged 7; standalone subtitle, not counted in 982)")
nm.add("a", "contributed", "60 = 19 + 22 + 10 + 9", T_USE, "982 + 47 + 60 = 1089 closed loop")
nm.add("a", "app_scope_note", "9 dedicated sets: 7×X_SOURCE_OOC + 2×X_CELLLINE_IPSC (split out from the exclusion side)",
       T_USE, "of which 272480/310465 are iPSC-derived MSC")
nm.add("a", "discovery", "19 studies / 180 samples (174 evaluated; training models 17 / 166)", T_PC, "")
nm.add("a", "external", "22 studies / 176 samples", T_PC, "")
nm.add("a", "application", "other-tissue 65 (9 sets) · non-baseline 140 (15-study additional arm) · passaging 200 (10 studies)",
       T_PC, "")

# =====================================================================
# 1b Study-team relationship graph (41-node backbone; envelope boxes; edges routed around nodes; 1 edge and 1 node added)
# =====================================================================
axb = fig.add_subplot(outer[2])
axb.set_xlim(0, 1)
axb.set_ylim(0, 1.06)
axb.axis("off")


def _hex6(cx, cy, rx, ry):
    ang = np.deg2rad([90, 30, -30, -90, -150, 150])
    return [(cx + rx * np.cos(a), cy + ry * np.sin(a)) for a in ang]


def _grid(gs, xs, ys):
    pos = {}
    for i, g in enumerate(gs):
        r, c = divmod(i, len(xs))
        pos[g] = (xs[c], ys[r])
    return pos


# D: series group (hexagon, K6) + 13 singletons
g02_sorted = sorted(g02grp, key=lambda g: int(g.replace("GSE", "")))
pos_b = dict(zip(g02_sorted, _hex6(0.125, 0.62, 0.078, 0.235)))
d_singles = sorted(D - g02grp, key=lambda g: int(g.replace("GSE", "")))
pos_b.update(_grid(d_singles, [0.30, 0.355, 0.41, 0.462], [0.94, 0.76, 0.58, 0.40]))
# C: uc pair + 20 singletons
uc_sorted = sorted(ucpair)
pos_b.update({uc_sorted[0]: (0.605, 0.86), uc_sorted[1]: (0.668, 0.86)})
c_singles = sorted(C - ucpair, key=lambda g: int(g.replace("GSE", "")))
pos_b.update(_grid(c_singles, [0.565, 0.652, 0.739, 0.826, 0.913], [0.64, 0.44, 0.24, 0.05]))

# Envelope boxes (light color)
for cx0, cy0, cx1, cy1, lab in [
        (0.032, 0.27, 0.218, 0.965, f"Team {tid['GSE12274']} · 6 datasets"),
        (0.545, 0.77, 0.728, 0.965, f"Team {tid['GSE72597']} · 2 datasets")]:
    axb.add_patch(FancyBboxPatch((cx0, cy0), cx1 - cx0, cy1 - cy0,
                                 boxstyle="round,pad=0.006,rounding_size=0.015",
                                 fc="#F2F2F2", ec="#BBBBBB", lw=0.6, zorder=0))
    axb.text(cx0 + 0.004, cy1 + 0.008, lab, fontsize=6, color=DGREY, ha="left", va="bottom")

edges_main = [e for e in nt_edges if e[0] in main_nodes and e[1] in main_nodes] + [NET_FIX_EDGE]
for a, b in edges_main:
    (x1, y1), (x2, y2) = pos_b[a], pos_b[b]
    axb.plot([x1, x2], [y1, y2], color="#AAAAAA", lw=0.6, zorder=1, alpha=0.9)
for cat, colr in [("D", C_DIS), ("C", C_EXT)]:
    gs = [g for g in main_nodes if cat_of[g] == cat]
    axb.scatter([pos_b[g][0] for g in gs], [pos_b[g][1] for g in gs], s=26, c=colr,
                zorder=3)
for g in sorted(main_nodes, key=lambda g: num_all[g]):
    x, y = pos_b[g]
    axb.text(x, y + 0.022, str(num_all[g]), fontsize=6, ha="center", va="bottom",
             zorder=4)
axb.text(0.035, 1.048, "Discovery cohort (19 studies)", fontsize=6, color=DGREY)
axb.text(0.545, 1.048, "External validation cohort (22 studies)", fontsize=6, color=DGREY)
axb.text(0.5, 1.008, "Two shared-author groups (6 + 2 datasets); no links between cohorts",
         ha="center", va="bottom", fontsize=6, color=DGREY)
axb.legend(handles=[Patch(fc=C_DIS, label="Discovery cohort"),
                    Patch(fc=C_EXT, label="External validation cohort")],
           loc="lower center", bbox_to_anchor=(0.5, -0.015), fontsize=6, ncol=2,
           handlelength=1.0, columnspacing=1.4, handletextpad=0.4, borderpad=0.2)
panel_label(axb, "B", x=-0.003, y=1.03)
nm.add("b", "n_nodes", 41, T_NET, "backbone network (D 19 + C 22)")
nm.add("b", "n_edges", 16, T_NET, "series group K6 15 edges + uc pair 1 edge")
nm.add("b", "fix_edge", "GSE72597–GSE76295", T_USE,
       "missing from original network table; same data-layer cluster G20")
nm.add("b", "fix_node", "GSE222263", T_USE, "one of 22 external validation items; missing from original network table")
nm.add("b", "groups", f"series group=Team {tid['GSE12274']} (6); uc pair=Team {tid['GSE72597']} (2)", T_TEAM, "envelope box label")
for g in sorted(main_nodes, key=lambda g: num_all[g]):
    nm.add("b", f"node_{num_all[g]}", g, T_NET,
           {"D": "Discovery", "C": "External validation"}[cat_of[g]] +
           f"; Team {tid.get(g, 'NA')}")

# =====================================================================
# 1c Cohort × source × platform (formerly 1d; O_ooc platforms filled in from gdstype)
# =====================================================================
axd = fig.add_subplot(outer[4])
axd.set_xlim(0, 100)
axd.set_ylim(-7, 58)
axd.axis("off")


def _row(y0, lab, nsamp, nstud, src_counts, plat_counts):
    axd.text(1.5, y0 + 1.7, lab, fontsize=6, ha="left", va="center")
    x = 16.0
    total = sum(src_counts.values())
    for key, colr in [("AD", C_AD), ("BM", C_BM), ("UC", C_UC), ("other", C_OTH)]:
        v = src_counts.get(key, 0)
        if v <= 0:
            continue
        w = 34.0 * v / total
        axd.add_patch(Rectangle((x, y0), w, 3.4, fc=colr, ec="white", lw=0.3, zorder=2))
        if v / total >= 0.08:
            axd.text(x + w / 2, y0 + 1.7, f"{v}", ha="center", va="center", fontsize=6,
                     color=("white" if colr == C_BM else "black"))
        x += w
    x = 56.0
    tot_p = sum(plat_counts.values())
    for key, colr, tcol in [("microarray", C_ARR, "white"), ("RNA-seq", C_RNAS, "black")]:
        v = plat_counts.get(key, 0)
        if v <= 0:
            continue
        w = 34.0 * v / tot_p
        axd.add_patch(Rectangle((x, y0), w, 3.4, fc=colr, ec="white", lw=0.3, zorder=2))
        if v / tot_p >= 0.08:
            axd.text(x + w / 2, y0 + 1.7, f"{v}", ha="center", va="center", fontsize=6,
                     color=tcol)
        x += w
    axd.text(99.5, y0 + 1.7, f"{nsamp} · {nstud}", ha="right", va="center", fontsize=6,
             color=DGREY)


psg = pc.loc[pc["passage_series"].astype(str).str.upper().eq("TRUE")]
for i, (key, lab) in enumerate([("D", "Discovery cohort"), ("C", "External validation"),
                                ("P", "Passage-series studies"), ("O_arm", "Non-baseline arms"),
                                ("O_ooc", "Other tissues")]):
    sub = psg if key == "P" else pc[pc["layer"] == key]
    y0 = 44 - i * 11
    src = {"AD": int(sub["n_AD"].sum()), "BM": int(sub["n_BM"].sum()),
           "UC": int(sub["n_UC"].sum()), "other": int(sub["n_other"].sum())}
    pk = {}
    for _, r0 in sub.iterrows():
        k = plat_of(r0["platform_plot"]) or plat_from_gdstype(r0["gse"])
        assert k, r0["gse"]
        pk[k] = pk.get(k, 0) + int(r0["n_samples"])
    _row(y0, lab, int(sub["n_samples"].sum()), int(sub.shape[0]), src, pk)
    nm.add("c", f"src_{lab}",
           f"AD {src['AD']} / BM {src['BM']} / UC {src['UC']} / other {src['other']}",
           T_PC, f"{lab}; {int(sub['n_samples'].sum())} samples total ({int(sub.shape[0])} datasets)")
    nm.add("c", f"plat_{lab}",
           f"microarray {pk.get('microarray', 0)} / RNA-seq {pk.get('RNA-seq', 0)}",
           T_PC if key != "O_ooc" else T_R3,
           "by sample count" + ("; O_ooc platforms filled in from r3.gdstype" if key == "O_ooc" else ""))
axd.text(33, 55, "Source", ha="center", va="center", fontsize=6, color=DGREY)
axd.text(73, 55, "Platform", ha="center", va="center", fontsize=6, color=DGREY)
axd.legend(handles=[Patch(fc=C_AD, label="AD"), Patch(fc=C_BM, label="BM"),
                    Patch(fc=C_UC, label="UC"), Patch(fc=C_OTH, label="Other tissues"),
                    Patch(fc=C_ARR, label="Microarray"), Patch(fc=C_RNAS, label="RNA-seq")],
           loc="lower left", bbox_to_anchor=(-0.002, -0.02), fontsize=6, ncol=6,
           handlelength=1.0, columnspacing=1.1, handletextpad=0.4, borderpad=0.2)
axd.text(99.5, 55, "samples · datasets", ha="right", va="center", fontsize=6, color=DGREY)
panel_label(axd, "C", x=0.002, y=0.99)
nm.add("c", "note_nonbaseline", "140 samples from additional arms of 15 cohort studies", T_PC, "O_arm")
nm.add("c", "note_ooc", "65 samples from 9 out-of-cohort datasets (8 tissues, application-scope only)", T_PC, "O_ooc")

# =====================================================================
save_fig(fig, "Fig1_design")
nm.save()
os.makedirs(os.path.join(FIG, "qc"), exist_ok=True)
deutan_sim(os.path.join(FIG, "Fig1_design.png"), os.path.join(FIG, "qc", "Fig1_design_deutan.png"))
print("done")

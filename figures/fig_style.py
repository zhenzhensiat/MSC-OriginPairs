"""fig_style.py —— shared figure style and result-table reader.

Conventions:
- colours: sources use the Okabe-Ito colour-blind-safe palette (AD #D55E00 /
  BM #0072B2 / UC #009E73); outcome colours: replicated #4B3F72 / reversed
  #E1A23B / not replicated #C9C9C9.
- Python performs display calculations (counts/proportions) only; statistical
  tests are not recomputed here.
- export: figures/ as PDF (fonts embedded) + PNG 300 dpi; double column 180 mm.
- numbers_FigX.tsv records the source of every plotted number.
- deuteranopia simulation (Vienot 1999 approximation) for colour-readability QC.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # msc_v2_code/
FIG = os.path.join(ROOT, "figures")

# ---- colours (uniform across all figures)-------------------------------------------------
C_AD, C_BM, C_UC = "#D55E00", "#0072B2", "#009E73"
SRC_COLORS = {"AD": C_AD, "BM": C_BM, "UC": C_UC}
CL_REP, CL_REV, CL_NOT = "#4B3F72", "#E1A23B", "#C9C9C9"
CL_COLORS = {"replicated": CL_REP, "reversed": CL_REV, "not_replicated": CL_NOT}
CL_LABELS = {"replicated": "Replicated (reported direction)",
             "reversed": "Replicated (opposite direction)",
             "not_replicated": "Not replicated"}
GREY, DGREY, LGREY = "#808080", "#333333", "#C9C9C9"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "figure.dpi": 150, "savefig.dpi": 300,
})
MM = 1.0 / 25.4  # mm → inch


def read_table(rel_path, **kw):
    """Read a result TSV (na_values=NA, matching R write_tsv)."""
    return pd.read_csv(os.path.join(ROOT, rel_path), sep="\t",
                       keep_default_na=False, na_values=["NA"], **kw)


def panel_label(ax, letter, x=-0.14, y=1.05):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left")


class Numbers:
    """Recorder for figures/numbers_FigX.tsv."""
    def __init__(self, fig_name):
        self.fig_name = fig_name
        self.rows = []

    def add(self, panel, label, value, source_table, note=""):
        if isinstance(panel, str) and len(panel) == 1 and panel.isalpha():
            panel = panel.upper()
        self.rows.append({"figure": self.fig_name, "panel": panel, "label": label,
                          "value": str(value), "source_table": source_table, "note": note})

    def save(self):
        os.makedirs(FIG, exist_ok=True)
        path = os.path.join(FIG, f"numbers_{self.fig_name}.tsv")
        pd.DataFrame(self.rows).to_csv(path, sep="\t", index=False)
        print("numbers ->", path)
        return path


def save_fig(fig, name):
    # PDF raster layer 600 dpi; PNG converted to RGB (no alpha) for submission systems
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight", dpi=600)
    png = os.path.join(FIG, name + ".png")
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    from PIL import Image
    Image.open(png).convert("RGB").save(png, dpi=(300, 300))
    print("saved", name + ".pdf/.png")


def deutan_sim(png_in, png_out):
    """deuteranopia simulation (Vienot 1999 linear approximation; colour-readability QC only)."""
    from PIL import Image
    img = np.asarray(Image.open(png_in).convert("RGB")).astype(float) / 255.0

    def lin(c):
        return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)

    def srgb(c):
        return np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1 / 2.4) - 0.055)

    rgb = lin(img)
    M = np.array([[0.625, 0.375, 0.0], [0.70, 0.30, 0.0], [0.0, 0.30, 0.70]])
    out = np.clip(srgb(np.einsum("...c,dc->...d", rgb, M)), 0, 1)
    Image.fromarray((out * 255).round().astype(np.uint8)).save(png_out)
    print("deutan sim ->", png_out)

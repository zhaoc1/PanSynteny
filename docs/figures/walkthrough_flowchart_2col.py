#!/usr/bin/env python
"""PanSynteny Steps 1-3 workflow - two-column (L-shaped) layout.

Steps 1 and 2 run down the LEFT column; a diagonal connector carries the flow
from the bottom-left (end of Step 2) up to the TOP-RIGHT, where Step 3 runs down
the RIGHT column. Alternative layout to the single-spine walkthrough_flowchart.dot
(both are kept).

Renders:  docs/figures/walkthrough_flowchart_2col.pdf  and  .png
Run:       python docs/figures/walkthrough_flowchart_2col.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---- palette (mirrors walkthrough_flowchart.dot) ----
C_PROC = "#FFFFFF"   # process box
C_IN   = "#EDE7F6"   # input
C_DATA = "#E3F2FD"   # intermediate data
C_OUT  = "#FFF3E0"   # L1/L2/L3 outputs
C_EDGE = "#444444"
C_TXT  = "#222222"
P_S1   = "#F3F6FB"; P_S1E = "#9DB2CE"   # step-1 panel
P_S2   = "#F2FAF4"; P_S2E = "#9CC7A8"   # step-2 panel
P_S3   = "#FDF6F0"; P_S3E = "#D8B79A"   # step-3 panel
C_LINK = "#E65100"                       # diagonal step-2 -> step-3 connector

fig, ax = plt.subplots(figsize=(10.6, 8.6))
ax.set_xlim(0, 13.8); ax.set_ylim(3.0, 12.95); ax.axis("off")

cxL, wL = 3.0, 4.2     # left column centre / box width
cxR, wR = 10.3, 4.8    # right column centre / box width
BH = 0.62              # box height

def box(cx, cy, text, fill, w, edge=C_EDGE, dashed=False, fcolor=C_TXT, h=BH):
    p = FancyBboxPatch((cx-w/2, cy-h/2), w, h,
                       boxstyle="round,pad=0.015,rounding_size=0.09",
                       facecolor=fill, edgecolor=edge, linewidth=1.1,
                       linestyle="--" if dashed else "-", zorder=3)
    ax.add_patch(p)
    ax.text(cx, cy, text, ha="center", va="center", fontsize=8.0,
            color=fcolor, zorder=4)

def panel(x0, x1, y0, y1, fill, edge, label):
    ax.add_patch(FancyBboxPatch((x0, y0), x1-x0, y1-y0,
                 boxstyle="round,pad=0.02,rounding_size=0.14",
                 facecolor=fill, edgecolor=edge, linewidth=1.2, zorder=0))
    ax.text(x0+0.18, y1-0.26, label, fontsize=10.5, fontweight="bold",
            color=C_TXT, ha="left", va="center", zorder=1)

def varrow(cx, y_top, y_bot, color=C_EDGE, lw=1.1):
    ax.annotate("", xy=(cx, y_bot+BH/2), xytext=(cx, y_top-BH/2),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                shrinkA=1, shrinkB=1), zorder=2)

# ========================= PANELS =========================
panel(0.7, 5.3, 8.05, 11.0, P_S1, P_S1E, "Step 1  -  Per-focal neighborhoods")
panel(0.7, 5.3, 4.35, 7.45, P_S2, P_S2E, "Step 2  -  Per-genome assembly")
panel(7.75, 12.85, 4.05, 12.3, P_S3, P_S3E, "Step 3  -  Cross-genome consolidation")

# ========================= LEFT COLUMN =========================
# inputs (two side by side, above step-1 panel)
box(1.85, 12.1, "Per-focal\nneighbor TSVs", C_IN, 2.1, edge="#5E35B1")
box(4.15, 12.1, "Catalog\nc80 tables",      C_IN, 2.1, edge="#5E35B1")
# step 1
box(cxL, 10.2, "Extract focal windows;\nlabel neighbors by c80", C_PROC, wL)
box(cxL,  9.4, "Operon size = mode of the\nper-genome size distribution", C_PROC, wL)
box(cxL,  8.6, "gene_neighbors\n(focal x genome x position)", C_DATA, wL)
# step 2
box(cxL,  6.75, "Pool all focal windows per\ngenome into one directed graph", C_PROC, wL)
box(cxL,  5.95, "DFS maximal paths\n(overlap assembly)", C_PROC, wL)
box(cxL,  5.15, "path_df\n(per-genome maximal path)", C_DATA, wL)

# left-column arrows
ax.annotate("", xy=(cxL, 10.2+BH/2), xytext=(1.85, 12.1-0.31),
            arrowprops=dict(arrowstyle="-|>", color=C_EDGE, lw=1.0), zorder=2)
ax.annotate("", xy=(cxL, 10.2+BH/2), xytext=(4.15, 12.1-0.31),
            arrowprops=dict(arrowstyle="-|>", color=C_EDGE, lw=1.0), zorder=2)
varrow(cxL, 10.2, 9.4)
varrow(cxL, 9.4, 8.6)
varrow(cxL, 8.6, 6.75)   # step1 -> step2 (across panel gap)
varrow(cxL, 6.75, 5.95)
varrow(cxL, 5.95, 5.15)

# ========================= DIAGONAL CONNECTOR =========================
# bottom-left (end of Step 2) -> top-right (start of Step 3)
ax.annotate("", xy=(7.95, 11.55), xytext=(5.15, 4.95),
            arrowprops=dict(arrowstyle="-|>", color=C_LINK, lw=2.4,
                            shrinkA=3, shrinkB=4,
                            connectionstyle="arc3,rad=-0.12"), zorder=5)
ax.text(6.05, 8.45, "across\ngenomes", ha="center", va="center", fontsize=9,
        style="italic", color=C_LINK, zorder=6,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none"))

# ========================= RIGHT COLUMN (Step 3) =========================
yR = [11.55, 10.75, 9.95, 9.15, 8.35, 7.55, 6.75, 5.95, 5.15]
box(cxR, yR[0], "Collapse identical path shapes", C_PROC, wR)
box(cxR, yR[1], "Canonicalize direction (lex-min):\nmerge fwd + reverse; n_genomes cut", C_PROC, wR)
box(cxR, yR[2], "Build joint components + project", C_PROC, wR)
box(cxR, yR[3], "Orient within component\n(to longest reference)", C_PROC, wR)
box(cxR, yR[4], "L1  canonical operons", C_OUT, wR, edge="#E65100")
box(cxR, yR[5], "Expand length-variant isoforms", C_PROC, wR)
box(cxR, yR[6], "L2  per-isoform", C_OUT, wR, edge="#E65100")
box(cxR, yR[7], "Per-genome expansion (needs_flip)\n+ truncation / fragmentation flags", C_PROC, wR)
box(cxR, yR[8], "L3  per-genome", C_OUT, wR, edge="#E65100")
for a, b in zip(yR[:-1], yR[1:]):
    varrow(cxR, a, b)

# output (below the Step-3 panel)
box(cxR, 3.55, "Downstream: gggenes figures (Step 5),\ntrait-block extraction (Step 6)",
    "#F5F5F5", wR, edge="#777777", dashed=True, fcolor="#555555")
varrow(cxR, yR[8], 3.55)

plt.tight_layout(pad=0.3)
out = "docs/figures/walkthrough_flowchart_2col"
fig.savefig(out + ".pdf", bbox_inches="tight")
fig.savefig(out + ".png", dpi=300, bbox_inches="tight")
print("wrote", out + ".pdf", "and", out + ".png")

#!/usr/bin/env python3
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: miRNA_figures.py
# Description: Generates publication figures for microRNA target interaction networks.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

S5 = "CuminGRdb_Table_S5_miRNA_Level_Summary.csv"
S6 = "CuminGRdb_Table_S6_Target_Gene_Level_Summary.csv"
S8 = "CuminGRdb_Table_S8_Expectation_Threshold_Summary.csv"
INHIBITION = "CuminGRdb_miRNA_Inhibition_Mode_Summary.csv"


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUTDIR = Path("CuminGRdb_miRNA_Figures_600dpi")
OUTDIR.mkdir(exist_ok=True)


# ============================================================
# PUBLICATION SETTINGS
# ============================================================

mpl.rcParams.update({

    "font.family": "DejaVu Sans",

    "font.size": 10,

    "axes.labelsize": 10,

    "axes.titlesize": 12,

    "axes.linewidth": 0.9,

    "xtick.labelsize": 9,

    "ytick.labelsize": 9,

    "xtick.major.width": 0.8,

    "ytick.major.width": 0.8,

    "pdf.fonttype": 42,

    "ps.fonttype": 42
})


# ============================================================
# PUBLICATION COLOR PALETTE
# ============================================================

COLORS = {

    "teal": "#2A9D8F",

    "blue": "#4575B4",

    "orange": "#E76F51",

    "purple": "#7B6AA8",

    "gold": "#E9C46A",

    "green": "#5A9367",

    "red": "#C8553D",

    "dark": "#333333"
}


# ============================================================
# LOAD DATA
# ============================================================

mirna = pd.read_csv(S5)

target = pd.read_csv(S6)

threshold = pd.read_csv(S8)

inhibition = pd.read_csv(INHIBITION)


# ============================================================
# BASIC CLEANING
# ============================================================

threshold["Threshold"] = (
    threshold["Expectation_Threshold"]
    .str.replace("<=", "", regex=False)
    .astype(float)
)

threshold = threshold.sort_values("Threshold")


# ============================================================
# FUNCTION TO SAVE THREE FORMATS
# ============================================================

def save_figure(fig, filename):

    base = OUTDIR / filename

    # PNG
    fig.savefig(
        f"{base}.png",
        dpi=600,
        bbox_inches="tight",
        facecolor="white"
    )

    # JPG
    fig.savefig(
        f"{base}.jpg",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"quality": 95}
    )

    # TIFF
    fig.savefig(
        f"{base}.tiff",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"}
    )

    print(f"Saved: {base}.png")
    print(f"Saved: {base}.jpg")
    print(f"Saved: {base}.tiff")


# ============================================================
# FIGURE 1
# EXPECTATION-SCORE THRESHOLD DISTRIBUTION
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.0, 5.2)
)

bars = ax.bar(
    threshold["Threshold"].astype(str),
    threshold["Interactions"],
    color=[
        COLORS["teal"],
        COLORS["blue"],
        COLORS["green"],
        COLORS["orange"],
        COLORS["purple"]
    ],
    edgecolor="black",
    linewidth=0.6
)

ax.set_xlabel(
    "psRNATarget expectation-score threshold"
)

ax.set_ylabel(
    "Predicted interactions"
)

ax.set_title(
    "Distribution of predicted miRNA–target interactions",
    fontweight="bold"
)

ax.grid(
    axis="y",
    linestyle="--",
    linewidth=0.5,
    alpha=0.3
)

ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bar, value, percentage in zip(
    bars,
    threshold["Interactions"],
    threshold["Percentage_of_All_Interactions"]
):

    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height(),
        f"{value:,}\n({percentage:.2f}%)",
        ha="center",
        va="bottom",
        fontsize=8
    )

ax.set_ylim(
    0,
    threshold["Interactions"].max() * 1.15
)

save_figure(
    fig,
    "Figure_miRNA_Expectation_Thresholds"
)

plt.close(fig)


# ============================================================
# FIGURE 2
# INHIBITION MECHANISM
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.0, 4.5)
)

inhibition = inhibition.sort_values(
    "Percentage",
    ascending=True
)

bars = ax.barh(
    inhibition["Inhibition"],
    inhibition["Percentage"],
    color=[
        COLORS["purple"],
        COLORS["orange"]
    ],
    edgecolor="black",
    linewidth=0.6
)

ax.set_xlabel(
    "Predicted interactions (%)"
)

ax.set_ylabel(
    "Inhibition mode"
)

ax.set_title(
    "Predicted miRNA inhibition mechanisms",
    fontweight="bold"
)

ax.set_xlim(
    0,
    inhibition["Percentage"].max() * 1.18
)

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.5,
    alpha=0.3
)

ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bar, value in zip(
    bars,
    inhibition["Percentage"]
):

    ax.text(
        value + 1,
        bar.get_y() + bar.get_height() / 2,
        f"{value:.2f}%",
        va="center",
        fontsize=9
    )

save_figure(
    fig,
    "Figure_miRNA_Inhibition_Mechanisms"
)

plt.close(fig)


# ============================================================
# FIGURE 3
# TOP 20 miRNAs BY NUMBER OF TARGET GENES
# ============================================================

top_mirnas = (
    mirna
    .sort_values(
        "Target_Genes",
        ascending=False
    )
    .head(20)
    .sort_values(
        "Target_Genes"
    )
)

fig, ax = plt.subplots(
    figsize=(8.0, 7.0)
)

bars = ax.barh(
    top_mirnas["miRNA_Acc."],
    top_mirnas["Target_Genes"],
    color=COLORS["teal"],
    edgecolor="black",
    linewidth=0.5
)

ax.set_xlabel(
    "Number of unique target genes"
)

ax.set_ylabel(
    "miRNA"
)

ax.set_title(
    "Top 20 miRNAs by predicted target-gene breadth",
    fontweight="bold"
)

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.5,
    alpha=0.3
)

ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bar, value in zip(
    bars,
    top_mirnas["Target_Genes"]
):

    ax.text(
        value + 2,
        bar.get_y() + bar.get_height() / 2,
        f"{value:,}",
        va="center",
        fontsize=8
    )

ax.set_xlim(
    0,
    top_mirnas["Target_Genes"].max() * 1.12
)

save_figure(
    fig,
    "Figure_miRNA_Top20_Target_Breadth"
)

plt.close(fig)


# ============================================================
# FIGURE 4
# TOP 20 TARGET GENES BY NUMBER OF TARGETING miRNAs
# ============================================================

top_targets = (
    target
    .sort_values(
        "Targeting_miRNAs",
        ascending=False
    )
    .head(20)
    .sort_values(
        "Targeting_miRNAs"
    )
)

fig, ax = plt.subplots(
    figsize=(8.0, 7.0)
)

bars = ax.barh(
    top_targets["Target_Acc."],
    top_targets["Targeting_miRNAs"],
    color=COLORS["purple"],
    edgecolor="black",
    linewidth=0.5
)

ax.set_xlabel(
    "Number of targeting miRNAs"
)

ax.set_ylabel(
    "Target gene"
)

ax.set_title(
    "Top 20 genes by predicted miRNA regulation",
    fontweight="bold"
)

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.5,
    alpha=0.3
)

ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bar, value in zip(
    bars,
    top_targets["Targeting_miRNAs"]
):

    ax.text(
        value + 5,
        bar.get_y() + bar.get_height() / 2,
        f"{value:,}",
        va="center",
        fontsize=8
    )

ax.set_xlim(
    0,
    top_targets["Targeting_miRNAs"].max() * 1.12
)

save_figure(
    fig,
    "Figure_miRNA_Top20_Target_Genes"
)

plt.close(fig)


# ============================================================
# FIGURE 5
# HIGH-CONFIDENCE INTERACTIONS
# ============================================================

fig, ax = plt.subplots(
    figsize=(7.5, 5.0)
)

x = threshold["Threshold"]

y = threshold["Percentage_of_All_Interactions"]

ax.plot(
    x,
    y,
    marker="o",
    markersize=7,
    linewidth=2,
    color=COLORS["red"]
)

ax.fill_between(
    x,
    y,
    alpha=0.15,
    color=COLORS["red"]
)

ax.set_xlabel(
    "Expectation-score threshold"
)

ax.set_ylabel(
    "Interactions retained (%)"
)

ax.set_title(
    "Cumulative retention of predicted interactions",
    fontweight="bold"
)

ax.set_xticks(x)

ax.grid(
    linestyle="--",
    linewidth=0.5,
    alpha=0.3
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for xi, yi in zip(x, y):

    ax.annotate(
        f"{yi:.2f}%",
        (xi, yi),
        xytext=(0, 9),
        textcoords="offset points",
        ha="center",
        fontsize=8
    )

save_figure(
    fig,
    "Figure_miRNA_Cumulative_Confidence"
)

plt.close(fig)


# ============================================================
# FINISHED
# ============================================================

print("\n==============================================")
print("All publication figures generated")
print("Resolution: 600 dpi")
print("Formats: PNG, JPG, TIFF")
print("==============================================")

#!/usr/bin/env python3
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: ssr_plots.py
# Description: Generates SSR motif frequency and repeat length distribution plots.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

INPUT_FILE = "CuminGRdb_Table_S3_SSR_Summary.csv"

# ============================================================
# PUBLICATION SETTINGS
# ============================================================

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.linewidth"] = 1.0

mpl.rcParams["xtick.major.width"] = 0.8
mpl.rcParams["ytick.major.width"] = 0.8

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# ============================================================
# PUBLICATION COLOR PALETTE
# ============================================================
#
# One restrained color for each biological category.
# These are deliberately not overly saturated.
#

COLOR_REPEAT = "#2A9D8F"      # Teal
COLOR_LOCATION = "#E76F51"    # Coral/orange
COLOR_PRIMER = "#6C5B9E"      # Purple

# ============================================================
# READ DATA
# ============================================================

df = pd.read_csv(INPUT_FILE)

# ============================================================
# PANEL A — SSR REPEAT CLASSES
# ============================================================

repeat_classes = [
    "Mononucleotide",
    "Dinucleotide",
    "Trinucleotide",
    "Tetranucleotide",
    "Pentanucleotide",
    "Hexanucleotide"
]

repeat_df = df[
    df["Feature"].isin(repeat_classes)
].copy()

repeat_df["Feature"] = pd.Categorical(
    repeat_df["Feature"],
    categories=repeat_classes,
    ordered=True
)

repeat_df = repeat_df.sort_values("Feature")

# ============================================================
# PANEL B — GENOMIC DISTRIBUTION
# ============================================================

location_features = [
    "Genic",
    "Intergenic"
]

location_df = df[
    df["Feature"].isin(location_features)
].copy()

# ============================================================
# PANEL C — PRIMER DESIGN
# ============================================================

primer_features = [
    "SSR loci with successful primer pairs",
    "SSR loci without successful primer pairs"
]

primer_df = df[
    df["Feature"].isin(primer_features)
].copy()

# ============================================================
# FIGURE
# ============================================================

fig, axes = plt.subplots(
    1,
    3,
    figsize=(16, 6.5)
)

# ============================================================
# PANEL A
# SSR REPEAT CLASSES
# ============================================================

ax = axes[0]

bars = ax.bar(
    repeat_df["Feature"],
    repeat_df["Number_of_Loci"],
    color=COLOR_REPEAT,
    edgecolor="black",
    linewidth=0.6
)

ax.set_title(
    "A. SSR repeat classes",
    fontsize=12,
    fontweight="bold",
    pad=10
)

ax.set_ylabel(
    "Number of SSR loci",
    fontsize=10
)

ax.set_xlabel(
    "Repeat class",
    fontsize=10
)

ax.tick_params(
    axis="x",
    rotation=45,
    labelsize=8.5
)

ax.grid(
    axis="y",
    linestyle="--",
    linewidth=0.5,
    alpha=0.30
)

ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

max_repeat = repeat_df["Number_of_Loci"].max()

for bar, count in zip(
    bars,
    repeat_df["Number_of_Loci"]
):

    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + max_repeat * 0.015,
        f"{int(count):,}",
        ha="center",
        va="bottom",
        fontsize=8
    )

ax.set_ylim(
    0,
    max_repeat * 1.13
)

# ============================================================
# PANEL B
# GENOMIC DISTRIBUTION
# ============================================================

ax = axes[1]

bars = ax.bar(
    location_df["Feature"],
    location_df["Number_of_Loci"],
    color=COLOR_LOCATION,
    edgecolor="black",
    linewidth=0.6
)

ax.set_title(
    "B. Genomic distribution",
    fontsize=12,
    fontweight="bold",
    pad=10
)

ax.set_ylabel(
    "Number of SSR loci",
    fontsize=10
)

ax.set_xlabel(
    "Genomic location",
    fontsize=10
)

ax.tick_params(
    axis="x",
    labelsize=9
)

ax.grid(
    axis="y",
    linestyle="--",
    linewidth=0.5,
    alpha=0.30
)

ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

max_location = location_df["Number_of_Loci"].max()

for bar, count in zip(
    bars,
    location_df["Number_of_Loci"]
):

    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + max_location * 0.015,
        f"{int(count):,}",
        ha="center",
        va="bottom",
        fontsize=9
    )

ax.set_ylim(
    0,
    max_location * 1.13
)

# ============================================================
# PANEL C
# PRIMER DESIGN
# ============================================================

ax = axes[2]

bars = ax.bar(
    primer_df["Feature"],
    primer_df["Percentage_of_Total_SSRs"],
    color=COLOR_PRIMER,
    edgecolor="black",
    linewidth=0.6
)

ax.set_title(
    "C. Primer-design success",
    fontsize=12,
    fontweight="bold",
    pad=10
)

ax.set_ylabel(
    "SSR loci (%)",
    fontsize=10
)

ax.set_xlabel(
    "Primer-design outcome",
    fontsize=10
)

ax.tick_params(
    axis="x",
    rotation=25,
    labelsize=8.5
)

ax.grid(
    axis="y",
    linestyle="--",
    linewidth=0.5,
    alpha=0.30
)

ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

max_primer = primer_df["Percentage_of_Total_SSRs"].max()

for bar, value in zip(
    bars,
    primer_df["Percentage_of_Total_SSRs"]
):

    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + max_primer * 0.025,
        f"{value:.2f}%",
        ha="center",
        va="bottom",
        fontsize=9
    )

ax.set_ylim(
    0,
    max_primer * 1.18
)

# ============================================================
# MAIN TITLE
# ============================================================

fig.suptitle(
    r"Genome-wide simple sequence repeat resource in "
    r"$\it{Cuminum\ cyminum}$",
    fontsize=15,
    fontweight="bold",
    y=1.02
)

# ============================================================
# LAYOUT
# ============================================================

plt.tight_layout()

# ============================================================
# SAVE — 600 DPI
# ============================================================

plt.savefig(
    "CuminGRdb_SSR_Resource.png",
    dpi=600,
    format="png",
    bbox_inches="tight"
)

plt.savefig(
    "CuminGRdb_SSR_Resource.jpg",
    dpi=600,
    format="jpg",
    bbox_inches="tight"
)

plt.savefig(
    "CuminGRdb_SSR_Resource.tiff",
    dpi=600,
    format="tiff",
    bbox_inches="tight",
    pil_kwargs={"compression": "tiff_lzw"}
)

plt.show()

print("\n==============================================")
print("Publication figure generated")
print("==============================================")
print("PNG  : CuminGRdb_SSR_Resource.png")
print("JPG  : CuminGRdb_SSR_Resource.jpg")
print("TIFF : CuminGRdb_SSR_Resource.tiff")
print("Resolution: 600 dpi")
print("==============================================")

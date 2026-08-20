#!/usr/bin/env python3
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: new_plant_tf.py
# Description: Transcription factor family classification and TAIR ortholog distribution plotter.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

INPUT_FILE = "PlantTFDB.txt"

# Plot all TF families
TOP_N = 50

OUTPUT_PREFIX = f"CuminGRdb_TF_All{TOP_N}"

# ============================================================
# PUBLICATION SETTINGS
# ============================================================

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 10

mpl.rcParams["axes.linewidth"] = 1.0

mpl.rcParams["xtick.major.width"] = 0.8
mpl.rcParams["ytick.major.width"] = 0.8

mpl.rcParams["xtick.minor.width"] = 0.6
mpl.rcParams["ytick.minor.width"] = 0.6

# ============================================================
# READ PlantTFDB FILE
# ============================================================

records = []

with open(INPUT_FILE, "r", encoding="utf-8") as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        fields = line.split()

        if len(fields) < 2:
            continue

        gene_id = fields[0]
        tf_family = fields[1]

        records.append(
            [gene_id, tf_family]
        )

df = pd.DataFrame(
    records,
    columns=["GeneID", "TF_Family"]
)

# ============================================================
# REMOVE MALFORMED TF-FAMILY RECORDS
# ============================================================
#
# These are records where a CcGene identifier was incorrectly
# placed in the TF-family column.
#
# ============================================================

df = df[
    ~df["TF_Family"].str.startswith("CcGene_")
].copy()

# ============================================================
# COUNT TF FAMILIES
# ============================================================

family_counts = (
    df["TF_Family"]
    .value_counts()
    .reset_index()
)

family_counts.columns = [
    "TF_Family",
    "Gene_Count"
]

# ============================================================
# SELECT ALL TF FAMILIES
# ============================================================

plot_df = (
    family_counts
    .head(TOP_N)
    .sort_values(
        "Gene_Count",
        ascending=False
    )
    .reset_index(drop=True)
)

# ============================================================
# PRINT DATA
# ============================================================

print("\nTF families used for plotting:\n")

print(
    plot_df.to_string(index=False)
)

print("\nTotal genes represented:",
      plot_df["Gene_Count"].sum())

# ============================================================
# COLOUR PALETTE
# ============================================================
#
# A restrained categorical palette suitable for scientific
# figures. Colours are repeated only if necessary.
#
# ============================================================

colors = [
    "#4C78A8",
    "#F58518",
    "#54A24B",
    "#E45756",
    "#72B7B2",
    "#B279A2",
    "#FF9DA6",
    "#9D755D",
    "#BAB0AC",
    "#2E86AB",
    "#F2A541",
    "#5C946E",
    "#C44E52",
    "#8172B3",
    "#64B5CD",
    "#CC8C3C",
    "#6B8E23",
    "#A05195",
    "#D45087",
    "#3F8EAA",
]

# Repeat palette if required for all 50 families
bar_colors = [
    colors[i % len(colors)]
    for i in range(len(plot_df))
]

# ============================================================
# FIGURE SIZE
# ============================================================

fig, ax = plt.subplots(
    figsize=(10, 13)
)

# ============================================================
# BAR PLOT
# ============================================================
#
# IMPORTANT:
# ascending=False means:
# highest value = top
# lowest value  = bottom
#
# ============================================================

y_positions = np.arange(len(plot_df))

bars = ax.barh(
    y_positions,
    plot_df["Gene_Count"],
    color=bar_colors,
    edgecolor="black",
    linewidth=0.45,
    height=0.72
)

# ============================================================
# Y-AXIS LABELS
# ============================================================

ax.set_yticks(y_positions)

ax.set_yticklabels(
    plot_df["TF_Family"],
    fontsize=9
)

# ============================================================
# CRITICAL: HIGHEST AT TOP
# ============================================================

ax.invert_yaxis()

# ============================================================
# X-AXIS
# ============================================================

ax.set_xlabel(
    "Number of transcription factor genes",
    fontsize=11,
    labelpad=8
)

# ============================================================
# Y-AXIS
# ============================================================

ax.set_ylabel(
    "Transcription factor family",
    fontsize=11,
    labelpad=8
)

# ============================================================
# TITLE
# ============================================================

#ax.set_title(
#    r"Transcription factor families in "
#    r"$\it{Cuminum\ cyminum}$",
#    fontsize=14,
#    fontweight="bold",
#    pad=14
#)

# ============================================================
# GRID
# ============================================================

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.55,
    alpha=0.30
)

ax.set_axisbelow(True)

# ============================================================
# CLEAN SPINES
# ============================================================

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Keep left and bottom axes
ax.spines["left"].set_linewidth(0.9)
ax.spines["bottom"].set_linewidth(0.9)

# ============================================================
# TICKS
# ============================================================

ax.tick_params(
    axis="y",
    labelsize=9,
    length=3,
    width=0.8
)

ax.tick_params(
    axis="x",
    labelsize=9,
    length=3,
    width=0.8
)

# ============================================================
# ADD GENE COUNTS
# ============================================================

max_count = plot_df["Gene_Count"].max()

for bar, count in zip(
    bars,
    plot_df["Gene_Count"]
):

    ax.text(
        count + max_count * 0.012,
        bar.get_y() + bar.get_height() / 2,
        f"{int(count):,}",
        va="center",
        ha="left",
        fontsize=8.5
    )

# ============================================================
# X-AXIS LIMIT
# ============================================================

ax.set_xlim(
    0,
    max_count * 1.13
)

# ============================================================
# LAYOUT
# ============================================================

plt.tight_layout()

# ============================================================
# EXPORT — 600 DPI PNG
# ============================================================

plt.savefig(
    f"{OUTPUT_PREFIX}.png",
    dpi=600,
    format="png",
    bbox_inches="tight"
)

# ============================================================
# EXPORT — 600 DPI JPG
# ============================================================

plt.savefig(
    f"{OUTPUT_PREFIX}.jpg",
    dpi=600,
    format="jpg",
    bbox_inches="tight"
)

# ============================================================
# EXPORT — 600 DPI TIFF
# ============================================================

plt.savefig(
    f"{OUTPUT_PREFIX}.tiff",
    dpi=600,
    format="tiff",
    bbox_inches="tight",
    pil_kwargs={
        "compression": "tiff_lzw"
    }
)

# ============================================================
# DISPLAY
# ============================================================

plt.show()

# ============================================================
# SUMMARY
# ============================================================

print("\n==============================================")
print("Publication figures generated successfully")
print("==============================================")

print(f"Total TF families : {len(family_counts)}")
print(f"Families plotted  : {len(plot_df)}")
print(f"Largest family    : {plot_df.iloc[0]['TF_Family']} "
      f"({plot_df.iloc[0]['Gene_Count']})")
print(f"Smallest family   : {plot_df.iloc[-1]['TF_Family']} "
      f"({plot_df.iloc[-1]['Gene_Count']})")

print("\nOutput files:")

print(f"{OUTPUT_PREFIX}.png")
print(f"{OUTPUT_PREFIX}.jpg")
print(f"{OUTPUT_PREFIX}.tiff")

print("==============================================")

#!/usr/bin/env python3
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: plant_tf.py
# Description: Analyzes transcription factor family distributions across 50 plant families.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

INPUT_FILE = "PlantTFDB.txt"

# Choose 20 or 30
TOP_N = 30

OUTPUT_PREFIX = f"CuminGRdb_TF_Top{TOP_N}"

# ============================================================
# PUBLICATION SETTINGS
# ============================================================

mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.linewidth"] = 1.0

mpl.rcParams["xtick.major.width"] = 0.8
mpl.rcParams["ytick.major.width"] = 0.8

# ============================================================
# READ PlantTFDB
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
# REMOVE MALFORMED FAMILY RECORDS
#
# Some records contain a CcGene ID in the TF-family field.
# These are excluded from the family distribution.
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
# SELECT TOP N
# ============================================================

plot_df = (
    family_counts
    .head(TOP_N)
    .sort_values(
        "Gene_Count",
        ascending=True
    )
)

# ============================================================
# PRINT DATA
# ============================================================

print("\nTF families used for plotting:\n")

print(
    plot_df
    .sort_values(
        "Gene_Count",
        ascending=False
    )
    .to_string(index=False)
)

# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 10)
)

bars = ax.barh(
    plot_df["TF_Family"],
    plot_df["Gene_Count"],
    edgecolor="black",
    linewidth=0.6
)

# ============================================================
# LABELS
# ============================================================

ax.set_xlabel(
    "Number of transcription factor genes",
    fontsize=11
)

ax.set_ylabel(
    "Transcription factor family",
    fontsize=11
)

ax.set_title(
    f"Top {TOP_N} transcription factor families in "
    r"$\it{Cuminum\ cyminum}$",
    fontsize=14,
    fontweight="bold",
    pad=12
)

# ============================================================
# GRID
# ============================================================

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.6,
    alpha=0.35
)

ax.set_axisbelow(True)

# ============================================================
# CLEAN SPINES
# ============================================================

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ============================================================
# TICKS
# ============================================================

ax.tick_params(
    axis="y",
    labelsize=9.5
)

ax.tick_params(
    axis="x",
    labelsize=9
)

# ============================================================
# ADD COUNTS
# ============================================================

max_count = plot_df["Gene_Count"].max()

for bar, count in zip(
    bars,
    plot_df["Gene_Count"]
):

    ax.text(
        count + max_count * 0.015,
        bar.get_y() + bar.get_height() / 2,
        f"{int(count):,}",
        va="center",
        ha="left",
        fontsize=9
    )

ax.set_xlim(
    0,
    max_count * 1.12
)

# ============================================================
# LAYOUT
# ============================================================

plt.tight_layout()

# ============================================================
# EXPORT — 600 DPI
# ============================================================

plt.savefig(
    f"{OUTPUT_PREFIX}.png",
    dpi=600,
    format="png",
    bbox_inches="tight"
)

plt.savefig(
    f"{OUTPUT_PREFIX}.jpg",
    dpi=600,
    format="jpg",
    bbox_inches="tight"
)

plt.savefig(
    f"{OUTPUT_PREFIX}.tiff",
    dpi=600,
    format="tiff",
    bbox_inches="tight",
    pil_kwargs={"compression": "tiff_lzw"}
)

plt.show()

# ============================================================
# SUMMARY
# ============================================================

print("\n==========================================")
print("Publication figure generated successfully")
print("==========================================")
print(f"Total TF families : {len(family_counts)}")
print(f"Families plotted  : {len(plot_df)}")
print("")
print(f"{OUTPUT_PREFIX}.png")
print(f"{OUTPUT_PREFIX}.jpg")
print(f"{OUTPUT_PREFIX}.tiff")
print("==========================================")

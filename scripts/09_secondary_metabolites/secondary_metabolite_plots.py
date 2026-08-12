#!/usr/bin/env python3
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: secondary_metabolite_plots.py
# Description: Generates secondary metabolite category and KEGG pathway distribution plots.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

INPUT = "secondary_metabolites.csv"
OUTDIR = "secondary_metabolite_figures"

DPI = 600

# Publication-friendly palette.
# Explicit colors are used because a multi-category figure
# needs clear visual separation between functional classes.
COLORS = [
    "#3B82F6",  # blue
    "#F59E0B",  # amber
    "#10B981",  # green
    "#EF4444",  # red
    "#8B5CF6",  # purple
    "#06B6D4",  # cyan
    "#EC4899",  # pink
    "#84CC16",  # lime
    "#F97316",  # orange
    "#6366F1",  # indigo
]


# ============================================================
# Output directory
# ============================================================

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# Publication plotting settings
# ============================================================

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "figure.dpi": 150,
    "savefig.dpi": DPI,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
})


# ============================================================
# Helper function
# ============================================================

def save_figure(fig, basename):
    """
    Save the same figure in PNG, TIFF and JPEG at 600 dpi.
    TIFF uses LZW compression for a smaller lossless file.
    """

    png = os.path.join(OUTDIR, basename + ".png")
    tiff = os.path.join(OUTDIR, basename + ".tiff")
    jpg = os.path.join(OUTDIR, basename + ".jpg")

    fig.savefig(
        png,
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white"
    )

    fig.savefig(
        tiff,
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"}
    )

    fig.savefig(
        jpg,
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"quality": 95}
    )

    plt.close(fig)


# ============================================================
# Read data
# ============================================================

df = pd.read_csv(INPUT)

required_columns = [
    "Gene_ID",
    "Metabolite_Category",
    "Description",
    "SwissProt_Hit",
    "NR_Hit",
    "GOs",
    "KEGG_Pathway"
]

missing = [c for c in required_columns if c not in df.columns]

if missing:
    raise ValueError(
        f"Missing required columns: {', '.join(missing)}"
    )

print("\n============================================================")
print("CuminGRdb secondary metabolism annotation analysis")
print("============================================================")

print(f"Total records: {len(df):,}")


# ============================================================
# Clean missing values
# ============================================================

for col in required_columns:
    df[col] = df[col].fillna("-").astype(str).str.strip()


# ============================================================
# A. Metabolite category distribution
# ============================================================

category_counts = (
    df["Metabolite_Category"]
    .replace("-", np.nan)
    .dropna()
    .value_counts()
    .sort_values(ascending=True)
)

print("\nMetabolite categories:")
print(category_counts.sort_values(ascending=False))


fig, ax = plt.subplots(figsize=(7.2, 5.2))

y = np.arange(len(category_counts))

bars = ax.barh(
    y,
    category_counts.values,
    color=[
        COLORS[i % len(COLORS)]
        for i in range(len(category_counts))
    ],
    edgecolor="black",
    linewidth=0.5
)

ax.set_yticks(y)
ax.set_yticklabels(category_counts.index)

ax.set_xlabel("Number of genes")
ax.set_ylabel("Metabolite-associated category")

ax.set_title(
    "Secondary metabolism-associated gene categories",
    fontweight="bold"
)

ax.grid(
    axis="x",
    linestyle="--",
    linewidth=0.5,
    alpha=0.25
)

ax.set_axisbelow(True)

for bar, value in zip(bars, category_counts.values):
    ax.text(
        value + max(category_counts.values) * 0.015,
        bar.get_y() + bar.get_height() / 2,
        f"{value:,}",
        va="center",
        fontsize=8
    )

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()

save_figure(
    fig,
    "SecondaryMetabolism_Category"
)


# ============================================================
# B. KEGG pathway analysis
# ============================================================

pathway_records = []

for _, row in df.iterrows():

    value = row["KEGG_Pathway"]

    if value in ["-", "", "nan"]:
        continue

    pathways = value.split(",")

    for pathway in pathways:

        pathway = pathway.strip()

        if not pathway:
            continue

        # Remove map/ko duplication.
        # ko00904 and map00904 represent the same KEGG pathway.
        pathway = re.sub(r"^map", "", pathway)
        pathway = re.sub(r"^ko", "", pathway)

        if pathway.isdigit():
            pathway_records.append(
                (row["Gene_ID"], pathway)
            )

pathway_df = pd.DataFrame(
    pathway_records,
    columns=["Gene_ID", "KEGG_ID"]
)

if not pathway_df.empty:

    pathway_counts = (
        pathway_df
        .drop_duplicates()
        .groupby("KEGG_ID")["Gene_ID"]
        .nunique()
        .sort_values(ascending=False)
        .head(20)
        .sort_values()
    )

    print("\nTop KEGG pathways:")
    print(pathway_counts.sort_values(ascending=False))

    fig, ax = plt.subplots(figsize=(7.5, 6.2))

    y = np.arange(len(pathway_counts))

    bars = ax.barh(
        y,
        pathway_counts.values,
        color=COLORS[1],
        edgecolor="black",
        linewidth=0.5
    )

    ax.set_yticks(y)
    ax.set_yticklabels(
        ["KEGG " + x for x in pathway_counts.index]
    )

    ax.set_xlabel("Number of associated genes")
    ax.set_ylabel("KEGG pathway")

    ax.set_title(
        "Top KEGG pathways represented by secondary metabolism-associated genes",
        fontweight="bold"
    )

    ax.grid(
        axis="x",
        linestyle="--",
        linewidth=0.5,
        alpha=0.25
    )

    ax.set_axisbelow(True)

    for bar, value in zip(bars, pathway_counts.values):

        ax.text(
            value + max(pathway_counts.values) * 0.015,
            bar.get_y() + bar.get_height() / 2,
            f"{value:,}",
            va="center",
            fontsize=8
        )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()

    save_figure(
        fig,
        "SecondaryMetabolism_KEGG"
    )


# ============================================================
# C. Functional annotation evidence
# ============================================================

def has_annotation(series):
    return (
        ~series.isin(["-", "", "nan", "None"])
    )


evidence = {
    "Swiss-Prot": has_annotation(df["SwissProt_Hit"]).sum(),
    "NR": has_annotation(df["NR_Hit"]).sum(),
    "GO": has_annotation(df["GOs"]).sum(),
    "KEGG": has_annotation(df["KEGG_Pathway"]).sum(),
}

evidence_df = pd.DataFrame(
    list(evidence.items()),
    columns=["Annotation", "Genes"]
)

print("\nAnnotation evidence:")
print(evidence_df)


fig, ax = plt.subplots(figsize=(6.8, 4.8))

bars = ax.bar(
    evidence_df["Annotation"],
    evidence_df["Genes"],
    color=COLORS[:4],
    edgecolor="black",
    linewidth=0.6
)

ax.set_ylabel("Number of genes")
ax.set_xlabel("Annotation source")

ax.set_title(
    "Functional annotation evidence",
    fontweight="bold"
)

ax.grid(
    axis="y",
    linestyle="--",
    linewidth=0.5,
    alpha=0.25
)

ax.set_axisbelow(True)

for bar, value in zip(
    bars,
    evidence_df["Genes"]
):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + max(evidence_df["Genes"]) * 0.02,
        f"{value:,}",
        ha="center",
        va="bottom",
        fontsize=8
    )

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()

save_figure(
    fig,
    "SecondaryMetabolism_AnnotationEvidence"
)


# ============================================================
# D. GO annotation breadth
# ============================================================

def count_go_terms(value):

    if value in ["-", "", "nan"]:
        return 0

    return len([
        x for x in value.split(",")
        if x.strip()
    ])


df["GO_Count"] = df["GOs"].apply(count_go_terms)

go_counts = df.loc[
    df["GO_Count"] > 0,
    "GO_Count"
]

fig, ax = plt.subplots(figsize=(7, 4.8))

bins = np.arange(
    go_counts.min(),
    go_counts.max() + 2
)

ax.hist(
    go_counts,
    bins=bins,
    color=COLORS[2],
    edgecolor="black",
    linewidth=0.4
)

ax.set_xlabel("Number of GO terms per gene")
ax.set_ylabel("Number of genes")

ax.set_title(
    "Gene Ontology annotation breadth",
    fontweight="bold"
)

ax.grid(
    axis="y",
    linestyle="--",
    linewidth=0.5,
    alpha=0.25
)

ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()

save_figure(
    fig,
    "SecondaryMetabolism_GO_Breadth"
)


# ============================================================
# E. Combined publication figure
# ============================================================

fig = plt.figure(figsize=(11, 8.5))

gs = fig.add_gridspec(
    2,
    2,
    hspace=0.42,
    wspace=0.32
)


# ------------------------------------------------------------
# Panel A
# ------------------------------------------------------------

ax1 = fig.add_subplot(gs[0, 0])

y = np.arange(len(category_counts))

bars = ax1.barh(
    y,
    category_counts.values,
    color=[
        COLORS[i % len(COLORS)]
        for i in range(len(category_counts))
    ],
    edgecolor="black",
    linewidth=0.45
)

ax1.set_yticks(y)
ax1.set_yticklabels(category_counts.index)

ax1.set_xlabel("Number of genes")
ax1.set_title(
    "A. Metabolite-associated categories",
    loc="left",
    fontweight="bold"
)

ax1.grid(
    axis="x",
    linestyle="--",
    linewidth=0.4,
    alpha=0.25
)

ax1.set_axisbelow(True)

for bar, value in zip(bars, category_counts.values):
    ax1.text(
        value + max(category_counts.values) * 0.015,
        bar.get_y() + bar.get_height() / 2,
        str(value),
        va="center",
        fontsize=7
    )

ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)


# ------------------------------------------------------------
# Panel B
# ------------------------------------------------------------

ax2 = fig.add_subplot(gs[0, 1])

if not pathway_df.empty:

    bars = ax2.barh(
        np.arange(len(pathway_counts)),
        pathway_counts.values,
        color=COLORS[1],
        edgecolor="black",
        linewidth=0.45
    )

    ax2.set_yticks(
        np.arange(len(pathway_counts))
    )

    ax2.set_yticklabels(
        ["KEGG " + x for x in pathway_counts.index],
        fontsize=7
    )

    for bar, value in zip(
        bars,
        pathway_counts.values
    ):
        ax2.text(
            value + max(pathway_counts.values) * 0.015,
            bar.get_y() + bar.get_height() / 2,
            str(value),
            va="center",
            fontsize=7
        )

ax2.set_xlabel("Number of genes")

ax2.set_title(
    "B. Top KEGG pathways",
    loc="left",
    fontweight="bold"
)

ax2.grid(
    axis="x",
    linestyle="--",
    linewidth=0.4,
    alpha=0.25
)

ax2.set_axisbelow(True)

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)


# ------------------------------------------------------------
# Panel C
# ------------------------------------------------------------

ax3 = fig.add_subplot(gs[1, 0])

bars = ax3.bar(
    evidence_df["Annotation"],
    evidence_df["Genes"],
    color=COLORS[:4],
    edgecolor="black",
    linewidth=0.5
)

for bar, value in zip(
    bars,
    evidence_df["Genes"]
):
    ax3.text(
        bar.get_x() + bar.get_width() / 2,
        value + max(evidence_df["Genes"]) * 0.02,
        f"{value:,}",
        ha="center",
        fontsize=7
    )

ax3.set_ylabel("Number of genes")

ax3.set_title(
    "C. Annotation evidence",
    loc="left",
    fontweight="bold"
)

ax3.grid(
    axis="y",
    linestyle="--",
    linewidth=0.4,
    alpha=0.25
)

ax3.set_axisbelow(True)

ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)


# ------------------------------------------------------------
# Panel D
# ------------------------------------------------------------

ax4 = fig.add_subplot(gs[1, 1])

ax4.hist(
    go_counts,
    bins=np.arange(
        go_counts.min(),
        go_counts.max() + 2
    ),
    color=COLORS[2],
    edgecolor="black",
    linewidth=0.4
)

ax4.set_xlabel("GO terms per gene")
ax4.set_ylabel("Number of genes")

ax4.set_title(
    "D. GO annotation breadth",
    loc="left",
    fontweight="bold"
)

ax4.grid(
    axis="y",
    linestyle="--",
    linewidth=0.4,
    alpha=0.25
)

ax4.set_axisbelow(True)

ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)


# ------------------------------------------------------------
# Overall figure
# ------------------------------------------------------------

fig.suptitle(
    "Genome-wide secondary metabolism-associated gene resource in Cuminum cyminum",
    fontsize=13,
    fontweight="bold",
    y=0.98
)

save_figure(
    fig,
    "SecondaryMetabolism_Figure"
)


# ============================================================
# Save summary tables
# ============================================================

category_counts.sort_values(
    ascending=False
).rename(
    "Gene_Count"
).to_csv(
    os.path.join(
        OUTDIR,
        "SecondaryMetabolism_Category_Summary.csv"
    )
)

evidence_df.to_csv(
    os.path.join(
        OUTDIR,
        "SecondaryMetabolism_Annotation_Evidence.csv"
    ),
    index=False
)

if not pathway_df.empty:

    pathway_counts.sort_values(
        ascending=False
    ).rename(
        "Gene_Count"
    ).to_csv(
        os.path.join(
            OUTDIR,
            "SecondaryMetabolism_KEGG_Summary.csv"
        )
    )


print("\n============================================================")
print("Figures generated successfully")
print("============================================================")
print(f"Output directory: {OUTDIR}")
print("Resolution       : 600 dpi")
print("Formats          : PNG, JPEG, TIFF")
print("============================================================\n")

#!/usr/bin/env python3
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: miRNA_summary.py
# Description: Summarizes microRNA target statistics and expectation cutoffs (E <= 3.0).
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

INPUT_FILE = "psRNATargetJob-1785848631371569.txt"

# ============================================================
# READ psRNATarget OUTPUT
# ============================================================

df = pd.read_csv(
    INPUT_FILE,
    sep="\t",
    comment="#",
    dtype=str
)

# Clean column names
df.columns = df.columns.str.strip()

# ============================================================
# CONVERT NUMERIC COLUMNS
# ============================================================

numeric_columns = [
    "Expectation",
    "UPE",
    "miRNA_start",
    "miRNA_end",
    "Target_start",
    "Target_end"
]

for col in numeric_columns:
    df[col] = pd.to_numeric(
        df[col],
        errors="coerce"
    )

# Clean text columns
for col in df.columns:
    if df[col].dtype == "object":
        df[col] = (
            df[col]
            .fillna("")
            .astype(str)
            .str.strip()
        )

# ============================================================
# BASIC DATASET STATISTICS
# ============================================================

total_interactions = len(df)

unique_mirnas = df["miRNA_Acc."].nunique()

unique_targets = df["Target_Acc."].nunique()

# ============================================================
# 1. COMPLETE INTERACTION DATASET
# ============================================================

df.to_csv(
    "CuminGRdb_Table_S4_miRNA_Target_Interactions.csv",
    index=False
)

# ============================================================
# 2. miRNA-LEVEL SUMMARY
#
# One row = one miRNA
# ============================================================

mirna_summary = (
    df.groupby("miRNA_Acc.")
    .agg(
        Target_Genes=(
            "Target_Acc.",
            "nunique"
        ),

        Total_Interactions=(
            "Target_Acc.",
            "size"
        ),

        Mean_Expectation=(
            "Expectation",
            "mean"
        ),

        Median_Expectation=(
            "Expectation",
            "median"
        ),

        Minimum_Expectation=(
            "Expectation",
            "min"
        )
    )
    .reset_index()
)

# ------------------------------------------------------------
# Expectation thresholds
# ------------------------------------------------------------

for threshold in [1, 2, 3]:

    counts = (
        df[df["Expectation"] <= threshold]
        .groupby("miRNA_Acc.")
        .size()
    )

    mirna_summary[
        f"Interactions_LE_{threshold}"
    ] = (
        mirna_summary["miRNA_Acc."]
        .map(counts)
        .fillna(0)
        .astype(int)
    )

# ------------------------------------------------------------
# Inhibition mechanism
# ------------------------------------------------------------

cleavage = (
    df[
        df["Inhibition"]
        .str.lower()
        .eq("cleavage")
    ]
    .groupby("miRNA_Acc.")
    .size()
)

translation = (
    df[
        df["Inhibition"]
        .str.lower()
        .str.contains("translation")
    ]
    .groupby("miRNA_Acc.")
    .size()
)

mirna_summary["Cleavage_Interactions"] = (
    mirna_summary["miRNA_Acc."]
    .map(cleavage)
    .fillna(0)
    .astype(int)
)

mirna_summary[
    "Translation_Repression_Interactions"
] = (
    mirna_summary["miRNA_Acc."]
    .map(translation)
    .fillna(0)
    .astype(int)
)

# Sort by regulatory breadth
mirna_summary = mirna_summary.sort_values(
    [
        "Target_Genes",
        "Total_Interactions"
    ],
    ascending=False
)

mirna_summary.to_csv(
    "CuminGRdb_Table_S5_miRNA_Level_Summary.csv",
    index=False
)

# ============================================================
# 3. TARGET-GENE-LEVEL SUMMARY
#
# One row = one CcGene
# ============================================================

target_summary = (
    df.groupby("Target_Acc.")
    .agg(
        Targeting_miRNAs=(
            "miRNA_Acc.",
            "nunique"
        ),

        Total_Interactions=(
            "miRNA_Acc.",
            "size"
        ),

        Mean_Expectation=(
            "Expectation",
            "mean"
        ),

        Median_Expectation=(
            "Expectation",
            "median"
        ),

        Minimum_Expectation=(
            "Expectation",
            "min"
        )
    )
    .reset_index()
)

# ------------------------------------------------------------
# Expectation thresholds
# ------------------------------------------------------------

for threshold in [1, 2, 3]:

    counts = (
        df[df["Expectation"] <= threshold]
        .groupby("Target_Acc.")
        .size()
    )

    target_summary[
        f"Interactions_LE_{threshold}"
    ] = (
        target_summary["Target_Acc."]
        .map(counts)
        .fillna(0)
        .astype(int)
    )

# ------------------------------------------------------------
# Inhibition mechanism
# ------------------------------------------------------------

cleavage = (
    df[
        df["Inhibition"]
        .str.lower()
        .eq("cleavage")
    ]
    .groupby("Target_Acc.")
    .size()
)

translation = (
    df[
        df["Inhibition"]
        .str.lower()
        .str.contains("translation")
    ]
    .groupby("Target_Acc.")
    .size()
)

target_summary["Cleavage_Interactions"] = (
    target_summary["Target_Acc."]
    .map(cleavage)
    .fillna(0)
    .astype(int)
)

target_summary[
    "Translation_Repression_Interactions"
] = (
    target_summary["Target_Acc."]
    .map(translation)
    .fillna(0)
    .astype(int)
)

# ------------------------------------------------------------
# Add target description
# ------------------------------------------------------------

target_description = (
    df[
        [
            "Target_Acc.",
            "Target_Desc."
        ]
    ]
    .drop_duplicates(
        subset=["Target_Acc."]
    )
)

target_summary = target_summary.merge(
    target_description,
    on="Target_Acc.",
    how="left"
)

# Sort by number of targeting miRNAs
target_summary = target_summary.sort_values(
    [
        "Targeting_miRNAs",
        "Total_Interactions"
    ],
    ascending=False
)

target_summary.to_csv(
    "CuminGRdb_Table_S6_Target_Gene_Level_Summary.csv",
    index=False
)

# ============================================================
# 4. HIGH-CONFIDENCE DATASET
#
# Expectation <= 3
# ============================================================

high_confidence = df[
    df["Expectation"] <= 3
].copy()

high_confidence = high_confidence.sort_values(
    [
        "Expectation",
        "miRNA_Acc.",
        "Target_Acc."
    ]
)

high_confidence.to_csv(
    "CuminGRdb_Table_S7_High_Confidence_Expectation_LE3.csv",
    index=False
)

# ============================================================
# 5. EXPECTATION THRESHOLD SUMMARY
# ============================================================

threshold_rows = []

for threshold in [1, 2, 3, 4, 5]:

    subset = df[
        df["Expectation"] <= threshold
    ]

    n_interactions = len(subset)

    threshold_rows.append({

        "Expectation_Threshold":
            f"<= {threshold}.0",

        "Interactions":
            n_interactions,

        "Percentage_of_All_Interactions":
            round(
                n_interactions /
                total_interactions *
                100,
                2
            ),

        "Unique_miRNAs":
            subset[
                "miRNA_Acc."
            ].nunique(),

        "Unique_Target_Genes":
            subset[
                "Target_Acc."
            ].nunique()
    })

expectation_summary = pd.DataFrame(
    threshold_rows
)

expectation_summary.to_csv(
    "CuminGRdb_Table_S8_Expectation_Threshold_Summary.csv",
    index=False
)

# ============================================================
# 6. INHIBITION-MODE SUMMARY
# ============================================================

inhibition_summary = (
    df["Inhibition"]
    .value_counts()
    .rename_axis("Inhibition")
    .reset_index(
        name="Interactions"
    )
)

inhibition_summary["Percentage"] = (
    inhibition_summary["Interactions"]
    / total_interactions
    * 100
).round(2)

inhibition_summary.to_csv(
    "CuminGRdb_miRNA_Inhibition_Mode_Summary.csv",
    index=False
)

# ============================================================
# 7. OVERALL DATASET SUMMARY
# ============================================================

summary_rows = [

    [
        "Total interaction records",
        total_interactions
    ],

    [
        "Unique miRNAs",
        unique_mirnas
    ],

    [
        "Unique target genes",
        unique_targets
    ],

    [
        "Cleavage interactions",
        (
            df["Inhibition"]
            .str.lower()
            .eq("cleavage")
            .sum()
        )
    ],

    [
        "Translation repression interactions",
        (
            df["Inhibition"]
            .str.lower()
            .str.contains("translation")
            .sum()
        )
    ],

    [
        "Expectation <= 1",
        (
            df["Expectation"] <= 1
        ).sum()
    ],

    [
        "Expectation <= 2",
        (
            df["Expectation"] <= 2
        ).sum()
    ],

    [
        "Expectation <= 3",
        (
            df["Expectation"] <= 3
        ).sum()
    ]
]

dataset_summary = pd.DataFrame(
    summary_rows,
    columns=[
        "Metric",
        "Value"
    ]
)

dataset_summary.to_csv(
    "CuminGRdb_miRNA_Dataset_Summary.csv",
    index=False
)

# ============================================================
# FINAL REPORT
# ============================================================

print("\n")
print("=" * 60)
print("CuminGRdb miRNA target analysis")
print("=" * 60)

print(
    f"Total interactions : "
    f"{total_interactions:,}"
)

print(
    f"Unique miRNAs     : "
    f"{unique_mirnas:,}"
)

print(
    f"Unique targets    : "
    f"{unique_targets:,}"
)

print("\nExpectation thresholds:")

for threshold in [1, 2, 3, 4, 5]:

    n = (
        df["Expectation"] <= threshold
    ).sum()

    print(
        f"  <= {threshold}.0 : "
        f"{n:,} "
        f"({n / total_interactions * 100:.2f}%)"
    )

print("\nFiles generated:")

files = [
    "CuminGRdb_Table_S4_miRNA_Target_Interactions.csv",
    "CuminGRdb_Table_S5_miRNA_Level_Summary.csv",
    "CuminGRdb_Table_S6_Target_Gene_Level_Summary.csv",
    "CuminGRdb_Table_S7_High_Confidence_Expectation_LE3.csv",
    "CuminGRdb_Table_S8_Expectation_Threshold_Summary.csv",
    "CuminGRdb_miRNA_Inhibition_Mode_Summary.csv",
    "CuminGRdb_miRNA_Dataset_Summary.csv"
]

for f in files:
    print(" ", f)

print("=" * 60)

#!/usr/bin/env python3
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: generate_sm_figures.py
# Description: Visualization and distribution analysis of secondary metabolite pathway genes.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from collections import Counter
import os

# Fix matplotlib cache issue
os.environ['MPLCONFIGDIR'] = '/Users/akmahato/Documents/Cumin-genome/.mpl_cache'
os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)

out_dir = '/Users/akmahato/.gemini/antigravity/brain/8c33d28a-fe9a-4b91-b5cd-221ff3faf1bd'

df = pd.read_csv('Tabulated_Results/cumin_secondary_metabolites.tsv', sep='\t')

# Consolidate categories for cleaner visualization
def simplify_category(cat):
    if 'Terpene' in str(cat):
        return 'Terpene Synthase (TPS)'
    elif 'Other' in str(cat):
        return 'Other'
    return cat

df['Category_Clean'] = df['Metabolite_Category'].apply(simplify_category)

# ============================================================
# FIGURE A: Donut chart of secondary metabolite gene families
# ============================================================
cat_counts = df['Category_Clean'].value_counts()

# Merge 'Other' into TPS for visual clarity
labels = ['Cytochrome P450\n(CYPs)', 'Glycosyltransferase\n(UGTs)', 'Acyltransferase', 'Terpene Synthase\n(TPS)']
sizes = [123, 122, 38, 11]  # TPS = 7+2+1+1
colors = ['#E63946', '#457B9D', '#2A9D8F', '#E9C46A']
explode = (0.03, 0.03, 0.03, 0.06)

fig, ax = plt.subplots(figsize=(8, 8), facecolor='#0D1117')
wedges, texts, autotexts = ax.pie(
    sizes, labels=labels, autopct='%1.1f%%',
    colors=colors, explode=explode,
    startangle=140,
    pctdistance=0.78,
    wedgeprops=dict(width=0.45, edgecolor='#0D1117', linewidth=2.5),
    textprops={'color': '#C9D1D9', 'fontsize': 12, 'fontweight': 'bold'}
)
for at in autotexts:
    at.set_color('white')
    at.set_fontsize(11)
    at.set_fontweight('bold')

# Centre text
ax.text(0, 0.05, '294', fontsize=40, fontweight='bold', color='white', ha='center', va='center')
ax.text(0, -0.12, 'Genes', fontsize=14, color='#8B949E', ha='center', va='center')

ax.set_title('Secondary Metabolite Gene Families\nin Cuminum cyminum',
             fontsize=16, fontweight='bold', color='#C9D1D9', pad=20)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'SM_Figure_A_donut.png'), dpi=300, facecolor='#0D1117', bbox_inches='tight')
plt.close()
print('Figure A done.')

# ============================================================
# FIGURE B: Sub-classification bar chart
# ============================================================
# Glycosyltransferase sub-families
gt_descs = df[df['Category_Clean'] == 'Glycosyltransferase']['Description'].value_counts()
gt_families = {
    'UDP-GT': 0, 'GT family 8': 0, 'GT family 2\n(Cellulose synthase)': 0,
    'GT family 31': 0, 'GT family 1': 0, 'GT family 17': 0, 'GT family 20': 0, 'Other GT': 0
}
for desc, count in gt_descs.items():
    desc_lower = str(desc).lower()
    if 'udp' in desc_lower:
        gt_families['UDP-GT'] += count
    elif 'family 8' in desc_lower or 'glycosyltransferase 8' in desc_lower:
        gt_families['GT family 8'] += count
    elif 'family 2' in desc_lower or 'cellulose' in desc_lower:
        gt_families['GT family 2\n(Cellulose synthase)'] += count
    elif 'family 31' in desc_lower or 'glycosyltransferase 31' in desc_lower:
        gt_families['GT family 31'] += count
    elif 'family 1' in desc_lower and 'family 17' not in desc_lower:
        gt_families['GT family 1'] += count
    elif 'family 17' in desc_lower:
        gt_families['GT family 17'] += count
    elif 'family 20' in desc_lower:
        gt_families['GT family 20'] += count
    else:
        gt_families['Other GT'] += count

# Acyltransferase sub-families
at_descs = df[df['Category_Clean'] == 'Acyltransferase']['Description'].value_counts()
at_families = {
    'Glycerol-3-phosphate\nacyltransferase': 0,
    'Phospholipid:DAG\nacyltransferase': 0,
    'O-acyltransferase\n(WSD1-like)': 0,
    'Membrane-bound\nacyltransferase': 0,
    '2-oxoacid DH\nacyltransferase': 0,
    'DAG O-acyltransferase': 0,
    'Other AT': 0
}
for desc, count in at_descs.items():
    desc_lower = str(desc).lower()
    if 'glycerol-3' in desc_lower:
        at_families['Glycerol-3-phosphate\nacyltransferase'] += count
    elif 'phospholipid' in desc_lower:
        at_families['Phospholipid:DAG\nacyltransferase'] += count
    elif 'wsd1' in desc_lower or 'o-acyltransferase' in desc_lower:
        at_families['O-acyltransferase\n(WSD1-like)'] += count
    elif 'membrane' in desc_lower:
        at_families['Membrane-bound\nacyltransferase'] += count
    elif '2-oxoacid' in desc_lower:
        at_families['2-oxoacid DH\nacyltransferase'] += count
    elif 'diacylglycerol' in desc_lower:
        at_families['DAG O-acyltransferase'] += count
    else:
        at_families['Other AT'] += count

# TPS subfamilies
tps_families = {
    'General TPS': 7,
    'Monoterpene\n(Pinene synthase)': 1,
    'TPS-Glycosyl-\ntransferase': 2,
    'IPP condensation\n(Terpenoid backbone)': 1
}

fig, axes = plt.subplots(1, 3, figsize=(20, 7), facecolor='#0D1117')

# Panel 1: Glycosyltransferases
ax1 = axes[0]
gt_labels = list(gt_families.keys())
gt_vals = list(gt_families.values())
bars1 = ax1.barh(gt_labels, gt_vals, color='#457B9D', edgecolor='#0D1117', height=0.6)
ax1.set_xlabel('Number of Genes', color='#C9D1D9', fontsize=11)
ax1.set_title('Glycosyltransferase\nSub-families (122)', color='#C9D1D9', fontsize=13, fontweight='bold')
ax1.set_facecolor('#161B22')
ax1.tick_params(colors='#8B949E')
ax1.spines['bottom'].set_color('#30363D')
ax1.spines['left'].set_color('#30363D')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
for bar, val in zip(bars1, gt_vals):
    ax1.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2, str(val),
             va='center', color='#C9D1D9', fontsize=10, fontweight='bold')
ax1.invert_yaxis()

# Panel 2: Acyltransferases
ax2 = axes[1]
at_labels = list(at_families.keys())
at_vals = list(at_families.values())
bars2 = ax2.barh(at_labels, at_vals, color='#2A9D8F', edgecolor='#0D1117', height=0.6)
ax2.set_xlabel('Number of Genes', color='#C9D1D9', fontsize=11)
ax2.set_title('Acyltransferase\nSub-families (38)', color='#C9D1D9', fontsize=13, fontweight='bold')
ax2.set_facecolor('#161B22')
ax2.tick_params(colors='#8B949E')
ax2.spines['bottom'].set_color('#30363D')
ax2.spines['left'].set_color('#30363D')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
for bar, val in zip(bars2, at_vals):
    ax2.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2, str(val),
             va='center', color='#C9D1D9', fontsize=10, fontweight='bold')
ax2.invert_yaxis()

# Panel 3: Terpene Synthases
ax3 = axes[2]
tps_labels = list(tps_families.keys())
tps_vals = list(tps_families.values())
bars3 = ax3.barh(tps_labels, tps_vals, color='#E9C46A', edgecolor='#0D1117', height=0.6)
ax3.set_xlabel('Number of Genes', color='#C9D1D9', fontsize=11)
ax3.set_title('Terpene Synthase\nSub-families (11)', color='#C9D1D9', fontsize=13, fontweight='bold')
ax3.set_facecolor('#161B22')
ax3.tick_params(colors='#8B949E')
ax3.spines['bottom'].set_color('#30363D')
ax3.spines['left'].set_color('#30363D')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
for bar, val in zip(bars3, tps_vals):
    ax3.text(bar.get_width() + 0.15, bar.get_y() + bar.get_height()/2, str(val),
             va='center', color='#C9D1D9', fontsize=10, fontweight='bold')
ax3.invert_yaxis()

fig.suptitle('Sub-classification of Secondary Metabolite Gene Families in C. cyminum',
             fontsize=15, fontweight='bold', color='#C9D1D9', y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'SM_Figure_B_subfamilies.png'), dpi=300, facecolor='#0D1117', bbox_inches='tight')
plt.close()
print('Figure B done.')

# ============================================================
# FIGURE C: KEGG Pathway mapping for secondary metabolite genes
# ============================================================
# Curated KEGG pathway names relevant to secondary metabolism
kegg_map = {
    'ko00904': 'Diterpenoid biosynthesis',
    'ko00905': 'Brassinosteroid biosynthesis',
    'ko00908': 'Zeatin biosynthesis',
    'ko00140': 'Steroid hormone biosynthesis',
    'ko00380': 'Tryptophan metabolism',
    'ko00830': 'Retinol metabolism',
    'ko00980': 'Xenobiotics by CYP',
    'ko00561': 'Glycerolipid metabolism',
    'ko00564': 'Glycerophospholipid metabolism',
    'ko00590': 'Arachidonic acid metabolism',
    'ko00591': 'Linoleic acid metabolism',
}

kegg_genes = df[df['KEGG_Pathway'].notna() & (df['KEGG_Pathway'] != '-')]
pathway_category_counts = {}
for _, row in kegg_genes.iterrows():
    cat = simplify_category(row['Metabolite_Category'])
    pathways = str(row['KEGG_Pathway']).split(',')
    for p in pathways:
        p = p.strip()
        if p in kegg_map:
            key = kegg_map[p]
            if key not in pathway_category_counts:
                pathway_category_counts[key] = Counter()
            pathway_category_counts[key][cat] += 1

# Sort pathways by total count
pathway_totals = {k: sum(v.values()) for k, v in pathway_category_counts.items()}
sorted_pathways = sorted(pathway_totals, key=pathway_totals.get, reverse=True)

categories = ['Cytochrome P450', 'Glycosyltransferase', 'Acyltransferase', 'Terpene Synthase (TPS)']
cat_colors = {'Cytochrome P450': '#E63946', 'Glycosyltransferase': '#457B9D',
              'Acyltransferase': '#2A9D8F', 'Terpene Synthase (TPS)': '#E9C46A'}

fig, ax = plt.subplots(figsize=(14, 7), facecolor='#0D1117')
x = np.arange(len(sorted_pathways))
width = 0.2
for i, cat in enumerate(categories):
    vals = [pathway_category_counts[p].get(cat, 0) for p in sorted_pathways]
    ax.bar(x + i * width, vals, width, label=cat, color=cat_colors[cat], edgecolor='#0D1117')

ax.set_xticks(x + width * 1.5)
ax.set_xticklabels(sorted_pathways, rotation=40, ha='right', color='#C9D1D9', fontsize=10)
ax.set_ylabel('Number of Genes', color='#C9D1D9', fontsize=12)
ax.set_title('KEGG Pathway Distribution of Secondary Metabolite Genes',
             color='#C9D1D9', fontsize=14, fontweight='bold', pad=15)
ax.set_facecolor('#161B22')
ax.tick_params(colors='#8B949E')
ax.spines['bottom'].set_color('#30363D')
ax.spines['left'].set_color('#30363D')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(facecolor='#161B22', edgecolor='#30363D', labelcolor='#C9D1D9', fontsize=10, loc='upper right')

fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'SM_Figure_C_kegg.png'), dpi=300, facecolor='#0D1117', bbox_inches='tight')
plt.close()
print('Figure C done.')

# ============================================================
# FIGURE D: Functional annotation coverage (stacked bar)
# ============================================================
annotation_types = ['SwissProt Hit', 'GO Terms', 'KEGG Pathway']
annotated = []
not_annotated = []

for cat in ['Cytochrome P450', 'Glycosyltransferase', 'Acyltransferase', 'Terpene Synthase (TPS)']:
    if cat == 'Terpene Synthase (TPS)':
        sub = df[df['Category_Clean'] == cat]
    else:
        sub = df[df['Category_Clean'] == cat]
    total = len(sub)
    sp = len(sub[(sub['SwissProt_Hit'].notna()) & (sub['SwissProt_Hit'] != '-')])
    go = len(sub[(sub['GOs'].notna()) & (sub['GOs'] != '-')])
    kegg = len(sub[(sub['KEGG_Pathway'].notna()) & (sub['KEGG_Pathway'] != '-')])
    print(f'{cat}: Total={total}, SwissProt={sp}, GO={go}, KEGG={kegg}')

fig, ax = plt.subplots(figsize=(10, 6), facecolor='#0D1117')
cats_display = ['CYPs\n(123)', 'UGTs\n(122)', 'Acyltransferase\n(38)', 'TPS\n(11)']
sp_vals = [120, 122, 38, 9]
go_vals = [117, 41, 29, 7]
kegg_vals = [110, 25, 24, 8]

x = np.arange(len(cats_display))
width = 0.22

bars1 = ax.bar(x - width, sp_vals, width, label='SwissProt Hit', color='#E63946', edgecolor='#0D1117')
bars2 = ax.bar(x, go_vals, width, label='GO Terms', color='#457B9D', edgecolor='#0D1117')
bars3 = ax.bar(x + width, kegg_vals, width, label='KEGG Pathway', color='#2A9D8F', edgecolor='#0D1117')

for bars in [bars1, bars2, bars3]:
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 1, str(int(height)),
                ha='center', va='bottom', color='#C9D1D9', fontsize=10, fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(cats_display, color='#C9D1D9', fontsize=11)
ax.set_ylabel('Number of Genes Annotated', color='#C9D1D9', fontsize=12)
ax.set_title('Functional Annotation Coverage of Secondary Metabolite Genes',
             color='#C9D1D9', fontsize=14, fontweight='bold', pad=15)
ax.set_facecolor('#161B22')
ax.tick_params(colors='#8B949E')
ax.spines['bottom'].set_color('#30363D')
ax.spines['left'].set_color('#30363D')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(facecolor='#161B22', edgecolor='#30363D', labelcolor='#C9D1D9', fontsize=11)

fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'SM_Figure_D_annotation.png'), dpi=300, facecolor='#0D1117', bbox_inches='tight')
plt.close()
print('Figure D done.')

print('\nAll figures generated successfully!')

#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_getorganelle.sh
# Description: Organelle genome assembly (Chloroplast and Mitochondria) using GetOrganelle.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

# Step 3b: Organelle Assembly
# GetOrganelle

# Chloroplast and Mitochondria
get_organelle_from_reads.py -1 trimmed_reads/sample_R1_paired.fq.gz -2 trimmed_reads/sample_R2_paired.fq.gz \
    -o chloroplast_out -R 15 -k 21,45,65,85,105 -F embplant_pt --only-assembler

get_organelle_from_reads.py -1 trimmed_reads/sample_R1_paired.fq.gz -2 trimmed_reads/sample_R2_paired.fq.gz \
    -o mito_out -R 50 -k 21,45,65,85,105 -F embplant_mt --only-assembler

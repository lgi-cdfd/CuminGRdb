#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_masurca.sh
# Description: De novo hybrid genome assembly using MaSuRCA assembler.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

# Step 3: Genome Assembly
# MaSuRCA

# Generate configuration file for MaSuRCA
cat <<EOF > config.txt
DATA
PE= pe 350 50 trimmed_reads/sample_R1_paired.fq.gz trimmed_reads/sample_R2_paired.fq.gz
END
PARAMETERS
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 1
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.15
NUM_THREADS = 64
JF_SIZE = 13000000000
DO_HOMOPOLYMER_TRIM = 0
LHE_COVERAGE = 25
EOF

# Run MaSuRCA
masurca config.txt
./assemble.sh

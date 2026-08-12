#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_edta.sh
# Description: Transposable element and repeat annotation using EDTA (Extensive de-novo TE Annotator).
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

# Step 4a: Repeat Annotation
# EDTA

perl EDTA.pl --genome cumin_assembly.fasta --species others --step all --cds cumin_cds.fa --threads 64

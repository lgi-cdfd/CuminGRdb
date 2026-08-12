#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_functional_annotation.sh
# Description: Functional gene annotation and ontology mapping using OmicsBox / Blast2GO.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

# Step 4c: Functional Annotation
# Diamond (NR & Swiss-Prot) and EggNOG

# Diamond against NR
diamond blastp -q cumin_proteins.fa -d nr.dmnd -e 1e-5 --sensitive --max-target-seqs 1 --threads 64 --outfmt 6 -o diamond_nr.tsv

# Diamond against Swiss-Prot
diamond blastp -q cumin_proteins.fa -d swissprot.dmnd -e 1e-5 --sensitive --max-target-seqs 1 --threads 64 --outfmt 6 -o diamond_swissprot.tsv

# EggNOG-mapper
emapper.py -i cumin_proteins.fa --output cumin_eggnog -m diamond --cpu 64

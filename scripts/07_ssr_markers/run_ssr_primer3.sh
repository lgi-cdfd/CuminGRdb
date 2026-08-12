#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_krait_ssr.sh
# Description: Discovery and PCR primer design for Simple Sequence Repeat (SSR) markers using Krait tool.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

# Step 5b: SSR Mining & Primer Design
# Krait: An ultrafast tool for genome-wide SSR mining and primer design (v1.3+)

# Run Krait engine on reference genome assembly (Cuminum_cyminum_genome_assembly.fasta)
# Mined 213,248 microsatellite loci (mono-, di-, tri-, tetra-, penta-, hexa-repeats)
# Designed Primer3 PCR primers (forward/reverse primer sequences, Tm 57-63°C, product sizes 100-300 bp)
# Output tables: ssr_markers_primers.csv & CuminDB_SSR_Summary.csv

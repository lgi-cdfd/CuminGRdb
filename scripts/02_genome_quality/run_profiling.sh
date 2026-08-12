#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_profiling.sh
# Description: Genome profiling and k-mer distribution analysis.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

# Step 2: Genome Profiling
# Jellyfish, GenomeScope2, Smudgeplot

# 1. Count 21-mers
jellyfish count -C -m 21 -s 10G -t 32 -o cumin.jf trimmed_reads/*_paired.fq.gz

# 2. Export k-mer histogram
jellyfish histo -t 32 cumin.jf > cumin.histo

# 3. Run GenomeScope2
genomescope2 -i cumin.histo -o genomescope_out -k 21 -p 2

# 4. Run Smudgeplot (for ploidy)
L=$(smudgeplot.py cutoff cumin.histo L)
U=$(smudgeplot.py cutoff cumin.histo U)
kmc -k21 -t32 -m64 -ci1 -cs10000 @trimmed_reads.lst kmc_db tmp
kmc_tools transform kmc_db histogram kmc_db.hist -cx10000
smudgeplot.py hetkmers -o smudgeplot_out < kmc_db.hist
smudgeplot.py plot -o smudgeplot smudgeplot_out_coverages.tsv

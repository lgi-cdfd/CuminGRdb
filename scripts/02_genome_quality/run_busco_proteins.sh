#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_busco_proteins.sh
# Description: Genome assembly completeness evaluation using BUSCO plant lineage datasets.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

set -e

WORKDIR="/home/ajay/cumin-genome/Cumin_GigaScience/05_Gene_Annotation"
cd $WORKDIR

source ~/miniforge3/bin/activate busco6

PROTEINS="cumin_evm_proteins.fasta"

if [ ! -f "$PROTEINS" ]; then
    echo "ERROR: Protein FASTA ($PROTEINS) not found."
    exit 1
fi

echo "Running BUSCO on EVM predicted proteins..."
busco -i $PROTEINS -l embryophyta_odb10 -o busco_protein_evaluation -m proteins -c 60 -f

echo "BUSCO analysis on proteins complete!"

#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_evm.sh
# Description: Structural gene annotation and consensus model integration using EvidenceModeler (EVM).
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

# Step 4b: Structural Annotation
# EvidenceModeler (EVM)

# Create weights file
cat <<EOF > weights.txt
ABINITIO_PREDICTION	GeneMark	1
ABINITIO_PREDICTION	SNAP	1
PROTEIN	miniprot	5
TRANSCRIPT	transdecoder	10
EOF

evidence_modeler.pl --genome cumin_masked.fasta \
    --weights weights.txt \
    --gene_predictions abinitio.gff3 \
    --protein_alignments miniprot.gff3 \
    --transcript_alignments transcripts.gff3 \
    --segmentSize 100000 --overlapSize 10000 > evm.out

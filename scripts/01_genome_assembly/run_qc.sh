#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: run_qc.sh
# Description: Genome assembly and initial quality control (FastQC & Trimmomatic).
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

# Step 1: Quality Control
# FastQC and Trimmomatic

mkdir -p qc_results trimmed_reads

# Run FastQC on raw reads
fastqc raw_reads/*_R1.fastq.gz raw_reads/*_R2.fastq.gz -o qc_results/ -t 16

# Run Trimmomatic
java -jar trimmomatic-0.39.jar PE -threads 32 \
    raw_reads/sample_R1.fastq.gz raw_reads/sample_R2.fastq.gz \
    trimmed_reads/sample_R1_paired.fq.gz trimmed_reads/sample_R1_unpaired.fq.gz \
    trimmed_reads/sample_R2_paired.fq.gz trimmed_reads/sample_R2_unpaired.fq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Run MultiQC
multiqc qc_results/ -o qc_results/

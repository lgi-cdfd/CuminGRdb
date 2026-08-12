#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: setup_secondary_metabolism.sh
# Description: Mining and profile identification of secondary metabolite biosynthetic genes.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

set -e

WORKDIR="/home/ajay/cumin-genome/Cumin_GigaScience/09_Secondary_Metabolism"
mkdir -p $WORKDIR/HMM_profiles
cd $WORKDIR

echo "Setting up HMMER environment..."
source ~/miniforge3/bin/activate base
if ! conda env list | grep -q "hmmer_env"; then
    mamba create -n hmmer_env -c bioconda hmmer -y
fi

echo "Creating execution script..."
cat << 'EOF' > run_hmmsearch.sh
#!/bin/bash
set -e
source ~/miniforge3/bin/activate hmmer_env
PROTEINS="/home/ajay/cumin-genome/Cumin_GigaScience/05_Gene_Annotation/cumin_evm_proteins.fasta"

mkdir -p HMM_profiles
cd HMM_profiles
if [ ! -f "Pfam-A.hmm" ]; then
    echo "Downloading full Pfam-A database (this ensures we have all profiles)..."
    wget -qO Pfam-A.hmm.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    echo "Indexing Pfam-A..."
    hmmpress Pfam-A.hmm
    
    echo "Extracting specific profiles for Terpenes and Cytochrome P450s..."
    hmmfetch Pfam-A.hmm PF01397 > PF01397_Terpene_synthase_N.hmm
    hmmfetch Pfam-A.hmm PF03936 > PF03936_Terpene_synthase_C.hmm
    hmmfetch Pfam-A.hmm PF00067 > PF00067_p450.hmm
fi

cd ..
mkdir -p results
echo "Running hmmsearch against Cumin genome..."
hmmsearch --tblout results/TPS_N.tblout -E 1e-5 HMM_profiles/PF01397_Terpene_synthase_N.hmm $PROTEINS > /dev/null
hmmsearch --tblout results/TPS_C.tblout -E 1e-5 HMM_profiles/PF03936_Terpene_synthase_C.hmm $PROTEINS > /dev/null
hmmsearch --tblout results/CYP450.tblout -E 1e-5 HMM_profiles/PF00067_p450.hmm $PROTEINS > /dev/null

echo "Filtering hits..."
grep -v "^#" results/TPS_N.tblout | awk '{print $1}' | sort -u > results/TPS_N_genes.txt
grep -v "^#" results/TPS_C.tblout | awk '{print $1}' | sort -u > results/TPS_C_genes.txt
grep -v "^#" results/CYP450.tblout | awk '{print $1}' | sort -u > results/CYP450_genes.txt

echo "Secondary Metabolism HMM search completed successfully!"
EOF
chmod +x run_hmmsearch.sh

echo "Setup complete. Secondary metabolism pipeline scaffolded."

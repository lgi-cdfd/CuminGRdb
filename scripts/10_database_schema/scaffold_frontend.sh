#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: scaffold_frontend.sh
# Description: Web presentation layer setup and component scaffolding for CuminGRdb.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

set -e

source ~/miniforge3/bin/activate base
if ! conda env list | grep -q "web_env"; then
    echo "Creating Web Frontend conda environment..."
    mamba create -n web_env -c conda-forge nodejs=20 npm -y
fi
source ~/miniforge3/bin/activate web_env

FRONTEND_DIR="/home/ajay/cumin-genome/Cumin_GigaScience/10_CuminBase/frontend"
mkdir -p $FRONTEND_DIR
cd $FRONTEND_DIR

# Check options first
npx -y create-vite@latest --help

# Create the React project in the current directory, using the react-ts template non-interactively
if [ ! -f "package.json" ]; then
    echo "Scaffolding React + Vite app..."
    npx -y create-vite@latest ./ --template react-ts
    
    echo "Installing core dependencies..."
    # Axios for API calls, React Router for navigation, TailwindCSS for styling (Wait, guidelines say: "Avoid using TailwindCSS unless the USER explicitly requests it")
    # I will use Vanilla CSS for maximum flexibility!
    npm install axios react-router-dom lucide-react
fi

echo "Frontend scaffolding complete!"

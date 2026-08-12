#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: setup_cuminbase.sh
# Description: Database schema setup and FastAPI backend initialization for CuminGRdb.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

set -e

WORKDIR="/home/ajay/cumin-genome/Cumin_GigaScience/10_CuminBase"
cd $WORKDIR

echo "Scaffolding FastAPI Backend..."
mkdir -p backend/app/{models,routes,services,schemas,core}
touch backend/app/__init__.py

cat << 'EOF' > backend/app/main.py
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(title="CuminBase API", version="1.0.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def read_root():
    return {"message": "Welcome to CuminBase API (GigaScience Release)"}

@app.get("/api/v1/genes/{gene_id}")
def get_gene(gene_id: str):
    return {"gene_id": gene_id, "status": "Not implemented yet"}
EOF

cat << 'EOF' > backend/requirements.txt
fastapi==0.104.1
uvicorn==0.24.0
sqlalchemy==2.0.23
psycopg2-binary==2.9.9
pydantic==2.5.2
EOF

echo "Scaffolding React Frontend..."
mkdir -p frontend/src/{components,pages,services,utils,assets}

cat << 'EOF' > frontend/package.json
{
  "name": "cuminbase-frontend",
  "version": "1.0.0",
  "private": true,
  "dependencies": {
    "react": "^18.2.0",
    "react-dom": "^18.2.0",
    "react-router-dom": "^6.20.0",
    "axios": "^1.6.2",
    "@jbrowse/react-linear-genome-view": "^2.10.1"
  }
}
EOF

echo "CuminBase structure scaffolded successfully."

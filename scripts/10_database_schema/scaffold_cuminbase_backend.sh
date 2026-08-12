#!/bin/bash
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: scaffold_cuminbase_backend.sh
# Description: Backend ORM models and RESTful API route scaffolding for CuminGRdb.
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================

set -e

WORKDIR="/home/ajay/cumin-genome/Cumin_GigaScience/10_CuminBase/backend"
cd $WORKDIR

mkdir -p app/core app/models

# Core Database Setup
cat << 'EOF' > app/core/database.py
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base, sessionmaker

SQLALCHEMY_DATABASE_URL = "sqlite:///./cuminbase.db"

engine = create_engine(
    SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False}
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
EOF

# Gene Models
cat << 'EOF' > app/models/gene.py
from sqlalchemy import Column, String, Integer, ForeignKey, Text
from sqlalchemy.orm import relationship
from app.core.database import Base

class Gene(Base):
    __tablename__ = "genes"
    id = Column(String, primary_key=True, index=True)
    chromosome = Column(String, index=True)
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(String)

    transcripts = relationship("Transcript", back_populates="gene")
    ssr_markers = relationship("SSR", back_populates="overlapping_gene")

class Transcript(Base):
    __tablename__ = "transcripts"
    id = Column(String, primary_key=True, index=True)
    gene_id = Column(String, ForeignKey("genes.id"))
    
    gene = relationship("Gene", back_populates="transcripts")
    annotation = relationship("Annotation", back_populates="transcript", uselist=False)
    tf_family = relationship("TranscriptionFactor", back_populates="transcript", uselist=False)

class Annotation(Base):
    __tablename__ = "annotations"
    transcript_id = Column(String, ForeignKey("transcripts.id"), primary_key=True)
    go_terms = Column(Text)
    interpro_domains = Column(Text)
    kegg_pathways = Column(Text)
    ec_number = Column(String)
    description = Column(Text)

    transcript = relationship("Transcript", back_populates="annotation")

class TranscriptionFactor(Base):
    __tablename__ = "transcription_factors"
    transcript_id = Column(String, ForeignKey("transcripts.id"), primary_key=True)
    family = Column(String, index=True)

    transcript = relationship("Transcript", back_populates="tf_family")
EOF

# SSR Models
cat << 'EOF' > app/models/ssr.py
from sqlalchemy import Column, String, Integer, ForeignKey, Text
from sqlalchemy.orm import relationship
from app.core.database import Base

class SSR(Base):
    __tablename__ = "ssrs"
    marker_id = Column(String, primary_key=True, index=True)
    chromosome = Column(String, index=True)
    start = Column(Integer)
    end = Column(Integer)
    motif = Column(String)
    repeats = Column(Integer)
    overlapping_gene_id = Column(String, ForeignKey("genes.id"), nullable=True)

    overlapping_gene = relationship("Gene", back_populates="ssr_markers")
    primers = relationship("SSRPrimer", back_populates="ssr", uselist=False)

class SSRPrimer(Base):
    __tablename__ = "ssr_primers"
    marker_id = Column(String, ForeignKey("ssrs.marker_id"), primary_key=True)
    forward_sequence = Column(String)
    forward_tm = Column(String)
    reverse_sequence = Column(String)
    reverse_tm = Column(String)
    product_size = Column(String)

    ssr = relationship("SSR", back_populates="primers")
EOF

# Include models in __init__
cat << 'EOF' > app/models/__init__.py
from app.models.gene import Gene, Transcript, Annotation, TranscriptionFactor
from app.models.ssr import SSR, SSRPrimer
from app.core.database import Base
EOF

echo "SQLAlchemy models created."

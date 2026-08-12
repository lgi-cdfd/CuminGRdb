# CuminGRdb: Cuminum cyminum L. Genome Project

[![Data License: CC BY 4.0](https://img.shields.io/badge/Data_License-CC_BY_4.0-blue.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Code License: MIT](https://img.shields.io/badge/Code_License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![NCBI BioProject](https://img.shields.io/badge/NCBI_BioProject-PRJNA1469524-007ec6.svg)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1469524/)
[![NCBI SRA](https://img.shields.io/badge/NCBI_SRA-SRR10076761-059669.svg)](https://www.ncbi.nlm.nih.gov/sra/SRR10076761)
[![Python](https://img.shields.io/badge/Python-v3.10%2B-3776AB.svg)](https://www.python.org/)
[![Database Portal](https://img.shields.io/badge/Database_Portal-CuminGRdb-0284c7.svg)](https://database.cdfd.org.in/cumingrdb/)

This repository contains the complete bioinformatics processing scripts, data models, and analytical pipelines used for the draft genome assembly, structural/functional annotation, and downstream multi-omics analysis of cumin (_Cuminum cyminum_ L.), as presented in **CuminGRdb**.

---

## 👥 Authors & Affiliations

**Ajay Kumar Mahato**<sup>1*</sup>, **Ramesh Eerapagula**<sup>1</sup>, **Rakesh Singh**<sup>2</sup>, **Avinash Mishra**<sup>3</sup>, **Lakshmi Devi**<sup>1</sup>, **Priyanka Kushwaha**<sup>1</sup>, **Ankit Bhagat**<sup>1</sup>, **Bishun Deo Prasad**<sup>4</sup>, **Sangita Sahni**<sup>5</sup>, **Nagendra Kumar Singh**<sup>6</sup>

1. **Laboratory of Genome Informatics**, Centre for DNA Fingerprinting and Diagnostics (CDFD), Hyderabad, Telangana, India
2. **Division of Genomic Resources**, ICAR–National Bureau of Plant Genetic Resources (ICAR-NBPGR), Pusa Campus, New Delhi 110012, India
3. **Applied Phycology and Biotechnology Division**, CSIR–Central Salt & Marine Chemicals Research Institute (CSIR-CSMCRI), Bhavnagar, Gujarat, India
4. **Department of Agricultural Biotechnology and Molecular Biology**, College of Basic Sciences and Humanities, Dr. Rajendra Prasad Central Agricultural University, Pusa, Samastipur, Bihar 848125, India
5. **Department of Plant Pathology**, T.C.A., Dholi, Dr. Rajendra Prasad Central Agricultural University, Pusa, Bihar, India
6. **National Institute for Plant Biotechnology** (ICAR-NIPB), Pusa Campus, New Delhi 110012, India

📧 **Corresponding Author**: Dr. Ajay Kumar Mahato (`akmahato@cdfd.org.in`)  
🌐 **Web Portal**: [https://database.cdfd.org.in/cumingrdb/](https://database.cdfd.org.in/cumingrdb/)

---

## 📦 Reference Datasets & Download Links (`data/`)

The compressed FASTA sequence reference files are provided in `data/`:

| Dataset File | Format | Size | Download / Repository Access |
| :--- | :---: | :---: | :--- |
| **Genome Assembly FASTA** (147,524 contigs) | `.fasta` | 1.25 GB | [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1469524/) \| [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra/SRR10076761) \| [ENA Browser](https://www.ebi.ac.uk/ena/browser/view/SRR10076761) \| [Portal Direct Download](https://database.cdfd.org.in/cumingrdb/api/download/assembly) |
| **Predicted CDS Sequences** | `.fasta.gz` | 11 MB | Included in `data/raw/Cuminum_cyminum_predicted_cds.fasta.gz` |
| **Predicted Proteome Sequences** | `.fasta.gz` | 6.5 MB | Included in `data/raw/Cuminum_cyminum_predicted_protein.fasta.gz` |
| **Chloroplast Genome** | `.fasta.gz` | 73 KB | Included in `data/organelle/Cuminum_cyminum_chloroplast.fasta.gz` |
| **Mitochondria Genome** | `.fasta.gz` | 2.9 MB | Included in `data/organelle/Cuminum_cyminum_mitochondria.fasta.gz` |

---

## 📊 Standardized Results Datasets (`results/tables/`)

All processed result files follow standardized publication-grade `CuminGRdb_` naming conventions:

- **`CuminGRdb_SSR_Markers_and_PCR_Primers.csv`** (66 MB, `.gz` 6.4 MB): 213,248 microsatellite loci mined via **Krait** with Primer3 PCR primers.
- **`CuminGRdb_SSR_Distribution_Summary.csv`**: Summary of SSR motif classes (mono-, di-, tri-, tetra-, penta-, hexa-repeats).
- **`CuminGRdb_PlantTFDB_Predictions.txt`**: 1,362 transcription factor predictions from **PlantTFDB v5.0**.
- **`CuminGRdb_TF_Family_Distribution.csv`**: Family breakdown across 50 plant TF families.
- **`CuminGRdb_miRNA_Target_Interactions.csv`** (88 MB, `.gz` 14 MB): 887,911 **psRNATarget** microRNA target interaction records.
- **`CuminGRdb_miRNA_High_Confidence_Expectation_LE3.csv`** (5.7 MB): 57,454 high-confidence target interactions ($E \le 3.0$).
- **`CuminGRdb_miRNA_Level_Summary.csv`**: 8,185 unique miRNAs summary statistics.
- **`CuminGRdb_Target_Gene_Level_Summary.csv`**: 30,181 target genes summary statistics.
- **`CuminGRdb_Functional_Annotations.csv`** & **`CuminGRdb_Combined_Functional_Annotation.tsv`**: 33,595 OmicsBox / Blast2GO functional gene annotations.
- **`CuminGRdb_Secondary_Metabolite_Genes.csv`** & **`CuminGRdb_Secondary_Metabolite_Pathways.tsv`**: Secondary metabolite pathway gene catalog.
- **`CuminGRdb_TE_and_Repeats_Summary.csv`**: Transposable element & repeat annotation breakdown (occupying 99.44% genome coverage).
- **`BUSCO_Embryophyta_Full_Table.tsv`** & **`BUSCO_Embryophyta_Missing_Genes.tsv`**: Assembly completeness scores.

---

## 📂 Complete Pipeline Scripts Inventory (`scripts/`)

```
ScientificData_Project-Github/scripts/
├── 01_genome_assembly/                       (FastQC, Trimmomatic, MaSuRCA, GetOrganelle)
├── 02_genome_quality/                        (GenomeScope profiling & BUSCO completeness)
├── 03_repeat_annotation/                     (EDTA transposable element annotation)
├── 04_structural_annotation/                 (EvidenceModeler EVM & PASA gene modeling)
├── 05_functional_annotation/                 (Diamond NR/Swiss-Prot & OmicsBox annotations)
├── 06_transcription_factors/                 (PlantTFDB v5.0 web-tool TF classification)
├── 07_ssr_markers/                           (Krait ultrafast SSR mining & Primer3 primers)
├── 08_mirna_targets/                         (psRNATarget web-tool miRNA target prediction)
├── 09_secondary_metabolites/                 (HMMER profile search & pathway visualization)
└── 10_database_schema/                       (FastAPI/Express schema & frontend scaffolding)
```

---

## 🚀 Bioinformatic Script Header Standard

All scripts inside `scripts/` follow a uniform, standardized header block format:

```bash
#!/bin/bash (or #!/usr/bin/env python3)
# ==============================================================================
# Project: CuminGRdb (Cuminum cyminum L. Genome Project)
# Script: <script_name>
# Description: <description>
# Author: Dr. Ajay Kumar Mahato
# Affiliation: Laboratory of Genome Informatics, BRIC-CDFD, India
# Repository: https://database.cdfd.org.in/cumingrdb/
# License: MIT License
# ==============================================================================
```

---

## 📖 Citation

If you use these scripts, datasets, or the CuminGRdb portal in your research, please cite:

```bibtex
@article{Mahato2026CuminGRdb,
  author = {Mahato, Ajay Kumar and Eerapagula, Ramesh and Singh, Rakesh and Mishra, Avinash and Devi, Lakshmi and Kushwaha, Priyanka and Bhagat, Ankit and Prasad, Bishun Deo and Sahni, Sangita and Singh, Nagendra Kumar},
  title = {A high-quality draft genome assembly, comprehensive annotation, and integrated breeding database (CuminGRdb) for cumin (Cuminum cyminum L.)},
  journal = {Scientific Data},
  year = {2026},
  publisher = {Nature Publishing Group}
}
```

---

## 📜 License

- **Source Code & Pipeline Scripts**: [MIT License](https://opensource.org/licenses/MIT)
- **Genomic & Functional Datasets**: Creative Commons Attribution 4.0 International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/))

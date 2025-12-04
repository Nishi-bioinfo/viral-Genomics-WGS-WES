# viral-Genomics-WGS-WES

WGS sequences the entire viral genome to identify mutations, track variants, and understand viral evolution. It provides a complete view of the virus for surveillance, outbreak investigation, and vaccine/drug development. WES focuses on protein-coding regions to detect changes that may affect viral infectivity, immune escape, or drug resistance.


This repository contains a small example pipeline and project layout for performing viral whole-genome (WGS) / whole-exome (WES) short-read analysis: download raw reads, run QC, align to a reference, and call variants. The sample dataset used in this repository is `ERR5743893` and the reference is `MN908947` (SARS-CoV-2).

**Goals**
- Provide an auditable set of commands (see `commands.txt`) used to produce the processed files in this repo.
- Show a minimal folder layout for alignment, QC reports and variant calls.
- Offer reproducibility guidance (packages and an example `conda` environment).

**Important:** This repo stores example outputs (BAM, VCF, QC HTML). If you plan to re-run the pipeline, follow the commands in `commands.txt` and the `Quickstart` below.

**Repository Structure**
- `Raw_Data/`: Raw input FASTQ files (example: `ERR5743893_1.fastq.gz`, `ERR5743893_2.fastq.gz`).
- `QC_Reports/`: FastQC HTML reports (example: `ERR5743893_1_fastqc.html`).
- `Alignment/`: Alignment outputs (SAM/BAM/BAI). Example: `ERR5743893.sorted.bam`, `ERR5743893.sorted.bam.bai`.
- `Reference/`: Reference FASTA and index files (example: `MN908947.fasta` and its index files).
- `Variant_Calling/`: Variant call files (VCF), e.g. `ERR5743893.vcf`.
- `commands.txt`: The shell command history used to generate the files in this repository (canonical commands are copied below).

**Prerequisites**
- Operating system: Linux (commands shown use bash).
- Tools used (examples): `wget`, `fastqc`, `bwa`, `samtools`, `freebayes` (or other variant callers), `bcftools`.
- Preferred install method: `conda` (Bioconda) or system packages (apt for quick installs). Using conda is recommended to pin versions and avoid permission issues.

Quicknote: the example `commands.txt` uses both `apt` and `conda` lines; prefer one package manager for a reproducible environment.

**Quickstart (example commands)**
The commands below were executed to produce the example outputs and are taken from `commands.txt`. Run them from the repository root (adjust thread counts and paths to your environment):

```bash
# Create folders
mkdir -p Raw_Data QC_Reports Alignment Variant_Calling Reference

# Download example reads
cd Raw_Data
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_2.fastq.gz

# Run FastQC (outputs to ../QC_Reports)
fastqc ERR5743893_1.fastq.gz --outdir ../QC_Reports
fastqc ERR5743893_2.fastq.gz --outdir ../QC_Reports

# (optional) decompress if your pipeline needs uncompressed FASTQ
gunzip ERR5743893_1.fastq.gz
gunzip ERR5743893_2.fastq.gz
cd ..

# Download reference (MN908947) into Reference/
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" -O Reference/MN908947.fasta

# Index reference
bwa index Reference/MN908947.fasta
samtools faidx Reference/MN908947.fasta

# Align reads, convert, sort and index
bwa mem Reference/MN908947.fasta Raw_Data/ERR5743893_1.fastq Raw_Data/ERR5743893_2.fastq > Alignment/ERR5743893.sam
samtools view -@ 20 -S -b Alignment/ERR5743893.sam > Alignment/ERR5743893.bam
samtools sort -@ 4 -o Alignment/ERR5743893.sorted.bam Alignment/ERR5743893.bam
samtools index Alignment/ERR5743893.sorted.bam

# Call variants with freebayes (example)
freebayes -f Reference/MN908947.fasta Alignment/ERR5743893.sorted.bam > Variant_Calling/ERR5743893.vcf
```

**Reproducibility / Suggested conda environment**
Create an isolated environment with pinned versions (example):

```yaml
name: viral-wgs-env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - bwa=0.7.17
  - samtools=1.17
  - fastqc=0.11.9
  - freebayes=1.3.5
  - bcftools=1.17
  - wget
  - gzip
```

Create and activate:
```bash
conda env create -f environment.yml   # or `conda create -n viral-wgs-env ...`
conda activate viral-wgs-env
```

**Notes & suggestions**
- Choose a variant caller appropriate for viral data and low-frequency variants (FreeBayes, LoFreq, iVar, etc.).
- For paired-end trimming and filtering consider `fastp` or `trimmomatic` before alignment.
- When working with many samples, wrap this workflow in a workflow manager (Snakemake, Nextflow) and add containerization (Docker/Singularity) for portability.

**Results in this repo**
- `Alignment/ERR5743893.sorted.bam` — mapped reads aligned to `MN908947.fasta`.
- `QC_Reports/ERR5743893_1_fastqc.html` — example FastQC report.
- `Variant_Calling/ERR5743893.vcf` — example variant calls produced by `freebayes`.

**License**
- No license is included. 

---

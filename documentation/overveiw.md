## NGS Variant Calling Pipeline — Overview

*Project Documentation*

This project implements a Next‑Generation Sequencing (NGS) analysis pipeline for variant detection from paired‑end FASTQ reads. The workflow includes quality control, reference genome preparation, read alignment, BAM processing, and variant calling using standard bioinformatics tools.
*Purpose:* Reproducible SARS‑CoV‑2 WGS pipeline for a single paired‑end sample (ERR5743893). This document lists the commands used, short descriptions, prerequisites, expected outputs, and notes for reproducibility.

# Analysis Goals

1. Evaluate raw read quality

- Tool: FastQC

2. Quantify alignment coverage and mapping rates

- Tools: BWA, SAMtools

3. Discover single nucleotide polymorphisms (SNPs) and insertions/deletions (indels)

- Relative to the reference genome

- Tool: FreeBayes

# Pipeline Workflow Summary

Raw FASTQ → Quality Control → Read Alignment → BAM Processing → Variant Calling → VCF Output

# Key Metrics Generated

FastQC reports:
Per-base quality scores, adapter content, duplication levels

Alignment statistics:
Mapping rate, coverage depth across the ~30 kb genome

Variant file:
High-confidence SNPs/indels with quality scores and allele frequencies

# Learning Objectives

This project demonstrates a standard viral genomics workflow applicable to:

- Surveillance

- Outbreak investigation

- Lineage tracking

# Skills gained include:

- Linux command-line bioinformatics

- Reference-based genome assembly

- Population variant calling core competencies

## Prerequisites
- OS: Linux (Ubuntu or similar).
- Disk: allow ~5–10 GB for reference, FASTQ, and intermediates for this example.
- Tools: fastqc, bwa, samtools, bamtools, bcftools, vcftools, freebayes, wget, (or install freebayes via conda -c bioconda).


## 🧪 Sample & Reference Information

| *Parameter* | *Description* |
|---|---|
| *Sample ID* | ERR5743893 |
| *Sequencing Type* | Paired‑end |
| *Data Source* | ENA (European Nucleotide Archive) |
| *Reference Genome* | SARS‑CoV‑2 (MN908947.3) |

## 🛠 Software Requirements

| *Tool* | *Purpose* |
|---|---|
| *FastQC* | Raw read quality control |
| *BWA* | Sequence alignment |
| *SAMtools* | SAM/BAM file processing |
| *BAMtools* | BAM file utilities |
| *VCFtools* | Variant file analysis |
| *BCFtools* | Variant file manipulation |
| *FreeBayes* | Variant calling |
| *Conda* | Package management (recommended for freebayes) |

## 📁 Directory Structure


.
├── Raw_Data/
│   ├── ERR5743893_1.fastq
│   └── ERR5743893_2.fastq
│
├── Output/
│   ├── fastq/
│   ├── sam/
│   ├── bam/
│   ├── vcf/
│   └── reference/
│
└── README.md


## ⚙ Pipeline Overview (step-by-step)

Below are the commands used in the pipeline with short explanations.

*Step 1: Create Project Directories*

bash
mkdir -p Raw_Data Output/{fastq,sam,bam,vcf,reference}


Organizes raw data and analysis outputs.

*Step 2: Install Required Tools*

bash
sudo apt update
sudo apt install fastqc bwa samtools bamtools vcftools bcftools
conda install -c bioconda freebayes


Installs bioinformatics dependencies (use conda/bioconda for freebayes where preferred).

*Step 3: Download Raw FASTQ Files*

bash
cd Raw_Data
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_2.fastq.gz


*Step 4: Quality Control of Raw Reads*

bash
cp ERR5743893_1.fastq.gz ../Output/fastq/
cp ERR5743893_2.fastq.gz ../Output/fastq/

fastqc ERR5743893_1.fastq.gz --outdir ../Output/fastq/
fastqc ERR5743893_2.fastq.gz --outdir ../Output/fastq/


Generates FastQC reports to assess sequencing quality.

*Step 5: Decompress FASTQ Files*

bash
gunzip ERR5743893_1.fastq.gz
gunzip ERR5743893_2.fastq.gz
cd ..


Prepares FASTQ files for alignment.

*Step 6: Download & Index Reference Genome*

bash
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" -O Output/reference/MN908947.fasta

bwa index Output/reference/MN908947.fasta
samtools faidx Output/reference/MN908947.fasta


Indexes the reference genome for alignment and variant calling.

*Step 7: Read Alignment (BWA‑MEM)*

bash
bwa mem Output/reference/MN908947.fasta \
	Raw_Data/ERR5743893_1.fastq \
	Raw_Data/ERR5743893_2.fastq > Output/sam/ERR5743893.sam


Aligns sequencing reads to the reference genome.

*Step 8: Convert, Sort & Index BAM File*

bash
samtools view -@ 4 -Sb Output/sam/ERR5743893.sam > Output/bam/ERR5743893.bam
samtools sort -@ 4 -o Output/bam/ERR5743893.sorted.bam Output/bam/ERR5743893.bam
samtools index Output/bam/ERR5743893.sorted.bam


Produces a sorted and indexed BAM file.

*Step 9: Variant Calling (FreeBayes)*

bash
freebayes -f Output/reference/MN908947.fasta Output/bam/ERR5743893.sorted.bam > Output/vcf/ERR5743893.vcf


Detects SNPs and small indels, producing a VCF file.

*Step 10: Save Command History*

bash
history > Output/commands_used.txt


Ensures pipeline reproducibility.

## 📊 Final Output Files

| *File* | *Description* |
|---|---|
| ERR5743893.sorted.bam | Sorted alignment file |
| ERR5743893.sorted.bam.bai | BAM index |
| ERR5743893.vcf | Variant call file |
| *_fastqc.html | Quality control reports |
| commands_used.txt | Executed commands log |

## Notes & troubleshooting
- Keep gzipped original FASTQ copies in Output/fastq as canonical raw data.
- Use -@ N or -t N flags for samtools/bwa to use multiple CPU threads.
- If VCF contains many low‑quality calls, filter using bcftools filter or vcffilter (e.g., QUAL > 20 and min depth).
- For amplicon or tiled primer data, trim primers (e.g., ivar) before variant calling.
- To generate a consensus sequence: filter VCF for high‑confidence variants and apply bcftools consensus.



## ✅ Summary

This pipeline represents a standard viral genome variant‑calling workflow using open‑source bioinformatics tools. It is suitable for research, teaching, and reproducible genomics analyses. To scale or productionize, consider wrapping this into a workflow manager (Snakemake/Nextflow) and containerizing the environment (Docker/Singularity).
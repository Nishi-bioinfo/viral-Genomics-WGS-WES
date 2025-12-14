#!/usr/bin/env bash
set -euo pipefail

echo "SARS-CoV-2 Viral Genomics Analysis Pipeline"

# ================================
# 1. Create project directories
# ================================
echo "Creating directory structure..."
mkdir -p Raw_Data Output/{fastq,sam,bam,vcf,reference}

# ================================
# 2. Install required tools
# ================================
echo "Installing bioinformatics tools..."

sudo apt update
sudo apt install -y fastqc bwa samtools bamtools vcftools bcftools
conda install -c bioconda freebayes

# ================================
# 3. Download raw sequencing data
# ================================
echo "Downloading FASTQ files..."
cd Raw_Data

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_2.fastq.gz

# Save copies of raw FASTQ files
cp ERR5743893_1.fastq.gz ../Output/fastq/
cp ERR5743893_2.fastq.gz ../Output/fastq/

# ================================
# 4. Quality control
# ================================
echo "Running FastQC..."
fastqc ERR5743893_1.fastq.gz --outdir ../Output/fastq/
fastqc ERR5743893_2.fastq.gz --outdir ../Output/fastq/

# ================================
# 5. Decompress FASTQ files
# ================================
echo "Decompressing FASTQ files..."
gunzip ERR5743893_1.fastq.gz
gunzip ERR5743893_2.fastq.gz

cd ..

# ================================
# 6. Download and index reference genome
# ================================
echo "Downloading reference genome..."
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" \
  -O Output/reference/MN908947.fasta

echo "Indexing reference genome..."
bwa index Output/reference/MN908947.fasta
samtools faidx Output/reference/MN908947.fasta

# ================================
# 7. Read alignment
# ================================
echo "Aligning reads with BWA-MEM..."
bwa mem Output/reference/MN908947.fasta \
  Raw_Data/ERR5743893_1.fastq \
  Raw_Data/ERR5743893_2.fastq \
  > Output/sam/ERR5743893.sam

# ================================
# 8. BAM processing
# ================================
echo "Converting SAM to BAM..."
samtools view -@ 4 -Sb Output/sam/ERR5743893.sam > Output/bam/ERR5743893.bam

echo "Sorting BAM file..."
samtools sort -@ 4 -o Output/bam/ERR5743893.sorted.bam Output/bam/ERR5743893.bam

echo "Indexing BAM file..."
samtools index Output/bam/ERR5743893.sorted.bam

# ================================
# 9. Variant calling
# ================================
echo "Calling variants with FreeBayes..."
freebayes -f Output/reference/MN908947.fasta \
  Output/bam/ERR5743893.sorted.bam \
  > Output/vcf/ERR5743893.vcf

# ================================
# 10. Save command history
# ================================
history > Output/commands_used.txt

echo "SARS-CoV-2 Viral Genomics Analysis Pipeline Completed Successfully!"


#!/bin/bash

echo "----------------------------------------------------"
echo " Viral NGS Pipeline | FASTQC + BWA + SAMTOOLS + FREEBAYES"
echo "----------------------------------------------------"

# =========================
# Set variables
# =========================

SAMPLE_NAME="ERR5743893"

RAW_DIR="Raw_Data"
REF_DIR="Reference"
OUTPUT_DIR="Output"

FASTQ1="${RAW_DIR}/${SAMPLE_NAME}_1.fastq.gz"
FASTQ2="${RAW_DIR}/${SAMPLE_NAME}_2.fastq.gz"

REFERENCE_GENOME="${REF_DIR}/MN908947.fasta"

THREADS=4

# =========================
# STEP 0: Directory setup
# =========================

echo "STEP 0: Creating directory structure"

mkdir -p ${RAW_DIR} ${REF_DIR}
mkdir -p ${OUTPUT_DIR}/{fastq,sam,bam,vcf,qc}

echo "Directories created ✔"
echo

# =========================
# STEP 1: Install tools
# =========================

echo "STEP 1: Updating system & installing tools"

sudo apt update
sudo apt install -y fastqc bwa samtools bamtools bcftools vcftools freebayes

echo "Tools installed ✔"
echo

# =========================
# STEP 2: Download FASTQ
# =========================

echo "STEP 2: Download paired-end FASTQ files"

cd ${RAW_DIR}

wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_2.fastq.gz

cp *.fastq.gz ../${OUTPUT_DIR}/fastq/

cd ..

echo "FASTQ files downloaded ✔"
echo

# =========================
# STEP 3: Quality Control
# =========================

echo "STEP 3: Running FastQC"

fastqc ${FASTQ1} ${FASTQ2} --outdir ${OUTPUT_DIR}/qc

echo "FastQC completed ✔"
echo

# =========================
# STEP 4: Unzip FASTQ
# =========================

echo "STEP 4: Unzipping FASTQ files"

gunzip -k ${FASTQ1} ${FASTQ2}

echo "FASTQ files unzipped ✔"
echo

# =========================
# STEP 5: Reference genome
# =========================

echo "STEP 5: Download reference genome"

cd ${REF_DIR}

wget -nc "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" -O MN908947.fasta

echo "Indexing reference genome"

bwa index MN908947.fasta
samtools faidx MN908947.fasta

cd ..

echo "Reference ready ✔"
echo

# =========================
# STEP 6: Alignment
# =========================

echo "STEP 6: Aligning reads using BWA-MEM"

bwa mem -t ${THREADS} \
${REFERENCE_GENOME} \
${RAW_DIR}/${SAMPLE_NAME}_1.fastq \
${RAW_DIR}/${SAMPLE_NAME}_2.fastq \
> ${OUTPUT_DIR}/sam/${SAMPLE_NAME}.sam

echo "Alignment completed ✔"
echo

# =========================
# STEP 7: SAM → BAM
# =========================

echo "STEP 7: Converting SAM to BAM"

samtools view -@ ${THREADS} -Sb \
${OUTPUT_DIR}/sam/${SAMPLE_NAME}.sam \
> ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bam

echo "Conversion completed ✔"
echo

# =========================
# STEP 8: Sort & Index BAM
# =========================

echo "STEP 8: Sorting BAM file"

samtools sort -@ ${THREADS} \
-o ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.sorted.bam \
${OUTPUT_DIR}/bam/${SAMPLE_NAME}.bam

samtools index ${OUTPUT_DIR}/bam/${SAMPLE_NAME}.sorted.bam

echo "Sorted & indexed BAM ✔"
echo

# =========================
# STEP 9: Variant Calling
# =========================

echo "STEP 9: Variant calling using FreeBayes"

freebayes -f ${REFERENCE_GENOME} \
${OUTPUT_DIR}/bam/${SAMPLE_NAME}.sorted.bam \
> ${OUTPUT_DIR}/vcf/${SAMPLE_NAME}.vcf

echo "Variant calling completed ✔"
echo

# =========================
# STEP 10: Save commands
# =========================

echo "STEP 10: Saving command history"

history > ${OUTPUT_DIR}/commands_used.txt

echo "----------------------------------------------------"
echo " Viral NGS Pipeline completed successfully 🎉"
echo " Results:"
echo "  QC   : ${OUTPUT_DIR}/qc"
echo "  BAM  : ${OUTPUT_DIR}/bam"
echo "  VCF  : ${OUTPUT_DIR}/vcf"
echo "----------------------------------------------------"

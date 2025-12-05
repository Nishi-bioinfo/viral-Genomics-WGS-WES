# Viral WGS/WES Pipeline Documentation

This document describes the complete pipeline workflow for processing viral whole-genome (WGS) and whole-exome (WES) sequencing data.

## Pipeline Overview

The pipeline follows these steps:
1. **Setup** — Create directory structure
2. **Quality Control** — FastQC on raw reads
3. **Reference Preparation** — Download and index reference genome
4. **Read Alignment** — Map reads to reference using BWA
5. **Sorting & Indexing** — Convert SAM to sorted BAM and create index
6. **Variant Calling** — Call variants using FreeBayes

---

## Step 1: Setup & Directory Structure

Create the folder hierarchy for organizing inputs and outputs:

```bash
mkdir Raw_Data QC_Reports Alignment Variant_Calling Reference
```

**Folder purposes:**
- `Raw_Data/` — Raw paired-end FASTQ files (`.fastq.gz` or `.fastq`)
- `QC_Reports/` — FastQC HTML reports
- `Alignment/` — SAM and BAM alignment files
- `Variant_Calling/` — VCF files with variant calls
- `Reference/` — Reference genome FASTA and index files

---

## Step 2: Download Raw Reads

Download paired-end FASTQ files from ENA/SRA:

```bash
cd ./Raw_Data/
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_2.fastq.gz
cd ..
```

**Notes:**
- Replace `ERR5743893` with your actual run accession.
- Use `-nc` (no-clobber) to avoid re-downloading if files exist.
- For multiple samples, create a loop or use the provided `tools/download_reads.sh` script.

---

## Step 3: Quality Control (FastQC)

Install FastQC (if not already installed):

```bash
sudo apt install fastqc
```

Run FastQC on both read files:

```bash
fastqc Raw_Data/ERR5743893_1.fastq.gz --outdir QC_Reports
fastqc Raw_Data/ERR5743893_2.fastq.gz --outdir QC_Reports
```

**Output:**
- `QC_Reports/ERR5743893_1_fastqc.html` — HTML report for read 1
- `QC_Reports/ERR5743893_2_fastqc.html` — HTML report for read 2

**Interpretation:**
- Check for adapter contamination, quality scores, and GC content.
- Look for per-base quality warnings; if severe, consider trimming with `fastp` or `trimmomatic`.

---

## Step 4: Reference Preparation

### Download Reference Genome

Download the reference (example: SARS-CoV-2, accession MN908947):

```bash
cd Reference
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" -O MN908947.fasta
cd ..
```

### Install Alignment & Variant Calling Tools

```bash
sudo apt install bwa samtools bamtools bcftools vcftools freebayes
```

**Tools:**
- `bwa` — Burrows-Wheeler Aligner (short-read mapping)
- `samtools` — SAM/BAM file manipulation
- `bcftools` — VCF file processing
- `vcftools` — VCF statistics and filtering
- `freebayes` — Variant caller

### Index Reference

Create indexes for fast access:

```bash
bwa index Reference/MN908947.fasta
samtools faidx Reference/MN908947.fasta
```

**Output:**
- `Reference/MN908947.fasta.bwt`, `.pac`, `.ann`, `.amb`, `.sa` — BWA indexes
- `Reference/MN908947.fasta.fai` — FASTA index

---

## Step 5: Read Alignment

### Align Reads to Reference

Use BWA MEM to align paired-end reads:

```bash
bwa mem Reference/MN908947.fasta Raw_Data/ERR5743893_1.fastq Raw_Data/ERR5743893_2.fastq > Alignment/ERR5743893.sam
```

**Notes:**
- This produces a SAM (Sequence Alignment Map) file.
- If FASTQs are gzipped, use `<(gzip -dc file.fastq.gz)` for streaming.
- Adjust thread count with `-t` if available (e.g., `bwa mem -t 8 ...`).

### Convert SAM to BAM

Convert text SAM to binary BAM format (much smaller):

```bash
samtools view -@ 20 -S -b Alignment/ERR5743893.sam > Alignment/ERR5743893.bam
```

**Flags:**
- `-@` — number of threads (20 in example)
- `-S` — assume input is SAM
- `-b` — output BAM format

### Sort BAM

Sort BAM by genomic coordinate (required for most downstream analyses):

```bash
samtools sort -@ 4 -o Alignment/ERR5743893.sorted.bam Alignment/ERR5743893.bam
```

**Output:** `Alignment/ERR5743893.sorted.bam`

### Index BAM

Create a BAM index for fast random access:

```bash
samtools index Alignment/ERR5743893.sorted.bam
```

**Output:** `Alignment/ERR5743893.sorted.bam.bai`

---

## Step 6: Variant Calling

Use FreeBayes to call variants:

```bash
freebayes -f Reference/MN908947.fasta Alignment/ERR5743893.sorted.bam > Variant_Calling/ERR5743893.vcf
```

**Output:** `Variant_Calling/ERR5743893.vcf`

**Notes:**
- FreeBayes is sensitive to low-frequency variants (good for viral data).
- For faster calls on large genomes, consider `bcftools mpileup | bcftools call` or `samtools mpileup | bcftools call`.
- Optional: filter variants by quality/depth (see **Filtering** below).

---

## (Optional) Variant Filtering

Filter VCF by quality and depth:

```bash
bcftools filter -O z -o Variant_Calling/ERR5743893_filtered.vcf.gz \
  -s LOWQUAL -e '%QUAL<20 | DP<10' Variant_Calling/ERR5743893.vcf
bcftools index Variant_Calling/ERR5743893_filtered.vcf.gz
```

---

## Reproducibility Recommendations

### Use Conda Environment

Instead of system `apt install`, create a reproducible conda environment:

```bash
conda env create -f envs/viral-wgs-env.yaml
conda activate viral-wgs-env
```

### Version Pinning

Specify exact tool versions in `environment.yml`:

```yaml
name: viral-wgs-env
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - bwa=0.7.17
  - samtools=1.17
  - fastqc=0.11.9
  - freebayes=1.3.5
  - bcftools=1.17
```

### Logging

Capture exact commands and outputs for reproducibility:

```bash
log_file="logs/pipeline_$(date +%Y%m%d%H%M%S).log"
mkdir -p logs
bwa mem ... 2>&1 | tee -a "$log_file"
```

---

## Performance Tuning

| Step | Bottleneck | Optimization |
|------|-----------|--------------|
| Alignment | I/O, CPU | Increase threads (`bwa mem -t N`), use SSD |
| SAM→BAM | I/O | Increase samtools threads (`-@ N`) |
| Sorting | Memory, disk | Increase sort buffer (`-m 4G`), use temp SSD |
| Variant calling | CPU | None (single-threaded in FreeBayes) |

### Example: Multi-threaded workflow

```bash
bwa mem -t 8 Reference/MN908947.fasta R1.fastq R2.fastq \
  | samtools view -@ 8 -b - \
  | samtools sort -@ 8 -m 4G -o sorted.bam -
samtools index sorted.bam
freebayes -f Reference/MN908947.fasta sorted.bam > variants.vcf
```


---

## Next Steps

- **For many samples:** wrap commands in a loop or use Snakemake (`Snakefile`).
- **For trimming:** add `fastp` before alignment to remove adapters.
- **For variant annotation:** use `snpEff` or `VEP`.
- **For visualization:** use `IGV` or `samtools tview` for BAM inspection.

---

## References

- BWA: http://bio-bwa.sourceforge.net/
- SAMtools: http://samtools.sourceforge.net/
- FreeBayes: https://github.com/freebayes/freebayes
- FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- BCFtools: http://samtools.github.io/bcftools/


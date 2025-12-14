# *Detailed Pipeline: SARS‑CoV‑2 NGS Variant Calling*

*Short overview:* This document provides a reproducible, step‑by‑step pipeline to process paired‑end Illumina reads (sample ERR5743893) through quality control, alignment to the SARS‑CoV‑2 reference (MN908947.3), BAM processing, and variant calling (FreeBayes). It includes recommended commands, post‑processing, and troubleshooting tips.

*Key outputs:* *sorted BAM, **BAM index (.bai), **VCF, **FastQC reports, and a **commands log*.

---

## Prerequisites & Environment

- OS: *Linux* (Ubuntu/Debian tested)
- Disk: reserve *~5–10 GB* for this example (may need more for multiple samples)
- CPU/RAM: 4–8 cores and 8–32 GB RAM recommended for moderate samples
- Package managers: apt for system packages, conda (bioconda) recommended for bioinformatics tools

Install example (apt + conda):

bash
sudo apt update
sudo apt install -y fastqc bwa samtools bamtools bcftools vcftools wget
conda install -c bioconda freebayes


Check versions (record for reproducibility):

bash
bwa 2>&1 | head -n 1
samtools --version | head -n 1
freebayes --version
fastqc --version


---

## Workspace layout (recommended)


. 
├── Raw_Data/                   # downloaded raw FASTQ (.fastq or .fastq.gz)
├── Output/
│   ├── fastq/                  # raw gzipped FASTQ copies + FastQC reports
│   ├── reference/              # reference FASTA + indices
│   ├── sam/                    # raw SAM outputs (optional)
│   ├── bam/                    # BAM, sorted BAM, indexes
│   └── vcf/                    # raw and filtered VCFs
└── Documentation/
	└── pipeline.md


---

## Step-by-step pipeline (commands + explanations)

| *Step* | *Command* | *Notes / Purpose* |
|---:|---|---|
| 1 | mkdir -p Raw_Data Output/{fastq,sam,bam,vcf,reference} | Create directories used by the pipeline. |
| 2 | cd Raw_Data | Move into raw data folder. |
| 3 | wget -nc ftp://.../ERR5743893_1.fastq.gz | Download R1 (use -nc to avoid re-download). |
| 4 | wget -nc ftp://.../ERR5743893_2.fastq.gz | Download R2. |
| 5 | cp *.fastq.gz ../Output/fastq/ | Archive gzipped FASTQ in Output/fastq for provenance. |
| 6 | fastqc ERR5743893_1.fastq.gz --outdir ../Output/fastq/ | Run *FastQC* for quality control (do for both reads). |
| 7 | gunzip ERR5743893_1.fastq.gz | Decompress FASTQ for alignment (keep gzipped copy in output). |
| 8 | wget "https://.../MN908947.3?download=true" -O Output/reference/MN908947.fasta | Download *reference* FASTA. |
| 9 | bwa index Output/reference/MN908947.fasta | Build BWA index files. |
| 10 | samtools faidx Output/reference/MN908947.fasta | Create FASTA index for samtools. |
| 11 | bwa mem -t 4 Output/reference/MN908947.fasta Raw_Data/ERR5743893_1.fastq Raw_Data/ERR5743893_2.fastq > Output/sam/ERR5743893.sam | Align reads with *BWA‑MEM* (use -t for threads). |
| 12 | samtools view -@ 4 -Sb Output/sam/ERR5743893.sam > Output/bam/ERR5743893.bam | Convert SAM → BAM (compressed). |
| 13 | samtools sort -@ 4 -o Output/bam/ERR5743893.sorted.bam Output/bam/ERR5743893.bam | Sort BAM for downstream (coordinate sort). |
| 14 | samtools index Output/bam/ERR5743893.sorted.bam | Index sorted BAM (.bai). |
| 15 | freebayes -f Output/reference/MN908947.fasta Output/bam/ERR5743893.sorted.bam > Output/vcf/ERR5743893.vcf | Call *variants* with *FreeBayes* (see filters below). |
| 16 | history > Output/commands_used.txt | Save commands for provenance. |

---

## Post-processing & filtering (recommended)

1. Normalize and index VCFs with bcftools:

bash
bcftools norm -f Output/reference/MN908947.fasta -m -both Output/vcf/ERR5743893.vcf -Oz -o Output/vcf/ERR5743893.norm.vcf.gz
bcftools index Output/vcf/ERR5743893.norm.vcf.gz


2. Filter low-quality calls (example thresholds — tune to dataset):

bash
bcftools filter -e 'QUAL<20 || DP<10' Output/vcf/ERR5743893.norm.vcf.gz -Oz -o Output/vcf/ERR5743893.filtered.vcf.gz
bcftools index Output/vcf/ERR5743893.filtered.vcf.gz


3. Optional: annotate/filter by allele frequency, strand bias, or read position using vcftools/bcftools or snpEff.

4. Generate *consensus sequence* from filtered VCF:

bash
bcftools consensus -f Output/reference/MN908947.fasta Output/vcf/ERR5743893.filtered.vcf.gz > Output/ERR5743893.consensus.fasta


---

## Annotation

- Use *snpEff* or *VEP* to add functional consequences to variants.
- Example (snpEff):

bash
snpEff ann -v NC_045512.2 Output/vcf/ERR5743893.filtered.vcf.gz > Output/vcf/ERR5743893.annotated.vcf


---

## Troubleshooting & Tips

- If *alignment* fails: ensure reference FASTA is valid and indexed (check .fai and BWA index files). Use samtools faidx to rebuild FASTA index.
- If *FastQC* shows adapter contamination: run trimming (e.g., trim_galore, cutadapt) before alignment.
- If many low‑quality variant calls: increase QUAL and DP thresholds, require minimum allele balance, or use additional callers (LoFreq, iVar).
- For *amplicon* datasets: trim primer sequences before variant calling (use ivar or primers BED-based trimming).
- Save history and record exact *tool versions* for reproducibility.

---

## Reproducibility & Scaling

- For many samples, implement this workflow in a workflow manager (recommended: *Snakemake* or *Nextflow*).
- Containerize the environment (Docker/Singularity) using images from *biocontainers* or create a custom Dockerfile with pinned tool versions.

Example Snakemake rule (alignment):

python
rule bwa_mem:
	input:
		ref='Output/reference/MN908947.fasta',
		r1='Raw_Data/{sample}_1.fastq',
		r2='Raw_Data/{sample}_2.fastq'
	output:
		'Output/sam/{sample}.sam'
	threads: 4
	shell:
		'bwa mem -t {threads} {input.ref} {input.r1} {input.r2} > {output}'


---

## Example full run (copy‑paste)

bash
mkdir -p Raw_Data Output/{fastq,sam,bam,vcf,reference}
cd Raw_Data
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_1.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/ERR5743893/ERR5743893_2.fastq.gz
cp *.fastq.gz ../Output/fastq/
fastqc *.fastq.gz --outdir ../Output/fastq/
gunzip *.fastq.gz
cd ..
wget "https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true" -O Output/reference/MN908947.fasta
bwa index Output/reference/MN908947.fasta
samtools faidx Output/reference/MN908947.fasta
bwa mem -t 4 Output/reference/MN908947.fasta Raw_Data/ERR5743893_1.fastq Raw_Data/ERR5743893_2.fastq > Output/sam/ERR5743893.sam
samtools view -@ 4 -Sb Output/sam/ERR5743893.sam > Output/bam/ERR5743893.bam
samtools sort -@ 4 -o Output/bam/ERR5743893.sorted.bam Output/bam/ERR5743893.bam
samtools index Output/bam/ERR5743893.sorted.bam
freebayes -f Output/reference/MN908947.fasta Output/bam/ERR5743893.sorted.bam > Output/vcf/ERR5743893.vcf
history > Output/commands_used.txt


---
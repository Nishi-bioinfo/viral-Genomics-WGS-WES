#!/usr/bin/env bash
set -euo pipefail

# viral-wgs-pipeline.sh
# Complete pipeline for viral WGS/WES analysis
# Downloads reads, runs QC, aligns, sorts, and calls variants

LOGDIR="logs"
mkdir -p "$LOGDIR"
LOGFILE="$LOGDIR/pipeline_$(date +%Y%m%d%H%M%S).log"

# Configuration defaults
SAMPLE="${1:-ERR5743893}"
REFERENCE="Reference/MN908947.fasta"
THREADS="${2:-4}"
DRY_RUN=false

usage(){
  cat <<EOF
Usage: $0 [SAMPLE] [THREADS] [OPTIONS]

Positional:
  SAMPLE      ENA run accession (default: ERR5743893)
  THREADS     CPU threads (default: 4)

Options:
  --dry-run   show commands without executing
  --no-download   skip FASTQ download
  -h, --help  show this help

Example:
  ./viral-wgs-pipeline.sh ERR5743893 8
  ./viral-wgs-pipeline.sh ERR5743893 4 --dry-run
EOF
}

log(){ echo "[INFO] $(date -Iseconds) - $*" | tee -a "$LOGFILE"; }
err(){ echo "[ERROR] $(date -Iseconds) - $*" | tee -a "$LOGFILE" >&2; }
warn(){ echo "[WARN] $(date -Iseconds) - $*" | tee -a "$LOGFILE"; }

run(){
  echo "+ $*" | tee -a "$LOGFILE"
  if ! $DRY_RUN; then
    eval "$@" 2>&1 | tee -a "$LOGFILE" || { err "Command failed: $*"; return 1; }
  fi
}

# Parse arguments
SKIP_DOWNLOAD=false
shift 2 || true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=true; shift;;
    --no-download) SKIP_DOWNLOAD=true; shift;;
    -h|--help) usage; exit 0;;
    *) err "Unknown arg: $1"; usage; exit 1;;
  esac
done

log "=== Starting viral WGS/WES pipeline ==="
log "Sample: $SAMPLE"
log "Reference: $REFERENCE"
log "Threads: $THREADS"
log "Log file: $LOGFILE"

# ============================================================================
# STEP 1: Setup directories
# ============================================================================
log ""
log "=== STEP 1: Creating directory structure ==="
run "mkdir -p Raw_Data QC_Reports Alignment Variant_Calling Reference"

# ============================================================================
# STEP 2: Download FASTQ files
# ============================================================================
log ""
log "=== STEP 2: Downloading raw reads ==="

R1="Raw_Data/${SAMPLE}_1.fastq.gz"
R2="Raw_Data/${SAMPLE}_2.fastq.gz"

if $SKIP_DOWNLOAD; then
  log "Skipping download (--no-download specified)"
  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    err "FASTQ files not found: $R1 or $R2"
    exit 1
  fi
else
  if [[ ! -f "$R1" ]]; then
    log "Downloading $SAMPLE read 1"
    run "wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/${SAMPLE}/${SAMPLE}_1.fastq.gz -P Raw_Data/"
  else
    log "Found existing: $R1"
  fi
  
  if [[ ! -f "$R2" ]]; then
    log "Downloading $SAMPLE read 2"
    run "wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR574/003/${SAMPLE}/${SAMPLE}_2.fastq.gz -P Raw_Data/"
  else
    log "Found existing: $R2"
  fi
fi

# ============================================================================
# STEP 3: Quality Control (FastQC)
# ============================================================================
log ""
log "=== STEP 3: Running FastQC ==="

# Check for fastqc
if ! command -v fastqc &> /dev/null; then
  warn "fastqc not found. Install with: apt install fastqc or conda install -c bioconda fastqc"
  err "Cannot continue without fastqc"
  exit 1
fi

run "fastqc $R1 --outdir QC_Reports"
run "fastqc $R2 --outdir QC_Reports"
log "FastQC reports saved to QC_Reports/"

# ============================================================================
# STEP 4: Reference preparation
# ============================================================================
log ""
log "=== STEP 4: Preparing reference genome ==="

if [[ ! -f "$REFERENCE" ]]; then
  log "Reference not found. Downloading MN908947 (SARS-CoV-2)..."
  run "wget 'https://www.ebi.ac.uk/ena/browser/api/fasta/MN908947.3?download=true' -O $REFERENCE"
else
  log "Reference found: $REFERENCE"
fi

# Check for alignment tools
for tool in bwa samtools freebayes; do
  if ! command -v "$tool" &> /dev/null; then
    warn "$tool not found. Install with: apt install $tool or conda install -c bioconda $tool"
    err "Cannot continue without $tool"
    exit 1
  fi
done

# Index reference if needed
if [[ ! -f "${REFERENCE}.bwt" ]]; then
  log "Creating BWA index for reference"
  run "bwa index $REFERENCE"
else
  log "BWA index exists"
fi

if [[ ! -f "${REFERENCE}.fai" ]]; then
  log "Creating FASTA index for reference"
  run "samtools faidx $REFERENCE"
else
  log "FASTA index exists"
fi

# ============================================================================
# STEP 5: Read alignment
# ============================================================================
log ""
log "=== STEP 5: Aligning reads ==="

SAM="Alignment/${SAMPLE}.sam"
BAM="Alignment/${SAMPLE}.bam"
SORTED_BAM="Alignment/${SAMPLE}.sorted.bam"
BAI="${SORTED_BAM}.bai"

log "Aligning with bwa mem"
run "bwa mem -t $THREADS $REFERENCE <(gzip -dc $R1) <(gzip -dc $R2) > $SAM"
log "Alignment saved to $SAM"

# ============================================================================
# STEP 6: Convert and sort BAM
# ============================================================================
log ""
log "=== STEP 6: Converting SAM to sorted BAM ==="

log "Converting SAM to BAM format"
run "samtools view -@ $THREADS -S -b $SAM > $BAM"

log "Sorting BAM by genomic coordinate"
run "samtools sort -@ $THREADS -o $SORTED_BAM $BAM"

log "Creating BAM index"
run "samtools index $SORTED_BAM"

log "Sorted BAM and index created: $SORTED_BAM, $BAI"

# Optional: remove intermediate SAM file
if ! $DRY_RUN; then
  log "Removing intermediate SAM file"
  rm "$SAM"
fi

# ============================================================================
# STEP 7: Variant calling
# ============================================================================
log ""
log "=== STEP 7: Calling variants with FreeBayes ==="

VCF="Variant_Calling/${SAMPLE}.vcf"

run "freebayes -f $REFERENCE $SORTED_BAM > $VCF"
log "Variants called and saved to $VCF"

# ============================================================================
# Summary
# ============================================================================
log ""
log "=== Pipeline complete ==="
log "Summary of outputs:"
log "  QC Reports:  QC_Reports/${SAMPLE}_*_fastqc.html"
log "  Alignment:   $SORTED_BAM"
log "  Index:       $BAI"
log "  Variants:    $VCF"
log ""
log "Next steps (optional):"
log "  - Filter variants:  bcftools filter -e '%QUAL<20' $VCF"
log "  - Annotate:         snpEff or VEP"
log "  - Visualize BAM:    samtools tview $SORTED_BAM $REFERENCE"
log "  - View VCF:         cat $VCF | head -50"
log ""
log "Logs saved to: $LOGFILE"


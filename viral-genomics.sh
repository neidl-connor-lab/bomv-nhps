#!/bin/bash
# sbatch options
#SBATCH -A MIP24002
#SBATCH -p skx
#SBATCH -N 8
#SBATCH -n 8
#SBATCH -o log-%x.out
#SBATCH -t 12:00:00

## setup -----------------------------------------------------------------------
# functions
mesg () { echo -e "[MSG] $@"; }
err () { echo -e "[ERR] $@"; exit 1; }
checkcmd () {
  if [ $? -eq 0 ]
  then
    mesg "$@ succeeded"
  else
    err "$@ failed"
  fi
}

# pre-set variables
VISUALIZATION="pipeline/visualization.r"
ANNOTATION="pipeline/mutationannotation.r"
K2DB="pipeline/kraken2db"
THRESHOLD="100"

# help message
HELP="usage: sbatch -J JOBNAME $(basename "$0") OPTIONS
VISUALIZATION AND ANNOTATION SCRIPTS MUST BE IN pipeline DIRECTORY

options (REQUIRED):
  -n NHP/host genome bowtie2 index
  -v viral genome bowtie2 index
  -f viral genome FASTA
  -g viral genome annotation GTF
  -o output directory
  -s sample ID
  -x forward FASTQ file
  -y [OPTIONAL] reverse FASTQ file
  -t [OPTIONAL] minimum aligned read depth (default: $THRESHOLD)
  -k [OPTIONAL] kraken2 database (default: $K2DB)
  -h show this message and exit
"
# parsing arguments
while getopts ":hn:v:f:g:o:s:x:y:t:k:" opt
do
  case ${opt} in
    n ) HOST="${OPTARG}"
      ;;
    v ) VIRAL="${OPTARG}"
      ;;
    f ) REFSEQ="${OPTARG}"
      ;;
    g ) GTF="${OPTARG}"
      ;;
    o ) ODIR="${OPTARG}"
      ;;
    s ) SAMPLE="${OPTARG}"
      ;;
    x ) R1="${OPTARG}"
      ;;
    y ) R2="${OPTARG}"
      ;;
    t ) THRESHOLD="${OPTARG}"
      ;;
    k ) K2DB="${OPTARG}"
      ;;
    h ) echo "$HELP" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

## logfile header -------------------------------------------------------------
echo "========================JOB INFO========================="
echo "Job ID: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Start date: $(date)"
echo "Current directory: $PWD"
echo "========================================================="
echo ""

## check inputs ----------------------------------------------------------------
mesg "STEP 0: CHECKING INPUTS"

# set localectl
export LC_ALL=C
export LANG=C

# load bowtie2
module unload xalt
module load biocontainers
module load bowtie2
module load samtools
module load kraken2
module load lofreq
module use /work/projects/singularity/rstudio/lua
module load Rstats/4.1.3

# double check that the visualization and annotations scripts exist
if [ ! -f "$VISUALIZATION" ]
then
  err "Visualization script not found: $VISUALIZATION"
fi
if [ ! -f "$ANNOTATION" ]
then
  err "SNV annotation script not found: $ANNOTATION"
fi

# double-check that kraken2 database exists
if [ -z "$K2DB" ]
then
  err "No Kraken2 database provided"
elif [ ! -d "$K2DB" ]
then
  err "Invalid Kraken2 database: $K2DB"
fi

# bowtie2 indicies should have 6 files
if [ -z "$HOST" ]
then
  err "No host index provided"
elif [ "$(ls -1 ${HOST}.* 2> /dev/null | wc -l)" -eq 6 ]
then
  mesg "Host index: $HOST"
else
  err "Invalid host index: $HOST"
fi
if [ -z "$VIRAL" ]
then
  err "No viral index provided"
elif [ "$(ls -1 ${VIRAL}.* 2> /dev/null | wc -l)" -eq 6 ]
then
  mesg "Viral index: $VIRAL"
else
  err "Invalid viral index: $VIRAL"
fi

# reference sequence
if [ -z "$REFSEQ" ]
then
  err "No viral reference FASTA provided"
elif [ -f "$REFSEQ" ]
then
  mesg "Reference FASTA: $REFSEQ"
else
  err "Invalid reference FASTA: $REFSEQ"
fi

# annotation GTF
if [ -z "$GTF" ]
then
  err "No viral annotation GTF provided"
elif [ -f "$GTF" ]
then
  mesg "Reference annotation: $GTF"
else
  err "Invalid viral annotation GTF: $GTF"
fi

# R1/R0 FASTQ file
if [ -z "$R1" ]
then
  err "No FASTQ file provided (-x)"
elif [ -f "$R1" ]
then
  mesg "First FASTQ file: $R1"
else
  err "Invalid first FASTQ file: $R1"
fi

# check for R2 FASTQ file
if [ -z "$R2" ]
then 
  mesg "Only 1 FASTQ file detected. Running as unpaired."
elif [ -f "$R2" ]
then
  mesg "Second FASTQ file: $R2"
else
  err "Invalid second FASTQ file: $R2"
fi

# output directory
if [ -z "$ODIR" ]
then
  err "No output directory provided"
elif [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi

# sample ID
if [ -z "$SAMPLE" ]
then
  err "No sample ID provided"
else
  mesg "Working with sample: $SAMPLE"
  ODIR="$ODIR/$SAMPLE"
  # make output directory
  mkdir -p "$ODIR"
fi

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## alignment to host ---------------------------------------------------------
mesg "STEP 1: DISCARD HOST READS"

# build command based on whether has paired reads
if [ -z "$R2" ] # unpaired
then
  CMD="bowtie2 --threads 8 --very-sensitive -x '$HOST' --un-conc-gz '$ODIR/filtered.fq.gz' -U '$R1' > /dev/null"
else
  CMD="bowtie2 --threads 8 --very-sensitive -x '$HOST' --un-conc-gz '$ODIR/filtered-r%.fq.gz' -1 '$R1' -2 '$R2' > /dev/null"
fi

# run alignment filter
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Host alignment" 
echo ""

# update R1 and R2
if [ -z "$R2" ] # unpaired
then
  R1="$ODIR/filtered.fq.gz"
else
  R1="$ODIR/filtered-r1.fq.gz"
  R2="$ODIR/filtered-r2.fq.gz"
fi

## metagenomic quantification on non-host reads ------------------------------
mesg "STEP 2: METAGENOMIC CLASSIFICATION"

# build command on whether has paired reads
CMD="kraken2 --threads 8 --db '$K2DB' --output - --report '$ODIR/metagenomics-raw.tsv' --use-names --gzip-compressed"
if [ -z "$R2" ] # unpaired
then
  CMD="$CMD '$R1'"
else 
  CMD="$CMD --paired '$R1' '$R2'"
fi

# run metagenomic classification
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Metagenomic classification"
echo ""

# only keep entries with % > 0.00
cat "$ODIR/metagenomics-raw.tsv" | awk '($1!=0.00) { print }' > "$ODIR/metagenomics.tsv"
rm "$ODIR/metagenomics-raw.tsv"

## alignment to virus --------------------------------------------------------
mesg "STEP 3: ALIGN TO VIRUS"

# build command based on whether has paired reads
if [ -z "$R2" ] # unpaired
then
  CMD="bowtie2 --threads 8 -x '$VIRAL' -U '$R1' > '$ODIR/alignment.sam'"
else # paired
  CMD="bowtie2 --threads 8 -x '$VIRAL' -1 '$R1' -2 '$R2' > '$ODIR/alignment.sam'"
fi

# run alignment
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Viral alignment"
echo ""

## process BAM -----------------------------------------------------------------
mesg "STEP 4: PROCESS ALIGNMENT"

# compress SAM to BAM
CMD="samtools view --threads 8 -b -h '$ODIR/alignment.sam' > '$ODIR/alignment-raw.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Compression"
rm "$ODIR/alignment.sam"

# fill in coordinates
CMD="samtools fixmate --threads 8 -m '$ODIR/alignment-raw.bam' '$ODIR/alignment-fixmate.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Filling in coordinates"
rm "$ODIR/alignment-raw.bam"

# sort by coordinate
CMD="samtools sort --threads 8 '$ODIR/alignment-fixmate.bam' > '$ODIR/alignment-sorted.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Sorting"
rm "$ODIR/alignment-fixmate.bam"

# mark duplicates
CMD="samtools markdup --threads 8 '$ODIR/alignment-sorted.bam' '$ODIR/alignment-dupmark.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Marking duplicates"
rm "$ODIR/alignment-sorted.bam"

# score indels to get final BAM
CMD="lofreq indelqual --dindel --ref '$REFSEQ' '$ODIR/alignment-dupmark.bam' > '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indelqual"
rm "$ODIR/alignment-dupmark.bam"

# index final BAM
CMD="samtools index '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indexing"
echo ""

## calculate coverage ----------------------------------------------------------
mesg "STEP 5: CALCULATE COVERAGE"

# coverage with samtools depth
CMD="samtools depth --threads 8 --excl-flags DUP -a -H '$ODIR/alignment.bam' > '$ODIR/coverage.tsv'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Coverage"
echo ""

## assemble consensus ----------------------------------------------------------
mesg "STEP 6: ASSEMBLE CONSENSUS"

# consensus with samtools
CMD="samtools consensus --threads 8 -a --excl-flags DUP --use-qual --min-depth $THRESHOLD --call-fract 0.5 --mode simple  --output '$ODIR/consensus.fa' '$ODIR/alignment.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Consensus"
echo ""

## quantify SNVs ---------------------------------------------------------------
mesg "STEP 7: QUANTIFY SNVs"

# lofreq doesn't have a DUP filter, so make a tmp BAM with DUPs removed
CMD="samtools view --threads 8 --bam --excl-flags DUP '$ODIR/alignment.bam' > '$ODIR/alignment-dedup.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Read deduplication"

# index tmp BAM
CMD="samtools index '$ODIR/alignment-dedup.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indexing deduplicated BAM"

# last check: there should be >0 mapped reads in the BAM
if [ $(samtools view -c -F 260 "$ODIR/alignment-dedup.bam" 2> /dev/null) -eq 0 ]
then
  mesg "No reads mapped to the viral genome; exiting..."

  # remove deduplicated BAM and BAI
  rm "$ODIR/alignment-dedup.bam"
  rm "$ODIR/alignment-dedup.bam.bai"

  # print modules and exit
  module list
  mesg "JOB COMPLETE!"
  echo ""
  exit 0
fi

# run lofreq
# keeping mapping quality parameters same between samtools and LoFreq
CMD="lofreq call-parallel --pp-threads 8 --call-indels --min-cov $THRESHOLD --ref '$REFSEQ' '$ODIR/alignment-dedup.bam' > '$ODIR/snvs.vcf'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "LoFreq"

# clean up deduplicated BAM and BAI
rm "$ODIR/alignment-dedup.bam"
rm "$ODIR/alignment-dedup.bam.bai"

# visualize coverage and SNV profiles
CMD="Rscript $VISUALIZATION --dir '$ODIR' --annotation '$GTF' --threshold $THRESHOLD"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Visualization"

# run mutation annotation script
CMD="Rscript $ANNOTATION --dir '$ODIR' --annotation '$GTF' --genome '$REFSEQ'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "SNV annotation"
echo ""

## package version -------------------------------------------------------------
mesg "Pipeline complete! Printing package versions..."
module list
mesg "JOB COMPLETE!"
echo ""


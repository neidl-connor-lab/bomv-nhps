#!/bin/bash -l

# qsub options
#$ -l h_rt=12:00:00
#$ -pe omp 8
#$ -j y
#$ -o log-$JOB_NAME.qlog

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
LOFREQ="pipeline/lofreq/lofreq"
# help message
HELP="usage: qsub -P PROJECT -N JOBNAME $0 -i INDEX -f FASTA -o ODIR -s SAMPLE -x R1 [-y R2]
Please submit the job from the pipeline directory!

arguments:
  -i bowtie2 index path and prefix
  -f reference FASTA
  -o output directory
  -s sample ID
  -x FASTQ file; R1 file if paired reads
  -y [OPTIONAL] R2 FASTQ file if paired reads
  -h print this message and exit
"

# parsing arguments
while getopts ":hi:f:o:x:y:s:" opt 
do 
  case ${opt} in 
    i ) IDX="${OPTARG}"
      ;;
    f ) REFSEQ="${OPTARG}"
      ;;
    o ) ODIR="${OPTARG}"
      ;;
    x ) R1="${OPTARG}"
      ;;
    y ) R2="${OPTARG}"
      ;;
    s) SAMPLE="${OPTARG}"
      ;;
    h ) echo "$HELP" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

## print job info for output log -----------------------------------------------
echo "=========================================================="
echo "Start date: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

## check inputs ----------------------------------------------------------------
mesg "STEP 0: CHECKING INPUTS"

# load bowtie2
module load bowtie2/2.4.2
checkcmd "Loading bowtie2/2.4.2"

# load samtools
module load htslib/1.18
module load samtools/1.18
checkcmd "Loading samtools/1.18"

# double-check that lofreq exists
if [ -z "$($LOFREQ version 2> /dev/null)" ]
then
  err "LoFreq error: $LOFREQ"
fi

# bowtie2 index should have 6 files
if [ -z "$IDX" ]
then
  err "No bowtie2 index provided"
elif [ "$(ls -1 ${IDX}.* 2> /dev/null | wc -l)" -eq 6 ]
then
  mesg "Bowtie2 index: $IDX"
else
  err "Invalid bowtie2 index; run setup.sh first: $IDX"
fi

# reference sequence
if [ -z "$REFSEQ" ]
then
  err "No reference FASTA provided"
elif [ -f "$REFSEQ" ]
then
  mesg "Reference FASTA: $REFSEQ"
else
  err "Invalid reference FASTA: $REFSEQ"
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

# sample ID
if [ -z "$SAMPLE" ]
then
  err "No sample ID provided"
else
  mesg "Working with sample: $SAMPLE"
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

# done checking inputs!
mesg "Done checking inputs!"
echo ""

## alignment -------------------------------------------------------------------
mesg "STEP 1: ALIGN TO GENOME"

# build command based on whether has paired reads
if [ -z "$R2" ] # unpaired
then
  CMD="bowtie2 --threads 8 -x '$IDX' -U '$R1' > '$ODIR/$SAMPLE.sam'"
else # paired
  CMD="bowtie2 --threads 8 -x '$IDX' -1 '$R1' -2 '$R2' > '$ODIR/$SAMPLE.sam'"
fi

# run alignment
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Alignment"

# compress SAM to BAM
CMD="samtools view --threads 8 -b -h '$ODIR/$SAMPLE.sam' > '$ODIR/$SAMPLE-raw.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Compression"
rm "$ODIR/$SAMPLE.sam"
echo ""

## process BAM -----------------------------------------------------------------
mesg "STEP 2: PROCESS ALIGNMENT"

# sort BAM
CMD="samtools sort --threads 8 '$ODIR/$SAMPLE-raw.bam' > '$ODIR/$SAMPLE-sorted.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Sorting"
rm "$ODIR/$SAMPLE-raw.bam"

# score indels to get final BAM
CMD="$LOFREQ indelqual --dindel --ref '$REFSEQ' '$ODIR/$SAMPLE-sorted.bam' > '$ODIR/$SAMPLE.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indelqual"
rm "$ODIR/$SAMPLE-sorted.bam"

# index final BAM
CMD="samtools index '$ODIR/$SAMPLE.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Indexing"
echo ""

## calculate coverage ----------------------------------------------------------
mesg "STEP 3: CALCULATE COVERAGE"

# coverage with samtools depth
CMD="samtools depth --threads 8 -a -H '$ODIR/$SAMPLE.bam' > '$ODIR/$SAMPLE.tsv'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Coverage"
echo ""

## assemble consensus ----------------------------------------------------------
mesg "STEP 4: ASSEMBLE CONSENSUS"

# consensus with samtools
CMD="samtools consensus --threads 8 -a --use-qual --min-depth 10 --call-fract 0.5 --output '$ODIR/$SAMPLE-tmp.fa' '$ODIR/$SAMPLE.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Consensus"

# update consensus header
mesg "Updating consensus header"
echo ">$SAMPLE" > "$ODIR/$SAMPLE.fa"
cat "$ODIR/$SAMPLE-tmp.fa" | grep "^[^>]" >> "$ODIR/$SAMPLE.fa"
rm "$ODIR/$SAMPLE-tmp.fa"
echo ""

## quantify SNVs ---------------------------------------------------------------
mesg "STEP 5: QUANTIFY SNVs"

# run lofreq
# keeping mapping quality parameters same between samtools and LoFreq
CMD="$LOFREQ call-parallel --pp-threads 8 --call-indels --min-cov 10 --ref '$REFSEQ' '$ODIR/$SAMPLE.bam' > '$ODIR/$SAMPLE.vcf'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "LoFreq"
echo ""

## package version -------------------------------------------------------------
mesg "Pipeline complete! Printing package versions..."
module list
echo "LoFreq"
$LOFREQ version
echo ""

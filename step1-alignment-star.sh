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

# default values and help message
ARGS="--runThreadN 8 --runMode alignReads --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMunmapped Within"
HELP="usage: qsub -P PROJECT -N JOBNAME $(basename "$0") OPTIONS

options (REQUIRED):
  -i index directory
  -o output directory
  -s sample ID
  -x forward FASTQ file
  -y reverse FASTQ file (optional)
  -h show this message and exit
"

# parsing arguments
while getopts ":hi:o:s:x:y:" opt 
do 
  case ${opt} in 
    i ) INDEX="${OPTARG}"
      ;;
    o ) ODIR="${OPTARG}"
      ;;
    s ) SAMPLE="${OPTARG}"
      ;;
    x ) R1="${OPTARG}"
      ;;
    y ) R2="${OPTARG}"
      ;;
    h ) echo "${HELP}" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

## job info --------------------------------------------------------------------
echo "=========================================================="
echo "Start date: $(date)"
echo "Running on node: $(hostname)"
echo "Current directory: $(pwd)"
echo "Job name: $JOB_NAME"
echo "Job ID: $JOB_ID"
echo "=========================================================="
echo ""

## check inputs ----------------------------------------------------------------
# index directory
if [ -z "$INDEX" ]
then
  err "No index directory provided"
elif [ -d "$INDEX" ]
then
  mesg "Valid index directory: $INDEX"
  ARGS="${ARGS} --genomeDir '${INDEX}'"
else 
  err "Invalid index directory: $INDEX"
fi

# sample ID
if [ -z "$SAMPLE" ]
then
  err "No sample ID provided"
else 
  mesg "Using sample ID: $SAMPLE"
fi

# R1 FASTQ
if [ -z "$R1" ]
then
  err "No forward FASTQ provided."
elif [ -f "$R1" ]
then
  mesg "Valid forward FASTQ: $R1"
else
  err "Invalid forward FASTQ: $R1"
fi

# R2 FASTQ
if [ -z "$R2" ]
then
  mesg "No reverse FASTQ file provided. Running in single-read mode."
elif [ -f "$R2" ]
then
  mesg "Valid reverse FASTQ: $R2"
else
  err "Invalid reverse FASTQ: $R2"
fi

# output directory
if [ -z "$ODIR" ]
then
  err "No output directory provided."
elif [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi

mesg "Done checking inputs."
echo ""

## set up and run STAR ---------------------------------------------------------
# load STAR
mesg "Loading STAR"
module load star/2.7.1a

# allow more virtual memory
ulimit -v 31428973073

# run STAR alignment
mesg "Running STAR alignment"
if [ -f "$R2" ] # paired
then
  CMD="STAR ${ARGS} --readFilesIn '$R1' '$R2' --outFileNamePrefix '$ODIR/$SAMPLE-'"
else
  CMD="STAR ${ARGS} --readFilesIn '$R1' --outFileNamePrefix '$ODIR/$SAMPLE-'"
fi
mesg "CMD: ${CMD}"
eval "${CMD}"
checkcmd "STAR alignment"

# rename BAM and log file
mv "$ODIR/$SAMPLE-Aligned.out.bam" "$ODIR/$SAMPLE.bam"
mv "$ODIR/$SAMPLE-Log.final.out" "$ODIR/${SAMPLE}Log.final.out"

# remove remaining files
rm -r $ODIR/$SAMPLE-*
echo ""

## done! -----------------------------------------------------------------------
module list
echo "STAR $(STAR --version)"
echo ""
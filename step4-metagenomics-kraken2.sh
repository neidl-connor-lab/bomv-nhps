#!/bin/bash -l

# qsub options
#$ -l h_rt=12:00:00
#$ -pe omp 8
#$ -l mem_total=1000G
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
KDB="pipeline/kraken2db"
ARGS="--threads 8 --output - --use-names --gzip-compressed"
HELP="usage: qsub -P project -N JOBNAME $(basename "$0") [-k KDB] -o ODIR -s SAMPLE -x R1 [-y R2]

options (default):
  -k Kraken2 database path ($KDB)
  -o output directory
  -s sample ID
  -x FASTQ file; R1 file if paired reads
  -y [OPTIONAL] R2 FASTQ file if paired reads
  -h show this message and exit
"

# parsing arguments
while getopts ":hk:o:s:x:y:" opt 
do 
  case ${opt} in 
    k ) KDB="${OPTARG}"
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
# check Kraken2 is in $PATH
if [ -z "$(which kraken2 2> /dev/null)" ]
then
  err "Kraken2 not found in your path"
fi

# check Kraken2 db
if [ -z "$KDB" ]
then
  err "Kraken2 database not provided"
elif [ -d "$KDB" ]
then
  mesg "Using Kraken2 DB: $KDB"
  ARGS="$ARGS --db '$KDB'"
else
  err "Invalid Kraken2 DB: $KDB"
fi

# sample ID
if [ -z "$SAMPLE" ]
then
  err "No sample ID provided"
else
  mesg "Working with $SAMPLE"
fi

# forward read file
if [ -z "$R1" ]
then
  err "No forward read FASTQ provided"
elif [ -f "$R1" ]
then 
  mesg "Forward FASTQ: $R1"
else
  err "Invalid forward FASTQ: $R1"
fi

# reverse read file
if [ -z "$R2" ]
then
  # running in unpaired mode
  mesg "No reverse read FASTQ provided; running in unpaired mode."
elif [ -f "$R2" ]
then 
  mesg "Reverse FASTQ: $R2"
  ARGS="$ARGS --paired"
else
  err "Invalid reverse FASTQ: $R2"
fi

# check output directory
if [ -z "$ODIR" ]
then
  err "Output directory not provided"
elif [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi

mesg "Done checking inputs."
echo ""

## set up and run command ------------------------------------------------------
mesg "Running Kraken2 on $SAMPLE"

# set up command whether paired/unpaired
if [ -z "$R2" ]
then
  # unpaired
  CMD="kraken2 $ARGS --report '$ODIR/$SAMPLE.tsv' '$R1'"
else
  # paired
  CMD="kraken2 $ARGS --report '$ODIR/$SAMPLE.tsv' '$R1' '$R2'"
fi

mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Kraken2"
echo ""

## done! -----------------------------------------------------------------------
module list
kraken2 --version
echo ""

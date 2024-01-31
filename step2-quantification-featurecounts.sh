#!/bin/bash -l

# qsub options
#$ -l h_rt=12:00:00
#$ -pe omp 8
#$ -j y
#$ -o log-$JOB_NAME.qlog

## setup -----------------------------------------------------------------------
# functions
mesg () { echo "[MSG] $@"; }
err () { echo "[ERR] $@"; exit 1; }
checkcmd () {
  if [ $? -eq 0 ]
  then
    mesg "$@ succeeded"
  else
    err "$@ failed"
  fi
}

# default values and help message
# -T 8: 8 threads
# -O: assign reads to all overlapping meta-features
# -M: count multi-mapping reads
ARGS="-T 8 -O -M -p"
FEAT="transcript"
META="gene_name"
HELP="usage: qsub -P PROJECT -N JOBNAME $(basename "$0") [-f FEAT] [-m META] -o OFILE -g GTF BAM [...BAM]

positional arguments:
  BAM path to alignment BAM files

arguments (default):
  -f feature to count ($FEAT)
  -m meta-feature for collaping counts ($META)
  -o output file
  -g GTF annotation file
  -h display this message and exit
"

# parsing arguments
while getopts ":ho:f:m:g:" opt
do
  case ${opt} in
    o ) OFILE="${OPTARG}"
      ;;
    f ) FEAT="${OPTARG}"
      ;;
    m ) META="${OPTARG}"
      ;;
    g ) GTF="${OPTARG}" 
      ;;
    h ) echo "${HELP}" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))
BAM="$@"

# job info
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

## check inputs ----------------------------------------------------------------
# check GTF
if [ -z "$GTF" ]
then
  err "GTF file required."
elif [ -f "$GTF" ]
then 
  mesg "Valid GTF file: $GTF"
  ARGS="$ARGS -a '$GTF'"
else
  err "Unable to find GTF: $GTF"
fi

# check output file given
if [ -z "$OFILE" ]
then
  err "No output filename provided"
fi
# check output directory
ODIR="$(dirname $OFILE)"
if [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi
ARGS="$ARGS -o '$OFILE'"

# set up -t and -g
mesg "Using feature: $FEAT"
mesg "Using meta-feature: $META"
ARGS="$ARGS -t $FEAT -g $META"

# get and check BAM files
if [ -z "$BAM" ]
then
  err "No BAM path provided."
fi
for i in $BAM
do
  if [ ! -f "$i" ]
  then
    err "Invalid BAM file: $i"
  fi
done
mesg "Valid BAM files: $BAM"
ARGS="$ARGS $BAM"
mesg "Done checking inputs."
echo ""

## run featureCounts -----------------------------------------------------------
# load subread
mesg "Loading featureCounts"
module load subread/1.6.2

# run featureCounts
mesg "Running featureCounts"
CMD="featureCounts $ARGS"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "featureCounts"
echo ""

## DONE! -----------------------------------------------------------------------
mesg "Counting complete!"
module list
featureCounts -v
echo ""
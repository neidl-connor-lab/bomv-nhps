#!/bin/bash -l

# qsub options
#$ -l h_rt=48:00:00
#$ -l mem_per_core=12G
#$ -pe omp 16
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
ARGS="--runThreadN 16 --runMode genomeGenerate --limitGenomeGenerateRAM 200000000000 --sjdbOverhang 100"
HELP="usage: qsub -P PROJECT -N JOBNAME $(basename "$0") -o ODIR -x FEATURE -g GTF -f FASTA 

options:
  -o output directory 
  -x feature for building transcripts (exon | transcript)
  -g GTF annotation file
  -f reference genome FASTA file
  -h show this message and exit"

# parsing arguments
while getopts ":ho:x:g:f:" opt
do
  case ${opt} in
    o ) ODIR="${OPTARG}"
      ;;
    x ) FEAT="${OPTARG}"
      ;;
    g ) GTF="${OPTARG}"
      ;;
    f ) FASTA="${OPTARG}"
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
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

## check inputs ----------------------------------------------------------------
mesg "STEP 0: CHECK INPUTS"

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
ARGS="$ARGS --genomeDir '$ODIR' --outFileNamePrefix '$ODIR'"

# check feature
if [ -z "$FEAT" ]
then
  err "No feature provided"
elif [ "$FEAT" == "exon" ]
then
  mesg "Using EXON as feature"
elif [ "$FEAT" == "transcript" ]
then
  mesg "Using TRANSCRIPT as feature"
else
  err "Provided feature not recognized: $FEAT"
fi
ARGS="$ARGS --sjdbGTFfeatureExon ${FEAT}"

# check GTF
if [ -z "$GTF" ]
then
  err "GTF file not provided!"
elif [ -f "$GTF" ]
then
  mesg "Valid GTF file: $GTF"
else
  err "Invalid GTF file: $GTF"
fi
ARGS="$ARGS --sjdbGTFfile '$GTF'"

# check genome
if [ -z "$FASTA" ]
then
  err "No genome FASTA file provided"
elif [ -f "$FASTA" ]
then
  mesg "Valid genome FASTA file: $FASTA"
else
  err "Invalid genome FASTA file: $FASTA"
fi
ARGS="$ARGS --genomeFastaFiles '$FASTA'"

# done checking inputs
mesg "Done checking inputs!"
echo ""

## run indexing command --------------------------------------------------------
mesg "STEP 1: BUILD INDEX"

# load STAR
module load star/2.7.1a

# run STAR index
mesg "Running STAR index"
CMD="STAR $ARGS"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "STAR index"
echo ""

## Done! output environment ----------------------------------------------------
module list
echo "STAR $(STAR --version)"

mesg "PIPELINE COMPLETE"
echo ""

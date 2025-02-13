#!/bin/bash
# sbatch options
#SBATCH -A MIP24002
#SBATCH -p skx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o log-%x.out
#SBATCH -t 08:00:00

## setup ----------------------------------------------------------------------
# functions
mesg () {
  echo -e "[MSG] $@"; 
}
err () { 
  echo -e "[ERR] $@"; exit 1; 
}
checkcmd () {
  if [ $? -eq 0 ]
  then 
    mesg "$@ succeeded"
  else 
    err "$@ failed"
  fi
}

# default values and help message
ARGS="--runMode genomeGenerate --runThreadN 1"
FEATURE="transcript" 
HELP="USAGE: sbatch -J JOBNAME $0 [-f FEATURE] -o ODIR -g GENOME -a ANNOTATION

arguments (default):
  -f GTF feature to use to build index ($FEATURE)
  -o index output directory
  -g genome FASTA file, not compressed
  -a annotation GTF file, not compressed
  -h show this message and exit
"

# parsing arguments
while getopts ":hf:o:g:a:" opt
do
  case $opt in
    f ) FEATURE="$OPTARG"
      ;;
    o ) ODIR="$OPTARG"
      ;;
    g ) GENOME="$OPTARG"
      ;;
    a ) ANNOTATION="$OPTARG"
      ;;
    h ) echo "$HELP" && exit 0
      ;;
    \? ) err "Invalid option: $opt \n $HELP"
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

## check inputs ---------------------------------------------------------------
mesg "STEP 1 OF 2: CHECKING INPUTS"

# genome FASTA
if [ -z "$GENOME" ]
then 
  err "No genome FASTA file provided."
elif [ -f "$GENOME" ]
then 
  mesg "Valid genome FASTA file: $GENOME"
else
  err "Invalid genome FASTA file: $GENOME"
fi

# annotation GTF
if [ -z "$ANNOTATION" ]
then
  err "No annotation GTF file provided."
elif [ -f "$ANNOTATION" ]
then 
  mesg "Valid annotation GTF file: $ANNOTATION"
else
  err "Invalid annotation GTF file: $ANNOTATION"
fi

# output directory
if [ -z "$ODIR" ]
then
  err "No index output directory provided."
elif [ -d "$ODIR" ]
then 
  mesg "Valid index output directory: $ODIR"
else 
  mesg "Creating index output directory: $ODIR"
  mkdir -p "$ODIR"
fi

# feature
if [ -z "$FEATURE" ]
then 
  err "No GTF feature provided"
else
  mesg "Using GTF feature: $FEATURE"
fi

# done checking inputs; add args
mesg "Done checking inputs!"
ARGS="$ARGS --genomeFastaFiles '$GENOME' --sjdbGTFfile '$ANNOTATION' --sjdbGTFfeatureExon $FEATURE --genomeDir '$ODIR'"
echo ""

## build and run command ------------------------------------------------------
mesg "STEP 2 OF 2: BUILD STAR INDEX"

# load star
module load biocontainers
module load star

# build command
CMD="STAR $ARGS"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Building STAR index"
echo ""

## print modules --------------------------------------------------------------
module list
mesg "JOB COMPLETE!"
echo ""


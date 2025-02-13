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
HELP="USAGE: sbatch -J JOBNAME $0 -g GENOME -b BASENAME

arguments:
  -g genome FASTA file, not compressed
  -b index output path and basename
  -h show this message and exit
"

# parsing arguments
while getopts ":hg:b:" opt
do
  case $opt in
    g ) GENOME="$OPTARG"
      ;;
    b ) BASENAME="$OPTARG"
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

# output directory and basename
if [ -z "$BASENAME" ]
then
  err "No index basename provided."
else
  mesg "Valid output index basename: $BASENAME"
fi

# make output directory if necessary
if [ ! -d "$(dirname $BASENAME)" ]
then
  mkdir -p "$(dirname $BASENAME)"
fi

# done checking inputs
mesg "Done checking inputs!"
echo ""

## build and run command ------------------------------------------------------
mesg "STEP 2 OF 2: BUILD BOWTIE2 INDEX"

# load bowtie2
module unload xalt
module load biocontainers
module load bowtie2

# build command
CMD="bowtie2-build '$GENOME' '$BASENAME'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "Building bowtie2 index"
echo ""

## print modules --------------------------------------------------------------
module list
mesg "JOB COMPLETE!"
echo ""


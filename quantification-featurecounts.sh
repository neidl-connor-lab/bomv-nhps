#!/bin/bash
# sbatch options
#SBATCH -A MIP24002
#SBATCH -p skx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o log-%x.out
#SBATCH -t 12:00:00

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
FEATURECOUNTS="pipeline/featurecounts.r"
FEATURE="transcript"
META="gene_name"
PAIRED=false
HELP="USAGE: sbatch -J JOBNAME $0 [-f FEATURE -m META -p] -a ANNOTATION -o OFILE BAM1 BAM2...BAMn

arguments (default):
  -f   GTF feature to use for quantification; e.g., exon or transcript ($FEATURE)
  -m   GTF meta-feature to use for quantification; e.g., gene_id or gene_name ($META)
  -p   libraries are paired-end ($PAIRED)
  -a   GTF annotation file (REQUIRED)
  -o   output TSV (REQUIRED)
  -h   show this message and exit
   BAM alignment file paths (REQUIRED)
"

# parsing arguments
while getopts ":hf:m:pa:o:" opt
do
  case $opt in
    f ) FEATURE="$OPTARG"
      ;;
    m ) META="$OPTARG"
      ;;
    p ) PAIRED=true
      ;;
    a ) ANNOTATION="$OPTARG"
      ;;
    o ) OFILE="$OPTARG"
      ;;
    h ) echo "$HELP" && exit 0
      ;;
    \? ) err "Invalid option: $opt \n $HELP"
      ;;
  esac
done
shift $((OPTIND -1))
BAMS="$@"

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

# check featurecounts script
if [ ! -f "$FEATURECOUNTS" ]
then
  err "featureCounts runner script not found: $FEATURECOUNTS"
fi

# annotation GTF
if [ -z "$ANNOTATION" ]
then
  err "No annotation GTF provided."
elif [ -f "$ANNOTATION" ]
then
  mesg "Valid annotation GTF: $ANNOTATION"
else
  err "Invalid annotation GTF: $ANNOTATION"
fi

# output file and directory
ODIR="$(dirname "$OFILE")"
if [ -z "$OFILE" ]
then
  err "No output filename provided."
elif [ -d "$ODIR" ]
then
  mesg "Valid output directory: $ODIR"
else
  mesg "Creating output directory: $ODIR"
  mkdir -p "$ODIR"
fi

# feature and meta-feature
if [ -z "$FEATURE" ]
then
  err "No feature provided."
elif [ -z "$META" ]
then 
  err "No meta-feature provided."
else
  mesg "Using feature: $FEATURE"
  mesg "Using meta-feature: $META"
fi

# paired- or single-end libraries?
if $PAIRED
then
  mesg "Running as paired-end libraries"
else
  mesg "Running as single-end libraries"
fi

# check that at least one filename is provided
if [ -z "$BAMS" ]
then
  err "No alignment BAMs provided."
fi
# loop through all BAMs and check that they exists
# if any BAM doesn't exist, kill the job
for i in $BAMS
do
  if [ ! -f "$i" ]
  then
    err "Invalid BAM file: $i"
  fi
done
# if we're still standing, report that all BAMs exist
mesg "Valid BAM file(s): $BAMS"

# done checking inputs; add args
mesg "Done checking inputs!"
echo ""

## build and run command ------------------------------------------------------
mesg "STEP 2 OF 2: RUN FEATURECOUNTS QUANTIFICATION"

# load R
module use /work/projects/singularity/rstudio/lua
module load Rstats/4.1.3

# build command (paired vs. unpaired)
if $PAIRED
then
  CMD="Rscript $FEATURECOUNTS -f $FEATURE -m $META -p -a '$ANNOTATION' -o '$OFILE' $BAMS"
else
  CMD="Rscript $FEATURECOUNTS -f $FEATURE -m $META -a '$ANNOTATION' -o '$OFILE' $BAMS"
fi
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "featureCounts quantification"
echo ""

## list modules ---------------------------------------------------------------
module list
mesg "JOB COMPLETE!"
echo ""


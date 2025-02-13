#!/bin/bash
# sbatch options
#SBATCH -A MIP24002
#SBATCH -p skx
#SBATCH -N 8
#SBATCH -n 8
#SBATCH -o log-%x.out
#SBATCH -t 02:00:00

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
ARGS="--runThreadN 8 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMunmapped Within"
HELP="USAGE: sbatch -J JOBNAME $0 -i INDEX -o ODIR -s SAMPLE -x R1 [-y R2]

arguments:
  -i path to genome index directory
  -o alignment output directory
  -s sample ID to use as file prefix
  -x R1 FASTQ file
  -y R2 FASTQ file if using paired reads (optional)
  -h show this message and exit
"

# parsing arguments
while getopts ":hi:o:s:x:y:" opt
do
  case $opt in
    i ) INDEX="$OPTARG"
      ;;
    o ) ODIR="$OPTARG"
      ;;
    s ) SAMPLE="$OPTARG"
      ;;
    x ) R1="$OPTARG"
      ;;
    y ) R2="$OPTARG"
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
mesg "STEP 1 OF 3: CHECKING INPUTS"

# set localectl
export LC_ALL=C
export LANG=C

# load modules
module unload xalt
module load biocontainers
module load samtools
module load star

# index directory
if [ -z "$INDEX" ]
then
  err "No genome index directory provided"
elif [ -d "$INDEX" ]
then
  mesg "Valid genome index: $INDEX"
else
  err "Invalid genome index: $INDEX"
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

# R1 file
if [ -z "$R1" ]
then
  err "No R1 FASTQ file provided."
elif [ -f "$R1" ]
then
  mesg "Valid R1 FASTQ file: $R1"
else
  err "Invalid R1 FASTQ file: $R1"
fi

# R2 file (optional) 
if [ -z "$R2" ]
then
  mesg "No R2 file provided; running as unpaired reads."
elif [ -f "$R2" ]
then
  mesg "Valid R2 FASTQ file: $R2"
else
  err "Invalid R2 FASTQ file: $R2"
fi

# sample ID
if [ -z "$SAMPLE" ]
then
  err "No sample ID provided."
else
  mesg "Using sample ID: $SAMPLE"
fi

# done checking inputs; add args
PREFIX="$ODIR/$SAMPLE-"
ARGS="$ARGS --genomeDir '$INDEX' --outFileNamePrefix '$PREFIX' --readFilesIn '$R1'"
# add R2 if present
if [ ! -z "$R2" ]
then
  ARGS="$ARGS '$R2'"
fi
mesg "Done checking inputs!"
echo ""

## build and run command ------------------------------------------------------
mesg "STEP 2 OF 3: RUN STAR ALIGNMENT"

# build command
CMD="STAR $ARGS"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "STAR alignment"
echo ""

# clean up intermediate files
mesg "Cleaning up intermediate files..."
rm "${PREFIX}Log.out"
rm "${PREFIX}Log.progress.out"
rm "${PREFIX}SJ.out.tab"
rm -r "${PREFIX}_STARtmp"

## mark duplicates ------------------------------------------------------------
mesg "STEP 3 OF 3: MARK DUPLICATE READS"

# fixmate
CMD="samtools fixmate --threads 8 -m '${PREFIX}Aligned.out.bam' '${PREFIX}fixmate.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "samtools fixmate"
rm "${PREFIX}Aligned.out.bam"
echo ""

# sort by coordinate
CMD="samtools sort '${PREFIX}fixmate.bam' > '${PREFIX}sorted.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "samtools sort"
rm "${PREFIX}fixmate.bam"
echo ""

# markdup
CMD="samtools markdup --threads 8 '${PREFIX}sorted.bam' '$ODIR/$SAMPLE.bam'"
mesg "CMD: $CMD"
eval "$CMD"
checkcmd "samtools markdup"
rm "${PREFIX}sorted.bam"
echo ""

## list modules ---------------------------------------------------------------
module list
mesg "JOB COMPLETE!"
echo ""


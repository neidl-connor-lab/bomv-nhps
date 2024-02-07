#!/bin/bash -l

# qsub options
#$ -P lasvchal
#$ -l h_rt=48:00:00
#$ -l mem_per_core=12G
#$ -pe omp 16
#$ -j y
#$ -o log-$JOB_NAME.qlog

# job info
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"
echo "=========================================================="
echo ""

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
OVERHANG="100"
OUTDIR="/restricted/projectnb/lasvchal/Jacquelyn/indices/star"
FEAT="exon"
ARGS="--runThreadN 16 --runMode genomeGenerate --limitGenomeGenerateRAM 200000000000"
HELP="usage: qsub -N JOBNAME $(basename "$0") [OPTIONS] GTF FASTA 

positional argument:
  GTF   combined GTF file
  FASTA combined genome FASTA file

optional arguments (default):
  -l read length - 1 (${OVERHANG}) 
  -o output directory (${OUTDIR})
  -f feature for building transcripts (${FEAT})
  -h show this message and exit"

# parsing arguments
while getopts ":hl:o:f:" opt
do
  case ${opt} in
    l ) OVERHANG="${OPTARG}"
      ;;
    o ) OUTDIR="${OPTARG}"
      ;;
    f ) FEAT="${OPTARG}"
      ;;
    h ) echo "${HELP}" && exit 0
      ;;
    \? ) err "Invalid option ${opt}\n${HELP}"
      ;;
  esac
done
shift $((OPTIND -1))

# get GTF and FASTA file(s)
GTF="$1"
FASTA="$2"

# check and add OUTDIR to ARGS
if [ -z "${OUTDIR}" ]
then
  err "No output directory provided"
elif [ -d "${OUTDIR}" ]
then
  mesg "Valid output directory: ${OUTDIR}"
else
  mesg "Creating output directory: ${OUTDIR}"
  mkdir -p "${OUTDIR}"
fi
ARGS="${ARGS} --genomeDir '${OUTDIR}' --outFileNamePrefix '${OUTDIR}'"

# check GTF
if [ -z "${GTF}" ]
then
  err "GTF file not provided!"
elif [ -f "${GTF}" ]
then
  mesg "Valid GTF file: ${GTF}"
else
  err "Invalid GTF file: ${GTF}"
fi
ARGS="${ARGS} --sjdbGTFfile '${GTF}'"

# check genome
if [ -z "${FASTA}" ]
then
  err "No genome FASTA file provided"
elif [ -f "${FASTA}" ]
then
  mesg "Valid genome FASTA file: ${FASTA}"
else
  err "Invalid genome FASTA file: ${FASTA}"
fi
ARGS="${ARGS} --genomeFastaFiles '${FASTA}'"

# check feature and overhang
mesg "Using feature: ${FEAT}"
ARGS="${ARGS} --sjdbGTFfeatureExon ${FEAT}"
mesg "Using overhang: ${OVERHANG}"
ARGS="${ARGS} --sjdbOverhang ${OVERHANG}"
mesg "Done checking inputs"
echo ""

# load STAR and show version
mesg "Loading STAR"
module load star/2.7.1a
echo "STAR $(STAR --version)"

# run STAR index
mesg "Running STAR index"
CMD="STAR ${ARGS}"
mesg "CMD: ${CMD}"
eval "${CMD}"
checkcmd "STAR index"


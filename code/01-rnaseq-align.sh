#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="01-rnaseq-align"
#SBATCH --output=log/%x_%A_%a.log
##################################################################

## RNA-seq Workflow - Alignment #############
##
## This workflow will carry out the following steps:
## 1 - Run QC
## 2 - Trim adapters and low-quality bases
## 3 - Alignment with STAR
## 4 - Indexing
## 5 - IGV formatting
## 6 - Featurecounts

## Get configuration file
if [ -f code/config.txt ]; then
    source code/config.txt
else
    echo "You're missing the config.txt file"
    exit
fi

# get starting time point
echo "Run start"
timestamp

## Job variables
N=${SLURM_ARRAY_TASK_ID}
NUM_CORES=8

## Running workflow
echo "## Starting RNA-seq Workflow ##########"

sed -n ${N}p $SAMPLES | while IFS="," read SAMPLENAME FASTQ1 FASTQ2 REST
do

    if [[ ${FASTQ2} == *"gz"* ]]
    then
      PAIRED=T
    else
      PAIRED=F
    fi
  
    echo "Processing sample = ${SAMPLENAME}"
    echo "Paired-end = ${PAIRED}"
    echo "Strandedness = ${STRANDED}"

    ## Generate output directories
    LIST_DIR=( 
        ${FASTQC_DIR}
        ${TRIMGALORE_DIR}
        ${STAR_DIR}
        ${FEATURECOUNTS_DIR}
    )
    
    for DIR in "${LIST_DIR[@]}"
    do
      [ ! -d "$DIR" ] && mkdir -p "${DIR}"
    done
    
    ## Run pipeline
    #   run_fastqc
    #   run_trimgalore
      run_star
      sam_index
      bam2bw
      run_featurecounts

done

# end time stamp
echo "Run complete"
timestamp
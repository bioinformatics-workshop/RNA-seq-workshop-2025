#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH --time=4:00:00
#SBATCH --job-name="02-rnaseq-diffexpr"
#SBATCH --output=log/%x_%j.log
##################################################################

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

## Generate output directories
[ ! -d "$DEG_DIR" ] && mkdir -p "${DEG_DIR}"

## Running workflow
echo "## Starting Differential abundance analysis ##########"

Rscript --vanilla code/differential_expression.R ${SAMPLE} ${FACTOR}
multiqc_summary

# end time stamp
echo "Run complete"
timestamp
#!/bin/bash -l
#SBATCH --partition=highmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100g
#SBATCH --time=1-12:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="00-build-index"
#SBATCH --output=log/%x_%j.log
##################################################################

## RNA-seq Workflow - Pre-processing #############
##
## This workflow will carry out the following steps:
## 1 - Indexing genome

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

## Job variable
NUM_CORES=8

# Run indexing
run_star_index

# end time stamp
echo "Run complete"
timestamp
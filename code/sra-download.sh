#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00     # 4 hours
#SBATCH --job-name="sra-download"
#SBATCH -p epyc # This is the default partition, you can use any of the following; epyc, short, intel, batch, highmem, gpu
#SBATCH --output=log/%x_%A_%a.log
#################################

# get starting time point
echo "Run start"
timestamp

## Job variables
N=${SLURM_ARRAY_TASK_ID}

SRR_INFO=metadata/metadata_subset_6mo.csv
# SRR_INFO=raw/PRJNA936435.metadata.tmp
RAW_DIR=raw

# current version 3.0.0 4/28/22
module load sratoolkit

sed -n ${N}p $SRR_INFO | cut -f8 | while IFS="," read SAMPLENAME FQ1 FQ2 SRR_ID REST
do
    fastq-dump --split-3 --gzip --outdir ${RAW_DIR} ${SRR_ID}
done

# get ending time point
echo "Run complete"
timestamp

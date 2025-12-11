#!/bin/bash -l
#SBATCH --partition=epyc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12g
#SBATCH --time=0-2:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="get-sra-info"
#SBATCH --output=log/%x_%A.log
##################################################################

## RNA-seq Workflow - Alignment #############
# NCBI tools
module load ncbi_edirect

# Get all SRA runs for a given BioProject
# https://ncbi-hackathons.github.io/EDirectCookbook/

# set parameters
DB="sra"
BIOPROJ="PRJNA936435"
RAW_DIR=raw

# create output directory
[ ! -d "$RAW_DIR" ] && mkdir -p "${RAW_DIR}"

# Search for sample metadata
esearch -db bioproject -query $BIOPROJ | \
elink -target sra | \
efetch -format native -mode xml | \
xtract -pattern EXPERIMENT_PACKAGE \
    -block SAMPLE -element VALUE -block RUN_SET -element PRIMARY_ID > \
    $RAW_DIR/${BIOPROJ}.metadata.tmp

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=6:00:00
#SBATCH --job-name="genome-download"
#SBATCH --output=log/%x_%j.log
#####################################
# Download reference genomes
######################################

# UCSC/GenBank assembly and annotations
# https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips
GENOME=https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
GTF=https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz
# REFGENE_GTF=https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz
OUT_DIR=genome

# download files
# for URL in TRANSCRIPTS GENOME GTF GFF3
# do
#   wget --directory-prefix=${OUT_DIR} ${!URL}
# done

for URL in GENOME GTF REFGENE_GTF
do
  wget --directory-prefix=${OUT_DIR} ${!URL}
done





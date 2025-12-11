#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
library(glue)

# Read in file
data <- read_tsv("raw/PRJNA936435.metadata.tmp", col_names = F, skip = 16)

# Clean up file
data_clean <- data %>%
       select(srr_id = X8, 
              sex = X2, 
              organ = X1,
              strain = X3,
              age = X5) %>%
       mutate(biorep = c(rep(1:4,28)),
              samplename = glue("{organ}_{age}_br{biorep}"),
              factor = glue("{organ}_{age}"),
              fq1 = glue("{srr_id}_1.fastq.gz"),
              fq2 = glue("{srr_id}_2.fastq.gz")) %>%
       select(samplename, fq1, fq2, srr_id, sex, strain, organ, age, biorep, factor)

# Export file
write_csv(data_clean, "metadata/metadata.csv")

data_subset_6mo <- data_clean |>
       filter(age == "6mo")
write_csv(data_subset_6mo, "metadata/metadata_subset_6mo.csv")
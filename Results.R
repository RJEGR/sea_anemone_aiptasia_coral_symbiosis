rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

#dir <- "C:\\wsl.localhos\Ubuntu\home\ricardo\SHORTSTACKS_OUTPUTS"

# path <- "C:/Users/Israel V/Documents/SHORTSTACKS_OUTPUTS/"
path <- "C:/Users/Israel V/Documents/OUTPUTS_SHORTSTACKS_MIRTRACE/"
setwd(path)

res_f <- list.files(path = path, pattern = "Results.txt", full.names = T)

library(tidyverse)

RESULTS <- read_tsv(res_f)

# Q: Number of novel and know miRs annotated

RESULTS %>% dplyr::count(MIRNA)

RES <- RESULTS %>% 
  select(Locus, Name, known_miRNAs, MIRNA) %>%
  mutate(known_miRNAs = ifelse(is.na(known_miRNAs) & MIRNA == "Y", Name, known_miRNAs)) %>%
  # mutate(MIRNA = ifelse(!is.na(known_miRNAs), "Y", MIRNA)) %>%
  drop_na(known_miRNAs) %>%
  mutate(Type = ifelse(grepl("Cluster_", known_miRNAs), "Novel", "known"))

# count number of true novel + homology-identify sRNA locus

RES %>% dplyr::count(Type)

table(substr(RES$Locus, 1,2))


RESULTS %>% dplyr::count(substr(Chrom, 1,8))

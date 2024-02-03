rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

#dir <- "C:\\wsl.localhos\Ubuntu\home\ricardo\SHORTSTACKS_OUTPUTS"

path <- "C:/Users/Israel V/Documents/SHORTSTACKS_OUTPUTS/"

setwd(path)

count_f <- list.files(path = path, pattern = "Counts.txt", full.names = T)

library(tidyverse)

COUNTS <- read_tsv(count_f)

colNames <- names(read_tsv(count_f))

names(COUNTS)[names(COUNTS) %in% colNames] <- gsub("_trimmed", "", colNames)

rowNames <- COUNTS$Name

COUNTS <- COUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(COUNTS) <- rowNames

head(COUNTS)

apply(COUNTS, 1, function(x) sum(x > 0)) %>% table()

prevelancedf = apply(COUNTS, 1, function(x) sum(x > 0))

data.frame(Prevalence = prevelancedf, 
           TotalAbundance = rowSums(COUNTS)) %>% # mean_se
  as_tibble(rownames = "Name") %>%
  arrange(desc(TotalAbundance)) -> prevelancedf

prevelancedf %>% 
  count(Prevalence) %>% 
  ggplot(aes(Prevalence, n)) + geom_col() +
  theme_classic(base_family = "GillSans") + 
  scale_y_continuous("Number of sRNAs", labels = scales::comma) +
  scale_x_continuous(breaks = 1:12) -> ps

                     
# COR heatmap
# USING VST FOR HEATMAP

head(DATA <- assay(vst)) # NEEDS DESEQ2

sample_cor = cor(DATA, method='pearson', use='pairwise.complete.obs')

h <- heatmap(sample_cor, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Rowv)

hc_order <- hc_samples$labels[h$rowInd]
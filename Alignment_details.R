# PLOT alignment_details.tsv

# alignment details ----

# mapping_type
# U: Uniquely mapped (not a multimapper).
# P: Multimapper placed using the method set by option --mmap.
# R: Multimapper placed at random.
# H: Very highly multi-mapped read (>=50 hits).
# N: Unmapped reads.

rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

library(tidyverse)

recode_to <- c(`U` = "Uniquely ", `P`= "Multimapper (mmap)",`R` = "Random mmap", `H` = "Highly mmap (>=50 hits) ", `N` = "Unmapped")

# path <- "C:/Users/Israel V/Documents/SHORTSTACKS_OUTPUTS/"

path <- "C:/Users/Israel V/Documents/OUTPUTS_SHORTSTACKS_MIRTRACE/"


setwd(path)

f <- list.files(path = path, pattern = "alignment_details.tsv", full.names = T)

alignment_details <- read_tsv(f)

readfile <- gsub("_trimmed.mirna.unknown.bam", "", basename(alignment_details$readfile))

alignment_details$LIBRARY_ID <- readfile

mtd_f <- list.files(path = path, pattern = "METADATA", full.names = T)

MTD <- read_tsv(mtd_f)


alignment_details <- alignment_details %>%
  left_join(MTD) %>% select(LIBRARY_ID, Design, mapping_type,read_length,count)

alignment_details %>% group_by(LIBRARY_ID) %>% tally(count) %>% view()# <- concordantly w/ initial proccessed reads


# mapping_type read_length   count
# ES
# recode_to <- c(`U` = "Único ", `P`= "Múltiple (mmap)",`R` = "mmap aleatorio", `H` = "mmap (>=50 hits)", `N` = "Sin mapeo")

# simplificar a mapeos unicos, multiples y sin mapeo
recode_to <- c(`U` = "Único ", `P`= "Múltiple",`H` = "Múltiple", `R` = "Múltiple", `N` = "Sin mapeo")

ylab <- "% Lecturas"
xlab <- "Longitud (Nucleótidos)"

read_lengthL <- c("<21", "21", "22", "23", "24", ">24")

alignment_details %>%
  dplyr::mutate(mapping_type = dplyr::recode_factor(mapping_type, !!!recode_to)) %>%
  group_by(Design, mapping_type, read_length) %>% 
  summarise(n = sum(count)) %>%
  group_by(read_length) %>%
  mutate(frac = n/sum(n)) %>%
  dplyr::mutate(read_length = factor(read_length, levels = read_lengthL)) %>%
  ggplot(aes(x = read_length, y = frac, fill = mapping_type)) +
  geom_col(width = 0.85) + # position="dodge"
  scale_y_continuous(ylab, labels = scales::percent) +
  labs(x = xlab) +
  guides(fill = guide_legend(title = "Mapeo")) + #  nrow = 5
  see::scale_fill_flat(reverse = T) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(strip.background = element_rect(fill = 'white', color = 'white'),
        legend.position = 'right',
        panel.border = element_blank(),
        panel.grid.minor = element_blank()) -> ps

# https://easystats.github.io/see/articles/seecolorscales.html#overview-of-palette-colors

ggsave(ps, filename = 'ALIGNMENT_DETAILS_ES.png', path = path, width = 4, height = 2, device = png)



# CALL FOR VARIANT AND PRECISION IDENTIFICATION:
which_vars <- c(c(as.character(21:24)))

f <- list.files(path = path, pattern = "Results.txt", full.names = T)


Treads <- read_tsv(f) %>% select_at(all_of(which_vars)) %>% rowSums()

read_tsv(f) %>%
  mutate(PRECISION = Treads/Reads) %>%
  ggplot(aes(x = MIRNA, y = PRECISION)) +
  geom_boxplot()

which_vars <- c("Short", c(as.character(21:30)),  "Long")

sum(read_tsv(f)$Reads == Treads)

read_tsv(f) %>% 
  mutate(PRECISION = Treads/Reads) %>%
  select_at(all_of(c(which_vars, 'DistinctSequences', 'MajorRNAReads', 'Reads', "FracTop", "PRECISION"))) %>%
  rstatix::cor_mat() %>%
  rstatix::cor_reorder() %>%
  rstatix::pull_lower_triangle() %>%
  rstatix::cor_plot(label = TRUE)

which_vars2 <- c('DicerCall','DistinctSequences', 'MajorRNAReads', 'Reads')

# ???

read_tsv(f) %>% select_at(all_of(c(which_vars, which_vars2))) %>%
  mutate(PRECISION = (DistinctSequences+MajorRNAReads)/Reads) %>%
  ggplot(aes(x = DicerCall, y = PRECISION)) +
  geom_boxplot()



read_tsv(f) %>% 
  mutate(PRECISION = Treads/Reads) %>%
  ggplot(aes(PRECISION, FracTop)) + geom_point()

# FracTop: Fraction of Reads aligned to the top genomic strand.
# Strand: Inferred strandednes of the locus, inferred from FracTop and the --strand_cutoff setting
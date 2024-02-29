rm(list = ls())

if(!is.null(dev.list())) dev.off()

options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)

path <- "D:/"

setwd(path)

count_f <- list.files(path = path, pattern = "Bre_All_expressed_miRNAs.xlsx", full.names = T)

mtd_f <- list.files(path = path, pattern = "METADATA", full.names = T)

library(tidyverse)

MTD <- read_tsv(mtd_f)

COUNTS <- readxl::read_xlsx(count_f)

rowNames <- COUNTS$`#miRNA`

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
  scale_x_continuous(breaks = 1:16)


# 1) Filter data by removing low-abundance genes ----

by_count <- 1; by_freq <- 2

keep <- rowSums(COUNTS > by_count) >= by_freq

sum(keep) # 215 sRNAs from 

nrow(COUNTS <- COUNTS[keep,])

COUNTS <- round(COUNTS)

# 2) COR heatmap -----
# USING VST FOR HEATMAP


# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)

head(DATA <- DESeq2::varianceStabilizingTransformation(COUNTS+1)) # NEEDS DESEQ2

sample_cor = cor(DATA, method='pearson', use='pairwise.complete.obs')

h <- heatmap(sample_cor, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Rowv)

hc_order <- hc_samples$labels[h$rowInd]

#

sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  mutate(group = sapply(strsplit(LIBRARY_ID, "_"), `[`, 1)) %>%
  # left_join(MTD) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order))) -> sample_cor_long

library(ggh4x)

P <- sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) +  
  geom_tile(linewidth = 0.2) +
  scale_x_discrete(position = 'bottom') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) +
  # ggh4x::scale_x_dendrogram(hclust = hc_samples, position = "top", labels = NULL) +
  theme_bw(base_size = 7, base_family = "GillSans") +
  ggsci::scale_fill_material("indigo") +
  labs(x = '', y = '') +
  theme(
    legend.position = "bottom",
    # axis.text.x = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
    strip.background = element_rect(fill = 'grey89', color = 'white'),
    panel.border = element_blank(),
    plot.title = element_text(hjust = 0),
    plot.caption = element_text(hjust = 0),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()) 

P <- P + guides(
  fill = guide_colorbar(barwidth = unit(1.5, "in"),
                        barheight = unit(0.05, "in"), label.position = "bottom",
                        alignd = 0.5,
                        ticks.colour = "black", ticks.linewidth = 0.5,
                        frame.colour = "black", frame.linewidth = 0.5,
                        label.theme = element_text(family = "GillSans", size = 7)))


# TOP

TOPDF <- sample_cor_long %>%
  distinct(LIBRARY_ID, group) %>%
  # dplyr::mutate(hpf = dplyr::recode_factor(hpf, !!!recode_to)) %>%
  # mutate(label = ifelse(pH %in% "Low", "*", "")) %>%
  mutate(y = 1)


topplot <- TOPDF %>%
  ggplot(aes(y = y, x = LIBRARY_ID, color = as.factor(group))) +
  geom_point(shape = 15, size = 2) +
  #geom_text(aes(label = label),  vjust = -0.7, hjust = 0.5, size = 1.5, family =  "GillSans", color = "#d73027") +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  # ggh4x::guide_dendro()
  # guides(x.sec = guide_axis_manual(labels = hc_order, label_size = 3.5)) +
  guides(color = guide_legend(title = "", ncol = 4, )) +
  theme_bw(base_family = "GillSans", base_size = 7) +
  # see::scale_color_pizza(name = "", reverse = T) +
  # scale_color_manual("", values = c("#DADADA", "#D4DBC2")) + # "#4575b4", "#d73027"
  theme(legend.position = 'top',
        panel.border = element_blank(), 
        plot.background = element_rect(fill='transparent', color = 'transparent'),
        plot.margin = unit(c(0,0,0,0), "pt"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank()) 

library(patchwork)

psave <- topplot/ plot_spacer() /P + plot_layout(heights = c(0.6, -0.5, 5))

psave

ggsave(psave, filename = 'SAMPLE_HEATMAP_MIR.png', path = path, width = 4, height = 4, device = png, dpi = 300)


# PCA
vst <- DESeq2::vst(COUNTS) # vst if cols > 10

PCA = prcomp(t(vst), center = T, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

psave <- PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  left_join(MTD) %>%
  ggplot(., aes(PC1, PC2)) +
  # coord_fixed(ratio = sd_ratio) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  ggforce::geom_mark_ellipse(aes(group = as.factor(Design), label = as.factor(Design)),
                             fill = 'grey', colour = 'white') +
  geom_point(size = 7, alpha = 0.7, aes(color = Design)) +
  # geom_text( family = "GillSans", mapping = aes(label = LIBRARY_ID), size = 2.5) +
  labs(caption = '', color = "") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

ggsave(psave, filename = 'SAMPLE_PCA_MIR.png', path = path, width = 7, height = 7, device = png, dpi = 300)


# https://github.com/RJEGR/Small-RNASeq-data-analysis/blob/master/E_expression_analysis/COUNT_TRANSFORMATION.R
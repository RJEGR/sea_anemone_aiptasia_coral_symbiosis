# Correlation analysis

Calculate Pearson correlation coeficient and then visualize sample-to-sample heatmap. A correlation heatmap is a graphical representation that displays the correlation coefficients between pairs of variables.

# Author

Dr. Ricardo Gomez-Reyes



# Set the session

Clean the canvas r memory

```R
rm(list = ls())

if(!is.null(dev.list())) dev.off()

```

Set the work directory

```R
path <- "wheremydatafind"

setwd(path)
```

List the absolute path were files find

```R
count_f <- list.files(path = path, pattern = "Bre_All_expressed_miRNAs.xlsx", full.names = T)

mtd_f <- list.files(path = path, pattern = "METADATA", full.names = T)

```

Load cran package

```R
library(tidyverse)
library(DESeq2)
library(ggh4x)
library(ggsci)
library(patchwork)


# (Optional) If extra options parameter needed use 
options(stringsAsFactors = FALSE, readr.show_col_types = FALSE)
```

Load data with different formats

```R
MTD <- readr::read_tsv(mtd_f)
COUNTS <- readxl::read_xlsx(count_f)
```

Convert dataframes to matrix prior to perform correlation analysis (format is requiered)

```R
rowNames <- COUNTS$`#miRNA`

COUNTS <- COUNTS %>% select_if(is.double) %>% as(., "matrix")

rownames(COUNTS) <- rowNames

head(COUNTS)
```

Filter data by removing low-abundance genes. Usually by count > 1 and sample frequency >= 2.

```R
by_count <- 1; by_freq <- 2

keep <- rowSums(COUNTS > by_count) >= by_freq

nrow(COUNTS <- COUNTS[keep,])

COUNTS <- round(COUNTS)
```

# Heatmap 

Calculating correlation coeficient

```R

# when samples > 10, use DESeq2::vst(COUNTS) that provides much faster estimation of the dispersion trend

head(DATA <- DESeq2::varianceStabilizingTransformation(COUNTS)) # use pseudo-count as COUNTS+1 if warning

sample_cor = cor(DATA, method='pearson', use='pairwise.complete.obs')
```

Reordering samples. Calculate left and top dendogram

```R
h <- heatmap(sample_cor, col = cm.colors(12), keep.dendro = T)

hc_samples <- as.hclust(h$Rowv)

hc_order <- hc_samples$labels[h$rowInd]
```

Prepare data for nice plotting

```R
sample_cor %>% 
  as_tibble(rownames = 'LIBRARY_ID') %>%
  pivot_longer(cols = colnames(sample_cor), values_to = 'cor') %>%
  mutate(group = sapply(strsplit(LIBRARY_ID, "_"), `[`, 1)) %>%
  # left_join(MTD) %>% 
  mutate(LIBRARY_ID = factor(LIBRARY_ID, levels = rev(hc_order))) -> sample_cor_long
```

Plot several plots

1) Core heatmap

```R
P <- sample_cor_long %>%
  ggplot(aes(x = LIBRARY_ID, y = name, fill = cor)) +  
  geom_tile(linewidth = 0.2) +
  scale_x_discrete(position = 'bottom') +
  ggh4x::scale_y_dendrogram(hclust = hc_samples, position = "left", labels = NULL) +
  guides(y.sec = guide_axis_manual(labels = hc_order, label_size = 5, label_family = "GillSans")) +
  theme_bw() +
  ggsci::scale_fill_material("indigo") +
  labs(x = '', y = '') +
  theme(
    legend.position = "bottom",
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
```

2. Top dendogram

```R
TOPDF <- sample_cor_long %>%
  distinct(LIBRARY_ID, group) %>%
  mutate(y = 1)


topplot <- TOPDF %>%
  ggplot(aes(y = y, x = LIBRARY_ID, color = as.factor(group))) +
  geom_point(shape = 15, size = 2) +
  ggh4x::scale_x_dendrogram(hclust = hc_samples, position = 'top', labels = NULL) +
  guides(color = guide_legend(title = "", ncol = 4, )) +
  theme_bw() +
  theme(legend.position = 'top',
        panel.border = element_blank(), 
        plot.background = element_rect(fill='transparent', color = 'transparent'),
        plot.margin = unit(c(0,0,0,0), "pt"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank()) 
```

3. paste plots and save

```R
psave <- topplot/ plot_spacer() /P + plot_layout(heights = c(0.6, -0.5, 5))

psave

ggsave(psave, filename = 'SAMPLE_HEATMAP.png', path = path, width = 4, height = 4)

```

# Principal component analysis

Quick visualization for dimensionality reduction.

```R
PCA = prcomp(t(DATA), center = T, scale. = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

sd_ratio <- sqrt(percentVar[2] / percentVar[1])

PCAdf <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2])

psave <- PCAdf %>%
  mutate(LIBRARY_ID = rownames(.)) %>%
  # left_join(MTD) %>%
  mutate(group = sapply(strsplit(LIBRARY_ID, "_"), `[`, 1)) %>%
  ggplot(., aes(PC1, PC2)) +
  coord_fixed(ratio = sd_ratio) +
  geom_abline(slope = 0, intercept = 0, linetype="dashed", alpha=0.5) +
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.5) +
  ggforce::geom_mark_ellipse(aes(group = as.factor(group), label = as.factor(group)),
                             fill = 'grey', colour = 'white') +
  geom_point(size = 7, alpha = 0.7, aes(color = group)) +
  # geom_text(mapping = aes(label = LIBRARY_ID), size = 2.5) +
  labs(caption = '', color = "") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme_bw(base_family = "GillSans", base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')

ggsave(psave, filename = 'SAMPLE_PCA_MIR.png', path = path, width = 7, height = 7)
```



![figure-pca](/Users/cigom/Documents/GitHub/sea_anemone_aiptasia_coral_symbiosis/Figures/figure-pca.jpeg)
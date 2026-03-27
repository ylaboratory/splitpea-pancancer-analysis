library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

cancer_tissue_splice <- read.csv("/results/cancer_alt_splice.csv")
colnames(cancer_tissue_splice) <- c("cancer_tissue_splicemb", "Cancer.Types")
cancer_tissue_splice$Cancer.Types <- gsub("_", "+", cancer_tissue_splice$Cancer.Types)
cancer_tissue_splice$Tissue.Types <- sub(".*\\+", "", cancer_tissue_splice$Cancer.Types)

cancers <- unique(cancer_tissue_splice$Cancer.Types)

cancer_tissue_splice$Cancer.Main <- sub("\\+.*", "", cancer_tissue_splice$Cancer.Types)

median <- as_tibble(cancer_tissue_splice) %>% 
  group_by(Cancer.Types) %>% 
  summarise(median(cancer_tissue_splicemb))
median <- median[match(cancers, median$Cancer.Types), ]
colnames(median) <- c("Cancer.Types", "median1")

tt <- as.data.frame(table(cancer_tissue_splice[2]))
tt <- tt[match(cancers, tt$Cancer.Types), ]

tcga.cohort.med = merge(tt, median, by.x = 'Cancer.Types', by.y = 'Cancer.Types')
tcga.cohort.med <- tcga.cohort.med[match(cancers, tcga.cohort.med$Cancer.Types), ]
colnames(tcga.cohort.med) = c('Cohort', 'Cohort_Size', 'Median_Mutations')
tcga.cohort.med <- tcga.cohort.med[order(tcga.cohort.med$Median_Mutations), ]
cancers <- unique(tcga.cohort.med$Cohort)

tcga.cohort = split(cancer_tissue_splice, as.factor(cancer_tissue_splice$Cancer.Types))
tcga.cohort <- tcga.cohort[cancers]

plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
  x = tcga.cohort[[i]]
  
  x = x[order(x$cancer_tissue_splicemb, decreasing = TRUE), ]
  
  result = data.table::data.table(
    rev(seq(i-1, i, length.out = nrow(x))),
    x$cancer_tissue_splicemb,
    x$Cancer.Types,
    x$Tissue.Types
  )
  
  return(result)
})

names(plot.dat) = names(tcga.cohort)
plot.dat <- data.table::rbindlist(plot.dat)
colnames(plot.dat) <- c("x", "y", "Cancer.Types", "Tissue.Types")

plot.dat$Cancer.Main <- sub("\\+.*", "", plot.dat$Cancer.Types)
plot.dat$Color <- factor(plot.dat$Tissue.Types, levels = names(tissue_colors))

tcga.cohort.med$Median_Mutations_log10 <- log10(tcga.cohort.med$Median_Mutations)

ggplot(plot.dat, aes(x = x, y = y, color = Color)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_segment(data = tcga.cohort.med, aes(x = seq_along(Cohort) - 1, xend = seq_along(Cohort), y = Median_Mutations, yend = Median_Mutations), color = 'black', size = 0.8) +
  scale_y_continuous(limits = c(0, max(plot.dat$y, na.rm = TRUE)), breaks = seq(0, max(plot.dat$y, na.rm = TRUE), 10000)) +
  scale_x_continuous(breaks = seq(0.5, length(cancers) - 0.5, 1), labels = cancers) +
  labs(
    title = "Skipped Exon Events Across Various Cancer Types", 
    x = "Cancer Types", 
    y = "Number of Skipped Exon Events"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5), # Center align and adjust title size
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.y = element_text(size = 10), # Smaller font size for y-axis label
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 0.3
  ) +
  scale_color_manual(values = tissue_colors) 





cancer_tissue_splice_corrected <- read.csv("/results/cancer_tissue_alt_splice.csv")
colnames(cancer_tissue_splice_corrected) <- c("cancer_tissue_splice_correctedmb", "Cancer.Types")
cancer_tissue_splice_corrected$Cancer.Types <- gsub("_", "+", cancer_tissue_splice_corrected$Cancer.Types)
cancer_tissue_splice_corrected$Tissue.Types <- sub(".*\\+", "", cancer_tissue_splice_corrected$Cancer.Types)

cancers <- unique(cancer_tissue_splice_corrected$Cancer.Types)

cancer_tissue_splice_corrected$Cancer.Main <- sub("\\+.*", "", cancer_tissue_splice_corrected$Cancer.Types)

cancer_colors <- c("OV" = "#1f78b4", "MESO" = "#7fcdbb", "TGCT" = "#fdb462",
                   "UCS" = "#fee391", "ACC" = "#41ab5d", "SKCM" = "#238b45",
                   "KICH" = "#ef3b2c", "BRCA" = "#fb6a4a", "STAD" = "#6a51a3",
                   "LIHC" = "#807dba", "THCA" = "#bcbddc", "GBM" = "#dadaeb",
                   "CESC" = "#ce1256", "LUAD" = "#f768a1", "READ" = "#969696",
                   "UCEC" = "#bdbdbd", "BLCA" = "#d9d9d9", "LUSC" = "#e5e5e5",
                   "KIRP" = "#41b6c4", "COAD" = "#a1dab4", "LGG" = "#1f78b4",
                   "ESCA" = "#7fcdbb", "PRAD" = "#fdb462", "PAAD" = "#fee391",
                   "LAML" = "#41ab5d")


tissue_colors <- c("ovary" = "#1f4f9c", "lung" = "#379f7c", "testis" = "#d95f0e",
                  "uterus" = "#fdae6b", "adrenal gland" = "#2c7f2b", "skin" = "#006d2c",
                  "kidney" = "#cb181d", "breast" = "#de2d26", "stomach" = "#54278f",
                  "cervix uteri" = "#6a51a3", "liver" = "#9e9ac8", "thyroid" = "#bcbddc",
                  "blood" = "#980043", "colon" = "#ae017e", "bladder" = "#fdb462",
                  "esophagus" = "#f768a1", "prostate" = "#fee391", "brain" = "#02818a",
                  "pancrease" = "#969696", "small intestine" = "9e1ba8"
                )

median <- as_tibble(cancer_tissue_splice_corrected) %>% 
  group_by(Cancer.Types) %>% 
  summarise(median(cancer_tissue_splice_correctedmb))
median <- median[match(cancers, median$Cancer.Types), ]
colnames(median) <- c("Cancer.Types", "median1")

tt <- as.data.frame(table(cancer_tissue_splice_corrected[2]))
tt <- tt[match(cancers, tt$Cancer.Types), ]

tcga.cohort.med = merge(tt, median, by.x = 'Cancer.Types', by.y = 'Cancer.Types')
tcga.cohort.med <- tcga.cohort.med[match(cancers, tcga.cohort.med$Cancer.Types), ]
colnames(tcga.cohort.med) = c('Cohort', 'Cohort_Size', 'Median_Mutations')
tcga.cohort.med <- tcga.cohort.med[order(tcga.cohort.med$Median_Mutations), ]
cancers <- unique(tcga.cohort.med$Cohort)

tcga.cohort = split(cancer_tissue_splice_corrected, as.factor(cancer_tissue_splice_corrected$Cancer.Types))
tcga.cohort <- tcga.cohort[cancers]

plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
  x = tcga.cohort[[i]]
  
  x = x[order(x$cancer_tissue_splice_correctedmb, decreasing = TRUE), ]
  
  result = data.table::data.table(
    rev(seq(i-1, i, length.out = nrow(x))),
    x$cancer_tissue_splice_correctedmb,
    x$Cancer.Types,
    x$Tissue.Types
  )
  
  return(result)
})

names(plot.dat) = names(tcga.cohort)
plot.dat <- data.table::rbindlist(plot.dat)
colnames(plot.dat) <- c("x", "y", "Cancer.Types", "Tissue.Types")

plot.dat$Cancer.Main <- sub("\\+.*", "", plot.dat$Cancer.Types)
plot.dat$Color <- factor(plot.dat$Tissue.Types, levels = names(tissue_colors))

tcga.cohort.med$Median_Mutations_log10 <- log10(tcga.cohort.med$Median_Mutations)

# Plot
ggplot(plot.dat, aes(x = x, y = y, color = Color)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_segment(data = tcga.cohort.med, aes(x = seq_along(Cohort) - 1, xend = seq_along(Cohort), y = Median_Mutations, yend = Median_Mutations), color = 'black', size = 0.8) +
  scale_y_continuous(limits = c(0, max(plot.dat$y, na.rm = TRUE)), breaks = seq(0, max(plot.dat$y, na.rm = TRUE), 1000)) +
  scale_x_continuous(breaks = seq(0.5, length(cancers) - 0.5, 1), labels = cancers) +
  labs(
    title = "Alternative Splicing Events Across Various Cancer & Tissue Types", 
    x = "Cancer & Tissue Types", 
    y = "Number of Alternative Splicing Events \n Corrected For Normal Tissue"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5), 
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.y = element_text(size = 10), 
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 0.3
  ) +
  scale_color_manual(values = tissue_colors) 



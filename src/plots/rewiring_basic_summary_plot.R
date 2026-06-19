library(ggplot2)
library(dplyr)
library(gridExtra)  
library(grid)      

# boxplot for total edge gains and losses across cancer types
data <- read.table('/results/tcga_splitpea_network.largest_ccs.txt', sep = '\t', header = TRUE)

data$cancer <- sapply(strsplit(as.character(data$cancer), '_'), `[`, 1)

gains <- data %>% filter(direction == 'positive')
losses <- data %>% filter(direction == 'negative')

gains$type <- 'Gains'
losses$type <- 'Losses'

combined_data <- rbind(gains, losses)

combined_median_sorted <- combined_data %>%
  group_by(cancer, type) %>%
  summarize(median_lcc = median(lcc_edges, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(cancer) %>%
  summarize(mean_median_lcc = mean(median_lcc, na.rm = TRUE)) %>%
  arrange(mean_median_lcc)

combined_data$cancer <- factor(combined_data$cancer, levels = combined_median_sorted$cancer)

ggplot(combined_data, aes(y = cancer, x = lcc_edges, fill = type)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(values = c('#1f78b4', '#e31a1c')) +
  labs(title = 'Total Edge Gains and Losses Across Cancer Types',
       y = 'Cancer Type',
       x = 'Number of Edges') +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = c(0.8, 0.2),  
    legend.title = element_text(face = "bold")
  )



chaos_lcc <- data %>% filter(direction == "chaos")
chaos_lcc$cancer <- gsub("_.*", "", chaos_lcc$cancer)

median_lcc <- chaos_lcc[, .(median_lcc_edges = median(lcc_edges)), by = cancer]

chaos_lcc$cancer <- factor(chaos_lcc$cancer, levels = median_lcc[order(median_lcc_edges), cancer])

ggplot(chaos_lcc, aes(x = cancer, y = lcc_edges)) +
  geom_boxplot() +
  labs(x = "Cancer Type", y = "total chaos gain", title = "Boxplot of LCC Edges in Chaos Direction by Cancer Type") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 



# scatter plot between the number of significant SE events and the number of rewired edges 
gained_edges <- data %>% filter(direction == 'positive')
lost_edges <- data %>% filter(direction == 'negative')

palette_colors <- scales::hue_pal()(length(unique(data$cancer_type)))

plot_gained <- ggplot(gained_edges, aes(x = alternative_se_sig_pd, y = lcc_edges, color = cancer_type)) +
  geom_point(alpha = 0.5, size = 1.5) +  
  scale_color_manual(values = palette_colors, , name = "Cancer Types") +
  labs(y = "Number of gained edges") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.x = element_blank(),  
    axis.title.y = element_text(size = 10), 
    axis.text = element_text(size = 10),  
    legend.position = c(0.1, 0.65),  
    legend.title = element_text(size = 10),  
    legend.text = element_text(size = 8),  
    legend.background = element_blank(), 
    legend.key.size = unit(0.5, "cm"),  
    plot.margin = margin(t = 10, r = 10, b = 0, l = 10)
  ) +
  guides(color = guide_legend(ncol = 2)) 

plot_lost <- ggplot(lost_edges, aes(x = alternative_se_sig_pd, y = lcc_edges, color = cancer_type)) +
  geom_point(alpha = 0.5, size = 1.5) + 
  scale_color_manual(values = palette_colors) +
  labs(x = "Number of significant alternative splicing events", y = "Number of lost edges") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 10),  
    axis.title.y = element_text(size = 10),  
    axis.text = element_text(size = 10),  
    legend.position = "none",  
    plot.margin = margin(t = 0, r = 10, b = 10, l = 10)
  )

grid.arrange(
  plot_gained, plot_lost, nrow = 2, 
  top = textGrob("Number of Significant Alternative Splicing Events Vs Number of Edges", 
                 gp = gpar(fontsize = 14))
)



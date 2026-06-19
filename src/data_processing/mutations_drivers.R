library(ggplot2)
library(dplyr)
library(rstatix)

boxplot_node_degree <- read.csv('/data/driver_non_mutation_df.csv')

b_df_fi <- boxplot_node_degree %>%
  select(GENE_SYMBOL, Participant, degree, label) %>%
  distinct()

# Boxplot: degree distribution per label 

sorted_labels <- b_df_fi %>%
  group_by(label) %>%
  summarise(mean_degree = mean(degree), .groups = "drop") %>%
  arrange(mean_degree) %>%
  pull(label)

p_box <- ggplot(b_df_fi, aes(x = factor(label, levels = sorted_labels), y = degree)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Degree by Label", x = "Label", y = "Degree") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12)
  ) +
  coord_cartesian(ylim = c(quantile(b_df_fi$degree, 0.05), quantile(b_df_fi$degree, 0.95)))

print(p_box)

# Bar plot: mean degree +/- SE per label 

barplot_data <- b_df_fi %>%
  group_by(label) %>%
  summarise(
    mean_degree = mean(degree, na.rm = TRUE),
    se_degree = sd(degree, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(label = gsub(" ", "\n", label))

sorted_labels_barplot <- barplot_data %>%
  arrange(mean_degree) %>%
  pull(label)

p_bar <- ggplot(barplot_data, aes(x = factor(label, levels = sorted_labels_barplot), y = mean_degree)) +
  geom_bar(stat = "identity", fill = "#d3d3d3", width = 0.5) +
  geom_errorbar(aes(ymin = mean_degree - se_degree, ymax = mean_degree + se_degree), width = 0.2) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_minimal() +
  labs(title = "Mean Degree by Label", x = "Label", y = "Mean Degree") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  )

print(p_bar)

# Wilcoxon tests 

b_df_fi_collapse <- b_df_fi %>%
  mutate(group = factor(
    ifelse(label %in% c("mutated driver", "non-mutated driver"), "driver", "non-driver"),
    levels = c("driver", "non-driver")
  ))

driver_vs_nondriver_wilcox <- wilcox_test(b_df_fi_collapse, degree ~ group)
driver_vs_nondriver_wilcox

driver_only <- b_df_fi %>%
  filter(label %in% c("mutated driver", "non-mutated driver")) %>%
  mutate(label = factor(label, levels = c("mutated driver", "non-mutated driver")))

mutdriver_vs_nonmutdriver_wilcox <- wilcox_test(driver_only, degree ~ label)
mutdriver_vs_nonmutdriver_wilcox

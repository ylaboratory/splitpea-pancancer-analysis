library(ggplot2)
library(dplyr)

data <- read.csv("/data/edge_stats.csv", stringsAsFactors = TRUE)

filtered_data <- data %>%
  filter((Positive_Count + Negative_Count) / Total_Samples >= 0.5)

ggplot(data, aes(x = Stability)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "lightblue", color = "black") +
  facet_wrap(~ cancer, ncol = 4) +
  labs(title = "Consistency Density Plots by Cancer Type", x = "Consistency", y = "Density") +
  theme_minimal()

ggplot(filtered_data, aes(x = Stability)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "lightblue", color = "black") +
  facet_wrap(~ cancer, ncol = 4) +
  labs(title = "Consistency For Edges in >= 0.5 samples in their cancer", x = "Consistency", y = "Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

quantile_breaks <- quantile(filtered_data$Stability, probs = seq(0, 1, length.out = 7), na.rm = TRUE)

ggplot(filtered_data, aes(x = Stability)) +
  geom_histogram(aes(y = ..density..), breaks = quantile_breaks, fill = "lightblue", color = "black") +
  facet_wrap(~ cancer, ncol = 4) +
  labs(title = "Consistency Density Plots by Cancer Type", x = "Consistency", y = "Density") +
  theme_minimal()

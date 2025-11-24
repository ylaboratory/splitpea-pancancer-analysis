library(ggplot2)
library(dplyr)
library(scales) 


res_melted <- read.csv("data/driver_splitpea_info.csv")
res_melted <- res_melted %>%
  mutate(driver_status = case_when(
    driver_status == "1_driver" ~ "1 driver",
    driver_status == "2_drivers" ~ "2 drivers",
    driver_status == "no_drivers_connected_to_driver" ~ "has driver neighbor",
    driver_status == "no_drivers_no_connection_to_driver" ~ "no driver neighbor",
    TRUE ~ driver_status 
  ))


data_order <- read.table('data/tcga_splitpea_network.largest_ccs.txt',
                           sep = '\t', header = TRUE)
data_order$cancer <- sapply(strsplit(as.character(data_order$cancer), '_'), `[`, 1)

gains_order <- data_order %>% filter(direction == 'positive')
losses_order <- data_order %>% filter(direction == 'negative')
gains_order$type <- 'Gains'
losses_order$type <- 'Losses'
combined_order <- rbind(gains_order, losses_order)

combined_median_sorted <- combined_order %>%
  group_by(cancer, type) %>%
  summarize(median_lcc = median(lcc_edges, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(cancer) %>%
  summarize(mean_median_lcc = mean(median_lcc, na.rm = TRUE)) %>%
  arrange(mean_median_lcc)

res_melted$cancer <- factor(res_melted$cancer, levels = combined_median_sorted$cancer)


df_summary <- res_melted %>%
  group_by(cancer, driver_status) %>%
  summarize(mean_value = mean(value, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(cancer) %>%
  mutate(total = sum(mean_value),
         proportion = mean_value / total) %>%
  ungroup()

df_summary <- res_melted %>%
  group_by(cancer, driver_status) %>%
  summarize(total_value = sum(value, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(cancer) %>%
  mutate(total = sum(total_value),
         proportion = total_value / total) %>%
  ungroup()

p <- ggplot(df_summary, aes(x = cancer, y = proportion, fill = driver_status)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Proportion of Edges by Cancer and Driver Status",
    x = "Cancer Type",
    y = "Proportion"
  )  +
  coord_flip()

print(p)

#ggsave("/plots/stacked_bar_driver.pdf", plot = p)
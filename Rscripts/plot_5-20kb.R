library(tidyverse)
library(patchwork)

deseq_data <- read.table(gzfile("results/RNAseq_log2_foldchange.txt.gz"), header = TRUE)

#filter conditions
filtered_data <- deseq_data %>%
  filter(padj < 0.01, abs(log2FoldChange) > 1)

annotations <- "results/annotated_peaks.txt"
annotation_data <- read_delim(file = annotations, delim = "\t") %>%
  as_tibble() %>%
  select(`Gene Name`, `Distance to TSS`)

joined_data <- right_join(annotation_data, filtered_data, by = c("Gene Name" = "genename"))

# status based on distance to TSS
calculate_status <- function(data, distance_threshold, regulated_name) {
  data %>%
    mutate(
      status = ifelse(abs(`Distance to TSS`) < distance_threshold & !is.na(`Distance to TSS`), "Bound", "Unbound"),
      regulated = ifelse(`log2FoldChange` < 0, paste0(regulated_name, " Down-regulated"), paste0(regulated_name, " Up-regulated"))
    )
}

# calc status for different distances
data5kb <- calculate_status(joined_data, 5000, "5kb")
data20kb <- calculate_status(joined_data, 20000, "20kb")

# merge data and make bound col
combined_data <- bind_rows(data5kb, data20kb)
combined_data$status <- factor(combined_data$status, levels = c("Unbound", "Bound"))

my_colors <- c("Bound" = "#FF3D00", "Unbound" = "#8B8878")

# create ggplot obj
p <- ggplot(combined_data, aes(x = regulated, fill = status)) +
  geom_bar(position = "fill") +
  geom_text(
    aes(label = after_stat(count)), 
    stat = "count", 
    position = position_fill(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  scale_fill_manual(values = my_colors, labels = c("Not Bound", "RUNX1 Bound")) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_x_discrete(limits = c("5kb Up-regulated", "5kb Down-regulated", "20kb Up-regulated", "20kb Down-regulated")) +
  labs(x = "5kb of TSS        20kb of TSS", y = "Percentage of Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save plot w/ lab
ggsave("plot5-20.png", p, width = 10, height = 6, units = "in")
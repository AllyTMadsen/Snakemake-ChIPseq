library(tidyverse)
library(patchwork)

deseq_data <- read.table(gzfile("results/RNAseq_log2_foldchange.txt.gz"), header = TRUE)

filtered_data <- deseq_data %>%
  filter(padj < 0.01, abs(log2FoldChange) > 1)

annotations <- "results/annotated_peaks.txt"
annotation_data <- read_delim(file = annotations, delim = "\t") %>%
  as_tibble() %>%
  select(`Gene Name`, `Distance to TSS`)

joined_data <- right_join(annotation_data, filtered_data, by = c("Gene Name" = "genename"))

calculate_status <- function(data, distance_threshold, regulated_name) {
  data %>%
    mutate(
      status = ifelse(abs(`Distance to TSS`) < distance_threshold & !is.na(`Distance to TSS`), "Bound", "Unbound"),
      regulated = ifelse(`log2FoldChange` < 0, paste0(regulated_name, " Down-regulated"), paste0(regulated_name, " Up-regulated"))
    )
}

data100kb <- calculate_status(joined_data, 100000, "100kb")

my_colors <- c("Bound" = "#8B8878", "Unbound" = "#FF3D00")

p <- ggplot(data100kb, aes(x = regulated, fill = status)) +
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
  scale_x_discrete(limits = c("100kb Up-regulated", "100kb Down-regulated")) +
  labs(x = "100kb of TSS", y = "Percentage of Genes")

ggsave("plot100.png", p, width = 10, height = 6, units = "in")
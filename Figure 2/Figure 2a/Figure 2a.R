#Figure 2a
library(ggplot2)
setwd("D:/PMI_Project/Figure 2/Figure 2a")

data_file <- "metabolite_matrix.csv"
metadata_file <- "Group.csv"
output_file <- "Figure2a_PCA.tiff"

data <- read.csv(
  "metabolite_matrix.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  "Group.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

common_samples <- intersect(rownames(data), rownames(metadata))
stopifnot(length(common_samples) > 0)

data <- data[common_samples, , drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

stopifnot(all(c("Day2", "Bodysite") %in% colnames(metadata)))

data_centered <- scale(data, center = TRUE, scale = FALSE)

pareto_scale_factors <- apply(data, 2, function(x) sqrt(sd(x)))
pareto_scale_factors[is.na(pareto_scale_factors) | pareto_scale_factors == 0] <- 1

data_pareto <- sweep(data_centered, 2, pareto_scale_factors, FUN = "/")

pca_result <- prcomp(data_pareto, center = FALSE, scale. = FALSE)

axes <- as.data.frame(pca_result$x[, 1:2, drop = FALSE])
colnames(axes) <- c("Axis.1", "Axis.2")
axes <- cbind(axes, metadata[rownames(axes), , drop = FALSE])

eigval <- round(summary(pca_result)$importance[2, ] * 100, 2)
day2_levels <- unique(axes$Day2)
bodysite_levels <- unique(axes$Bodysite)

custom_day2_colors <- c(
  "D3" = "#FF9FF3",
  "B1" = "#A3DEF4",
  "C2" = "#2A9D8F",
  "A0" = "#F4A261",
  "E5" = "#9C27B0",
  "F7" = "#FF6B6B"
)

stopifnot(all(day2_levels %in% names(custom_day2_colors)))

shape_values <- c(8, 16, 17, 15, 18)
stopifnot(length(bodysite_levels) <= length(shape_values))
bodysite_shapes <- setNames(shape_values[seq_along(bodysite_levels)], bodysite_levels)

axes$Day2 <- factor(axes$Day2, levels = day2_levels)
axes$Bodysite <- factor(axes$Bodysite, levels = bodysite_levels)

p <- ggplot(axes, aes(x = Axis.1, y = Axis.2, colour = Day2, shape = Bodysite)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(group = Day2), level = 0.95, linetype = 1) +
  scale_colour_manual(
    name = "PMI",
    values = custom_day2_colors,
    breaks = day2_levels
  ) +
  scale_shape_manual(
    name = "Bodysite",
    values = bodysite_shapes,
    breaks = bodysite_levels
  ) +
  labs(
    x = paste0("PC1 (", eigval[1], " %)"),
    y = paste0("PC2 (", eigval[2], " %)")
  ) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray50") +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  theme_bw(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave(
  filename = output_file,
  plot = p,
  device = "tiff",
  width = 10,
  height = 8,
  dpi = 300,
  compression = "lzw"
)

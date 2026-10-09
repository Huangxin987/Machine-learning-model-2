# Figure 3C: Bacteria PCoA
library(ape)
library(RColorBrewer)
library(ggplot2)
library(vegan)

setwd("")
output_dir <- file.path(getwd(), "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

data <- read.csv("ASV.csv", row.names = 1, stringsAsFactors = FALSE)
data <- t(data)

group <- read.csv("group.csv", row.names = 1, stringsAsFactors = FALSE)

common_samples <- intersect(rownames(data), rownames(group))
data <- data[common_samples, ]
group <- group[common_samples, ]

dis <- vegdist(data, "bray")
pcoa <- pcoa(dis)
axes <- as.data.frame(pcoa$vectors)
axes <- cbind(axes, group[rownames(axes), ])

eigval <- round(pcoa$values$Relative_eig * 100, digits = 2)

custom_day2_colors <- c(
  "D7" = "#E63946",
  "B3" = "#2A9D8F",
  "A2" = "#F4A261",
  "C5" = "#9C27B0"
)

bodysite_shapes <- c(8, 16, 17, 15, 18)
names(bodysite_shapes) <- unique(axes$Bodysite)

p <- ggplot(axes, aes(x = Axis.1, y = Axis.2, colour = Day2, shape = Bodysite)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_colour_manual(
    name = "PMI",
    values = custom_day2_colors,
    breaks = unique(axes$Day2)
  ) +
  scale_shape_manual(
    name = "Bodysite",
    values = bodysite_shapes,
    breaks = unique(axes$Bodysite)
  ) +
  xlab(paste0("PCo1 (", eigval[1], " %)")) +
  ylab(paste0("PCo2 (", eigval[2], " %)")) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray50") +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  stat_ellipse(aes(group = Day2), level = 0.95, linetype = 1) +
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

print(p)
ggsave("Figure3C_bacteria_PCoA.png", p, width = 8, height = 6, dpi = 300)
